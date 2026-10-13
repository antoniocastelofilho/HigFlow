#include "navierstokes2d.hxx"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>

namespace t8flow2d {

namespace {

constexpr double EPS = 1e-12;

//! Quantiza um ponto em meios-passos.  Exato para as posicoes de faceta e de
//! centro de uma malha uniforme, que e' a unica que este solver aceita.
unsigned long long chave_de_ponto(const Ponto2 &x, const Ponto2 &h)
{
  const long long qx = std::llround(x[0] / (h[0] / 2.0));
  const long long qy = std::llround(x[1] / (h[1] / 2.0));
  // DESLOCAMENTO para que negativos nao colidam com positivos ao truncar.
  const unsigned long long ux = static_cast<unsigned long long>(qx + (1LL << 30));
  const unsigned long long uy = static_cast<unsigned long long>(qy + (1LL << 30));
  return (ux << 32) | (uy & 0xFFFFFFFFULL);
}

}  // namespace

std::string nome_de(Tempo t)
{
  switch (t) {
    case Tempo::EulerExplicito: return "euler_explicito";
    case Tempo::EulerSemiImplicito: return "euler_semi_implicito";
  }
  return "?";
}

std::string nome_de(Projecao p)
{
  return p == Projecao::Incremental ? "incremental" : "nao_incremental";
}

std::string nome_de(Convectivo c)
{
  switch (c) {
    case Convectivo::Central: return "central";
    case Convectivo::Upwind1: return "upwind1";
    case Convectivo::Cubista: return "cubista";
    case Convectivo::Quick: return "quick";
  }
  return "?";
}

Tempo tempo_do_nome(const std::string &s)
{
  if (s == "euler_explicito" || s == "explicit_euler") return Tempo::EulerExplicito;
  if (s == "euler_semi_implicito" || s == "semi_implicit_euler")
    return Tempo::EulerSemiImplicito;
  throw std::invalid_argument(
      "tempo_do_nome: \"" + s + "\" nao e' valido.  Validos: euler_explicito, "
      "euler_semi_implicito (e os nomes do HiGFlow explicit_euler, "
      "semi_implicit_euler)");
}

Projecao projecao_do_nome(const std::string &s)
{
  if (s == "incremental") return Projecao::Incremental;
  if (s == "nao_incremental" || s == "non_incremental") return Projecao::NaoIncremental;
  throw std::invalid_argument("projecao_do_nome: \"" + s +
                              "\" nao e' valido.  Validos: incremental, "
                              "nao_incremental (HiGFlow: non_incremental)");
}

Convectivo convectivo_do_nome(const std::string &s)
{
  if (s == "central") return Convectivo::Central;
  if (s == "upwind1" || s == "first_order") return Convectivo::Upwind1;
  if (s == "cubista") return Convectivo::Cubista;
  if (s == "quick") return Convectivo::Quick;
  throw std::invalid_argument("convectivo_do_nome: \"" + s +
                              "\" nao e' valido.  Validos: central, upwind1, "
                              "cubista, quick");
}

NavierStokes::NavierStokes(const Malha &m, const DominioFacetas &fd,
                           const ContornoPorRegiao &bc_u0,
                           const ContornoPorRegiao &bc_u1,
                           const ContornoPorRegiao &bc_p)
    : m_(&m), fd_(&fd), u_(fd), ustar_(fd), p_(m), dp_(m)
{
  bcu_[0] = &bc_u0;
  bcu_[1] = &bc_u1;
  bcp_ = &bc_p;

  // O PARALELO E' RECUSADO, e isto e' limitacao declarada e nao defeito
  // escondido.  Faltam DUAS pecas, e nenhuma e' pequena:
  //
  //   A TROCA DE FRANJA DOS CAMPOS.  `CampoEscalar` cobre os ids LOCAIS, e o
  //   estencil da equacao le' celula e faceta do vizinho.  A malha ja' tem
  //   franja; os campos ainda nao a usam;
  //
  //   O PERCURSO POR BLOCOS.  Os lacos aqui vao de `0` a `num_blocos()`, e em
  //   paralelo um processo nao possui todas as arvores -- o que aparece como
  //   "arvore 2 nao e' local", longe da causa.
  //
  // Recusar e' melhor que rodar: com a franja ausente o estencil leria lixo
  // perto da fronteira de particao e a corrida daria numero.
  if (m.num_processos() > 1)
    throw std::invalid_argument(
        "NavierStokes: este solver ainda e' SERIAL (" +
        std::to_string(m.num_processos()) +
        " processos pedidos).  Faltam a troca de franja dos campos -- o "
        "estencil le' celula e faceta do vizinho, e `CampoEscalar` cobre so' "
        "os ids locais -- e o percurso por blocos, que aqui supoe todas as "
        "arvores locais.");

  // A MALHA ADAPTADA, NESTA ALTURA, SO' E' RECUSADA NO PASSO -- ver `passo`.
  //
  // O que ja' a trata sao os operadores de FLUXO: a divergencia e o Poisson da
  // projecao somam sub-facetas e usam a formula unica de face e distancia.  O
  // que ainda nao a trata e' o ESTENCIL DA QUANTIDADE DE MOVIMENTO, cuja
  // restricao e' diferente por componente.  Construir e medir o fluxo numa
  // malha adaptada tem de ser possivel; avancar no tempo ainda nao.
  //
  // O TEXTO ABAIXO FICA, porque e' a razao de o passo recusar:
  //
  // A transicao 2:1 numa faceta de VELOCIDADE nao e' o mesmo problema que numa
  // faceta de fluxo escalar: a componente normal se restringe por conservacao
  // de vazao e a tangencial por interpolacao, e usar a mesma regra para as duas
  // produz divergencia que nao fecha -- um defeito que aparece como o solver de
  // pressao nao convergindo, longe da causa.
  //
  // Rodar e dar numero seria pior que recusar: o numero sairia plausivel.

  if (bc_p.vazio())
    throw std::invalid_argument(
        "NavierStokes: o contorno da pressao esta' vazio.  Sem condicao em "
        "parede nenhuma o Poisson da projecao e' singular.");

  // O LADO DA FOLHA, e a conferencia de que ele e' o MESMO em todas.
  //
  // A chave do indice de facetas quantiza em `h/2`, o que so' e' exato se a
  // malha for uniforme.  Com folhas de tamanhos diferentes a quantizacao
  // colidiria e duas facetas distintas receberiam o mesmo indice -- e isso
  // nao falharia, produziria o valor da faceta errada.
  //
  // A transicao 2:1 ja' foi recusada acima, mas blocos com NIVEIS diferentes
  // passariam por ela: sao uniformes dentro de cada bloco e diferentes entre
  // si.  Por isso a conferencia e' sobre todas as folhas, e nao sobre o
  // contador de transicao.
  if (m.num_folhas() == 0)
    throw std::invalid_argument("NavierStokes: a malha nao tem folha local");
  h_ = m.lados_da_folha(0, 0);
  h_min_ = h_;
  uniforme_ = true;
  for (long b = 0; b < m.num_blocos(); b++)
    for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
      const Ponto2 hh = m.lados_da_folha(b, i);
      for (int d = 0; d < 2; d++) {
        if (hh[d] < h_min_[d]) h_min_[d] = hh[d];
        if (std::fabs(hh[d] - h_[d]) > 1e-12 * h_[d]) uniforme_ = false;
      }
    }

  // O INDICE DE FACETA POR PONTO.  Construido uma vez; a alternativa seria
  // varrer as facetas a cada consulta, e a montagem do passo consulta uma vez
  // por termo por faceta -- o mesmo quadratico que ja' custou caro no Poisson.
  for (int d = 0; d < 2; d++) {
    fac_por_ponto_[d].reserve(static_cast<size_t>(fd.num_facetas(d)) * 2);
    for (long f = 0; f < fd.num_facetas(d); f++) {
      const unsigned long long ch = chave_de_ponto(fd.centro(d, f), h_min_);
      const auto it = fac_por_ponto_[d].find(ch);
      if (it != fac_por_ponto_[d].end())
        throw std::runtime_error(
            "NavierStokes: duas facetas da direcao " + std::to_string(d) +
            " caem na mesma chave de indice -- a quantizacao nao serve para "
            "esta malha");
      fac_por_ponto_[d][ch] = f;
    }
  }

  // O POISSON DA PROJECAO NAO PODE SER TODO-NEUMANN.
  //
  // Com velocidade de Dirichlet em toda a volta -- a caixa com tampa
  // deslizante e' o caso classico --, a condicao natural da pressao e' Neumann
  // em todas as paredes, e ai' o operador tem a CONSTANTE no nucleo: a pressao
  // fica definida a menos de constante e a matriz e' singular.
  //
  // A SAIDA NAO E' FIXAR UMA CELULA INTERIOR por conta propria, e sim exigir
  // que o caso imponha Dirichlet em pelo menos um ponto de parede -- ali o
  // valor tem significado fisico e, numa solucao manufaturada, e' conhecido.
  // Fixar uma celula escolhida pelo solver poria a arbitrariedade dentro da
  // biblioteca, onde quem le' o resultado nao a ve'.
  //
  // E recusar e' melhor que resolver: o PETSc com matriz singular ou nao
  // converge ou converge para um elemento qualquer do subespaco, e o segundo
  // caso da' numero.
  p_tem_dirichlet_ = false;
  for (long b = 0; b < m.num_blocos() && !p_tem_dirichlet_; b++)
    for (long i = 0; i < m.num_folhas_do_bloco(b) && !p_tem_dirichlet_; i++) {
      const Ponto2 c = m.centroide(b, i);
      const Ponto2 hh = m.lados_da_folha(b, i);
      for (int d = 0; d < 2 && !p_tem_dirichlet_; d++)
        for (int s2 = -1; s2 <= 1; s2 += 2) {
          Ponto2 xn = c, xp = c;
          xn[d] += s2 * hh[d];
          xp[d] += s2 * hh[d] / 2.0;
          if (m.id_do_ponto(xn) >= 0) continue;      // face interior
          if (bc_p.em(xp).tipo == TipoContorno::Dirichlet) {
            p_tem_dirichlet_ = true;
            break;
          }
        }
    }
  if (!p_tem_dirichlet_)
    throw std::invalid_argument(
        "NavierStokes: a pressao tem Neumann em TODA a parede, e o Poisson da "
        "projecao fica singular -- a constante esta' no nucleo.  Imponha "
        "Dirichlet em pelo menos um ponto de parede; numa solucao manufaturada "
        "o valor ali e' conhecido.");

  // `(bloco, folha)` POR ID, e o DONO de cada faceta.  Os dois sao consultados
  // uma vez por termo por faceta; monta-los por varredura a cada consulta
  // daria o mesmo quadratico que ja' custou caro no Poisson.
  bi_por_id_.reserve(static_cast<size_t>(m.num_folhas()) * 2);
  for (long b = 0; b < m.num_blocos(); b++)
    for (long i = 0; i < m.num_folhas_do_bloco(b); i++)
      bi_por_id_[m.id_global(b, i)] = {b, i};

  for (int d = 0; d < 2; d++) dono_[d].assign(
      static_cast<size_t>(fd.num_facetas(d)), Dono{-1, -1, -1});
  for (long b = 0; b < m.num_blocos(); b++)
    for (long i = 0; i < m.num_folhas_do_bloco(b); i++)
      for (int d = 0; d < 2; d++)
        for (int lado = 0; lado < 2; lado++) {
          const long f = fd.id_faceta(b, i, d, lado);
          if (f >= 0) dono_[d][static_cast<size_t>(f)] = Dono{b, i, lado};
        }
  for (int d = 0; d < 2; d++)
    for (long f = 0; f < fd.num_facetas(d); f++)
      if (dono_[d][static_cast<size_t>(f)].b < 0)
        throw std::runtime_error(
            "NavierStokes: a faceta " + std::to_string(f) + " da direcao " +
            std::to_string(d) + " nao tem dono.  O estencil se orienta pela "
            "celula dona -- o passo ao longo da normal e' o lado DELA, e sem "
            "dono nao ha' passo.");

  u_.preenche(0.0);
  ustar_.preenche(0.0);
  p_.preenche(0.0);
  dp_.preenche(0.0);
}

void NavierStokes::reynolds(double re)
{
  if (!(re > 0.0))
    throw std::invalid_argument("NavierStokes::reynolds: " + std::to_string(re) +
                                " -- o numero de Reynolds e' positivo");
  re_ = re;
}

void NavierStokes::passo_de_tempo(double dt)
{
  if (!(dt > 0.0))
    throw std::invalid_argument("NavierStokes::passo_de_tempo: " +
                                std::to_string(dt) + " -- o passo e' positivo");
  dt_ = dt;
}

// ===========================================================================
// AS BUSCAS QUE O ESQUEMA PRECISA
// ===========================================================================

long NavierStokes::faceta_em(int d, const Ponto2 &x) const
{
  const auto it = fac_por_ponto_[d].find(chave_de_ponto(x, h_min_));
  return (it == fac_por_ponto_[d].end()) ? -1 : it->second;
}

// ===========================================================================
// OS TERMOS DA EQUACAO
// ===========================================================================

namespace {

//! O limitador em VARIAVEL NORMALIZADA (NVD).
//
// `fc` e' `(phi_C - phi_U) / (phi_D - phi_U)`, e devolve `phi_f` normalizado.
// FORA DE [0,1] TODOS CAEM PARA UPWIND, e isso nao e' detalhe: e' o que torna
// o esquema limitado.  Um "segunda ordem" que nao cai e' central com outro
// nome, e oscila onde o gradiente e' forte.
double nvd(Convectivo c, double fc)
{
  if (fc <= 0.0 || fc >= 1.0) return fc;        // upwind: limitado
  switch (c) {
    case Convectivo::Upwind1:
      return fc;
    case Convectivo::Central:
      return 0.5 * (fc + 1.0);                  // interpolacao linear, em NVD
    case Convectivo::Cubista:
      // Alves, Oliveira e Pinho (2003).  As tres faixas casam nos pontos de
      // quebra -- em 3/8 as duas primeiras dao 21/32, em 3/4 as duas ultimas
      // dao 15/16 --, e e' essa continuidade que o torna utilizavel.
      if (fc < 3.0 / 8.0) return (7.0 / 4.0) * fc;
      if (fc <= 3.0 / 4.0) return (3.0 / 4.0) * fc + 3.0 / 8.0;
      return (1.0 / 4.0) * fc + 3.0 / 4.0;
    case Convectivo::Quick:
      return (3.0 / 4.0) * fc + 3.0 / 8.0;
  }
  return fc;
}

}  // namespace

double NavierStokes::contorno_vizinho(int d, long fid, int e, int s) const
{
  // O FECHAMENTO DE CONTORNO, com a parede a MEIO LADO DA CELULA DONA.
  //
  // `h/2` da celula dona, e nao um `h` global: numa malha adaptada a faceta
  // junto a' parede pode ser fina ou grossa, e refletir pela distancia errada
  // poe a parede no lugar errado -- sem erro visivel, com a condicao valendo
  // em outro ponto.
  const Dono &dn = dono_[d][static_cast<size_t>(fid)];
  const Ponto2 hC = m_->lados_da_folha(dn.b, dn.i);
  const Ponto2 xf = fd_->centro(d, fid);
  const double aqui = u_.em(d, fid);
  if (e == d) return aqui;              // nao ocorre em faceta interior

  Ponto2 xp = xf;
  xp[e] += s * hC[e] / 2.0;
  const Contorno &c = bcu_[d]->em(xp);
  switch (c.tipo) {
    case TipoContorno::Dirichlet:
      return 2.0 * c.g(xp) - aqui;      // reflexao: a media da' o valor na parede
    case TipoContorno::Neumann:
      return aqui + s * hC[e] * c.g(xp);
    case TipoContorno::Robin: {
      double A = 0.0, B = 0.0;
      c.afim(xp, hC[e], &A, &B);
      return A * aqui + B;
    }
  }
  return aqui;
}

double NavierStokes::dist_ate_parede(int d, long fid, int e, int s) const
{
  (void) s;
  const Dono &dn = dono_[d][static_cast<size_t>(fid)];
  const Ponto2 hC = m_->lados_da_folha(dn.b, dn.i);
  // O FANTASMA esta' a' distancia de um lado inteiro: meio ate' a parede, meio
  // do outro lado dela.
  return hC[e];
}

double NavierStokes::comprimento_do_volume(int d, long fid, int e) const
{
  // O VOLUME DE CONTROLE DA FACETA: meia celula de cada lado ao longo da
  // normal, e a celula dona inteira na transversal.
  const Dono &dn = dono_[d][static_cast<size_t>(fid)];
  const Ponto2 hC = m_->lados_da_folha(dn.b, dn.i);
  if (e != d) return hC[e];

  const Ponto2 xf = fd_->centro(d, fid);
  double soma = 0.0;
  int n = 0;
  for (int lado = 0; lado < 2; lado++) {
    Ponto2 p = xf;
    p[d] += (lado == 0 ? -1.0 : 1.0) * h_min_[d] / 4.0;
    const long id = m_->id_do_ponto(p);
    const auto it = (id >= 0) ? bi_por_id_.find(id) : bi_por_id_.end();
    if (it == bi_por_id_.end()) continue;
    soma += m_->lados_da_folha(it->second.first, it->second.second)[d] / 2.0;
    n++;
  }
  return (n == 2) ? soma : hC[d];
}

double NavierStokes::velocidade_que_advecta(int d, long fid, int e, int s) const
{
  // A VELOCIDADE NA POSICAO DO FLUXO.  Para `e == d` essa posicao e' um CENTRO
  // DE CELULA e a componente que atravessa e' a propria `d`; para `e != d` e'
  // uma QUINA, e a que atravessa e' a componente `e`.
  //
  // Confundir as duas e' o erro classico do arranjo deslocado, e ele nao
  // falha: produz um termo convectivo plausivel e errado.
  const Dono &dn = dono_[d][static_cast<size_t>(fid)];
  const Ponto2 hC = m_->lados_da_folha(dn.b, dn.i);
  const Ponto2 cC = m_->centroide(dn.b, dn.i);
  const Ponto2 xf = fd_->centro(d, fid);
  bool ok = false;

  if (e == d) {
    Ponto2 p = xf;
    p[d] += s * h_min_[d] / 4.0;
    const long id = m_->id_do_ponto(p);
    const auto it = (id >= 0) ? bi_por_id_.find(id) : bi_por_id_.end();
    if (it == bi_por_id_.end()) return u_.em(d, fid);
    const double v = u_no_ponto(
        d, m_->centroide(it->second.first, it->second.second), &ok);
    return ok ? v : u_.em(d, fid);
  }

  Ponto2 quina = xf;
  quina[e] = cC[e] + s * hC[e] / 2.0;
  const double v = u_no_ponto(e, quina, &ok);
  return ok ? v : 0.0;
}

double NavierStokes::valor_na_face(int d, long fid, int e, int s,
                                   double V) const
{
  // QUEM E' `C`, `D` E `U` DEPENDE DO SINAL DE `V`, e nao da geometria: o
  // esquema olha para montante, e montante troca de lado quando o escoamento
  // troca.  Fixar `C` como "o daqui" faz o upwind virar downwind na metade do
  // dominio, e isso nao falha -- difunde.
  const double aqui = u_.em(d, fid);
  const Vizinha vz = vizinha_de(d, fid, e, s);
  const double alem = vz.dentro ? vz.valor : contorno_vizinho(d, fid, e, s);

  if (conv_ == Convectivo::Upwind1) return (V >= 0.0) ? aqui : alem;
  if (conv_ == Convectivo::Central) return 0.5 * (aqui + alem);

  double phiC, phiD, phiU;
  bool tem_u = false;
  if (V >= 0.0) {
    phiC = aqui;
    phiD = alem;
    // O PONTO A MONTANTE TEM DE EXISTIR DE VERDADE.  Fora do dominio ha' o
    // fantasma do contorno, e usa-lo aqui seria legitimo -- mas e' OUTRO
    // esquema, e trocar o esquema junto com o suporte 2:1 faria qualquer
    // diferenca medida ficar sem dono.  Fora do dominio, rebaixa e CONTA.
    const Vizinha vu = vizinha_de(d, fid, e, -s);
    tem_u = vu.dentro;
    phiU = vu.dentro ? vu.valor : phiC;
  } else {
    phiC = alem;
    phiD = aqui;
    // O PONTO A MONTANTE ESTA' DOIS PASSOS ADIANTE, e ali o fechamento de
    // contorno NAO vale: ele supoe a parede a meio caminho, e do segundo
    // passo em diante ela ja' ficou para tras.  A busca e' ESTRITA.
    Ponto2 x2 = fd_->centro(d, fid);
    x2[e] += 2.0 * s * m_->lados_da_folha(dono_[d][static_cast<size_t>(fid)].b,
                                          dono_[d][static_cast<size_t>(fid)].i)[e];
    const long f2 = faceta_em(d, x2);
    if (f2 >= 0) { phiU = u_.em(d, f2); tem_u = true; }
    else { phiU = phiC; tem_u = false; }
  }

  // SEM O PONTO A MONTANTE o esquema de segunda ordem nao existe, e a queda
  // para upwind e' CONTADA.  Perto da parede e' inevitavel; no meio do
  // dominio seria sintoma.
  if (!tem_u) { rebaixos_conv_++; return phiC; }

  const double den = phiD - phiU;
  if (std::fabs(den) < EPS) return phiC;          // campo localmente plano
  const double fc = (phiC - phiU) / den;
  return phiU + nvd(conv_, fc) * den;
}

double NavierStokes::conveccao(int d, long fid) const
{
  // FORMA DIVERGENTE, `sum_e d(u_e u_d)/dx_e`, e nao a advectiva.
  //
  // As duas coincidem no continuo porque a divergencia e' nula, e NAO
  // coincidem no discreto: a divergente conserva quantidade de movimento por
  // construcao, e a advectiva so' a conserva quando a divergencia discreta ja'
  // e' nula -- o que no meio do passo ela nao e'.
  double soma = 0.0;
  for (int e = 0; e < 2; e++) {
    const double Vp = velocidade_que_advecta(d, fid, e, +1);
    const double Vm = velocidade_que_advecta(d, fid, e, -1);
    const double fp = valor_na_face(d, fid, e, +1, Vp);
    const double fm = valor_na_face(d, fid, e, -1, Vm);
    const double L = comprimento_do_volume(d, fid, e);
    if (L > 0.0) soma += (Vp * fp - Vm * fm) / L;
  }
  return soma;
}

double NavierStokes::difusao(int d, long fid) const
{
  // A SEGUNDA DERIVADA COM ESPACAMENTO DESIGUAL, e nao a de tres pontos
  // igualmente espacados:
  //
  //     d2u/dx2 = 2/(dm+dp) . [ (up-u0)/dp - (u0-um)/dm ]
  //
  // Em malha uniforme `dm = dp` e ela vira a formula de sempre -- um caminho
  // so', e e' por isso que os casos uniformes verificam este codigo.
  //
  // Na transicao a formula e' de PRIMEIRA ordem, e nao de segunda: com
  // `dm != dp` o termo de terceira derivada nao cancela.  Isso e' propriedade
  // do estencil desigual, nao defeito -- e o preco aparece na norma do maximo
  // perto da interface, nao na norma 2.  Fica dito para nao ser atribuido ao
  // esquema depois.
  const double aqui = u_.em(d, fid);
  double lap = 0.0;
  for (int e = 0; e < 2; e++) {
    const Vizinha vp = vizinha_de(d, fid, e, +1);
    const Vizinha vm = vizinha_de(d, fid, e, -1);
    const double up = vp.dentro ? vp.valor : contorno_vizinho(d, fid, e, +1);
    const double um = vm.dentro ? vm.valor : contorno_vizinho(d, fid, e, -1);
    const double dp = vp.dentro ? vp.dist : dist_ate_parede(d, fid, e, +1);
    const double dm = vm.dentro ? vm.dist : dist_ate_parede(d, fid, e, -1);
    if (dp <= 0.0 || dm <= 0.0) continue;
    lap += 2.0 / (dm + dp) * ((up - aqui) / dp - (aqui - um) / dm);
  }
  return lap / re_;
}

double NavierStokes::grad_campo(const CampoEscalar &campo, int d,
                                long fid) const
{
  // A PRESSAO VIVE NOS CENTROS, e a faceta esta' entre dois deles -- mas em
  // malha adaptada NAO esta' no meio, e os centros NAO estao alinhados com a
  // normal.  Sao dois consertos, e nenhum e' opcional.
  //
  // A DISTANCIA e' a soma dos dois meios-lados, `(h_menos + h_mais)/2`, e nao
  // `h`.  Usar `h` da celula fina daria um gradiente duas vezes maior do lado
  // grosso.
  //
  // A CORRECAO TRANSVERSAL desloca o valor do lado GROSSO para a linha
  // transversal da faceta: o centro da grossa esta' deslocado em `t`, e a
  // diferenca bruta mede a derivada ao longo da LINHA DOS CENTROS, que tem
  // componente transversal.  Para a celula FINA esse erro NAO cancela, porque
  // ela ve' uma unica face contra a grossa.
  //
  // MEDIDO NO POISSON deste mesmo projeto, sem a correcao: o grau 1, exato em
  // malha uniforme, deu 3,9e-02 na adaptada, e a solucao trigonometrica teve
  // ordem NEGATIVA.  Nao era erro de O(h): era de O(1).
  const Ponto2 xf = fd_->centro(d, fid);
  const int t = 1 - d;

  long bi[2] = {-1, -1}, ii[2] = {-1, -1};
  for (int lado = 0; lado < 2; lado++) {
    Ponto2 p = xf;
    p[d] += (lado == 0 ? -1.0 : 1.0) * h_min_[d] / 4.0;
    const long id = m_->id_do_ponto(p);
    const auto it = (id >= 0) ? bi_por_id_.find(id) : bi_por_id_.end();
    if (it == bi_por_id_.end()) return 0.0;   // faceta de contorno: nao se resolve ali
    bi[lado] = it->second.first;
    ii[lado] = it->second.second;
  }

  double val[2], meio[2];
  for (int lado = 0; lado < 2; lado++) {
    const Ponto2 c = m_->centroide(bi[lado], ii[lado]);
    const Ponto2 hh = m_->lados_da_folha(bi[lado], ii[lado]);
    meio[lado] = hh[d] / 2.0;
    const double delta = xf[t] - c[t];        // zero quando os centros alinham
    val[lado] = campo.em(m_->id_global(bi[lado], ii[lado]));
    if (std::fabs(delta) > 1e-14 * hh[t])
      val[lado] += delta * grad_transversal(campo, bi[lado], ii[lado], t);
  }
  return (val[1] - val[0]) / (meio[0] + meio[1]);
}

// ===========================================================================
// O PASSO
// ===========================================================================

void NavierStokes::aplica_contorno_velocidade()
{
  // UMA FACETA E' DE CONTORNO quando uma das duas celulas que ela separa nao
  // existe.  Ali a velocidade nao se resolve: ela e' DADA (Dirichlet) ou
  // extrapolada (Neumann).
  for (int d = 0; d < 2; d++)
    for (long f = 0; f < fd_->num_facetas(d); f++) {
      const Ponto2 x = fd_->centro(d, f);
      const Dono &dn = dono_[d][static_cast<size_t>(f)];
      const Ponto2 hC = m_->lados_da_folha(dn.b, dn.i);
      Ponto2 a = x, b = x;
      a[d] -= h_min_[d] / 4.0;
      b[d] += h_min_[d] / 4.0;
      const long ca = m_->id_do_ponto(a), cb = m_->id_do_ponto(b);
      if (ca >= 0 && cb >= 0) continue;        // interior
      if (ca < 0 && cb < 0)
        throw std::runtime_error(
            "NavierStokes: a faceta " + std::to_string(f) + " da direcao " +
            std::to_string(d) + " nao tem celula de nenhum lado");

      const Contorno &c = bcu_[d]->em(x);
      if (c.tipo == TipoContorno::Dirichlet) {
        u_.em(d, f) = c.g(x);
      } else {
        // NEUMANN NA COMPONENTE NORMAL: `du/dn = g`, e o unico valor interior
        // disponivel esta' uma faceta adentro.
        const int s = (ca >= 0) ? -1 : +1;      // para dentro
        // O PASSO E' O LADO DA CELULA DONA, nao um `h` global: a faceta junto
        // a' parede pode ser fina ou grossa, e extrapolar pela distancia
        // errada poe a condicao valendo em outro ponto, sem erro visivel.
        Ponto2 xi = x;
        xi[d] += s * hC[d];
        const long fi = faceta_em(d, xi);
        if (fi < 0) continue;                   // dominio de uma celula so'
        u_.em(d, f) = u_.em(d, fi) - s * hC[d] * c.g(x);
      }
    }
}

void NavierStokes::preditor()
{
  // O PREDITOR NAO TOCA AS FACETAS DE CONTORNO: ali a velocidade vem da
  // condicao, e sobrescreve-la com a equacao e' como a condicao de contorno
  // deixa de valer sem que nada acuse.
  ustar_ = u_;
  rebaixos_conv_ = 0;

  for (int d = 0; d < 2; d++)
    for (long f = 0; f < fd_->num_facetas(d); f++) {
      const Ponto2 x = fd_->centro(d, f);
      Ponto2 a = x, b = x;
      a[d] -= h_min_[d] / 4.0;
      b[d] += h_min_[d] / 4.0;
      if (m_->id_do_ponto(a) < 0 || m_->id_do_ponto(b) < 0) continue;

      // `RHO` NA FACETA, e nao na celula: a equacao de quantidade de
      // movimento vive na faceta, e dividir pela densidade da celula errada
      // desloca a inercia meio passo -- o que, atraves de uma interface com
      // razao 10 de densidade, e' erro de dezenas por cento.
      const double rho_f = rho_na_faceta(d, f);

      // A DIVISAO POR `RHO` E' SO' DOS TERMOS DE FORCA, e nao do convectivo:
      // `div(uu)` ja' e' aceleracao.  Dividi-lo tambem seria dividir duas
      // vezes quando `rho` e' constante e uma vez de menos quando nao e'.
      double rhs = -conveccao(d, f) + g_[d]
                   + (tensao_viscosa(d, f) + (forca_campo_ ? forca_campo_->em(d, f) : 0.0))
                         / rho_f;

      // A PRESSAO ANTIGA ENTRA SO' NA INCREMENTAL, e e' nisso que as duas
      // diferem.  Sem ela o preditor ignora o gradiente que ja' se conhece, e
      // a projecao tem de reconstrui-lo inteiro a cada passo -- o que custa
      // uma ordem na pressao.
      if (proj_ == Projecao::Incremental) rhs -= grad_campo(p_, d, f) / rho_f;

      if (fonte_faceta_) rhs += fonte_faceta_(x, d, t_bc_);

      ustar_.em(d, f) = u_.em(d, f) + dt_ * rhs;
    }
}

void NavierStokes::resolve_pressao()
{
  // O POISSON DA PROJECAO: `lap(dp) = div(u*) / dt`.
  //
  // A CONDICAO DELE NAO E' A DA PRESSAO, e' DERIVADA dela -- confundir as duas
  // e' erro de formulacao, nao de implementacao.  Onde a pressao tem NEUMANN o
  // incremento tem Neumann HOMOGENEO; onde tem DIRICHLET, o incremento vale
  // `g` na nao incremental (ali `p = dp`) e `g - p^n_parede` na incremental --
  // Dirichlet homogeneo prenderia a VARIACAO e nao a pressao, e o nivel ficaria
  // no que o transiente deixou.
  const long n = m_->num_folhas();
  const long p0 = m_->primeiro_id();
  t8flow::Solver s(t8flow::Motor::Petsc, static_cast<int>(p0),
                   static_cast<std::size_t>(n), 9);
  double soma_b = 0.0;
  s.tolerancia(1e-12);

  for (long b = 0; b < m_->num_blocos(); b++)
    for (long i = 0; i < m_->num_folhas_do_bloco(b); i++) {
      const long id = m_->id_global(b, i);
      const Ponto2 c = m_->centroide(b, i);
      const Ponto2 hh = m_->lados_da_folha(b, i);
      const double vol = hh[0] * hh[1];

      std::map<long, double> linha;
      double rhs = 0.0;

      for (int d = 0; d < 2; d++) {
        const int t = 1 - d;
        for (int lado = 0; lado < 2; lado++) {
          const double sinal = (lado == 0) ? -1.0 : 1.0;

          // A DIVERGENCIA DE `u*` SOMA AS SUB-FACETAS, na mesma varredura: a
          // faceta e' a mesma, e percorre-la duas vezes abriria espaco para as
          // duas contas discordarem na borda.
          const std::vector<SubFaceta> sfs = sub_facetas(b, i, d, lado);
          for (const SubFaceta &sf : sfs) {
            const double fluxo = sinal * ustar_.em(d, sf.id) * sf.medida / dt_;
            rhs += fluxo;
            soma_b += fluxo;
          }

          const std::vector<long> vz = m_->ids_vizinhos(b, i, 2 * d + lado);
          if (vz.empty()) {
            // ---- parede ------------------------------------------------
            Ponto2 xp = c;
            xp[d] += sinal * hh[d] / 2.0;
            const Contorno &cp = bcp_->em(xp);
            if (cp.tipo != TipoContorno::Dirichlet) continue;   // Neumann: fluxo nulo
            double val = cp.g(xp);
            if (proj_ == Projecao::Incremental) {
              Ponto2 xcc = c;
              xcc[d] -= sinal * hh[d];
              const long cc = m_->id_do_ponto(xcc);
              const double p_parede =
                  (cc >= 0) ? 1.5 * p_.em(id) - 0.5 * p_.em(cc) : p_.em(id);
              val -= p_parede;
            }
            linha[id] -= 2.0 * hh[t] / hh[d];
            rhs -= 2.0 * hh[t] * val / hh[d];
            continue;
          }

          // ---- face interna, UNIFORME OU EM TRANSICAO 2:1, uma formula ---
          //
          //     |face| = o menor dos dois lados transversais
          //     dist   = (h_self + h_viz) / 2
          //     coef   = |face| / dist
          //
          // UMA formula e nao tres (uniforme, grossa->fina, fina->grossa):
          // tres divergem, e divergir aqui quebra a CONSERVACAO sem quebrar
          // nada visivel.  Ela e' conservativa por construcao -- o coeficiente
          // que a grossa usa para cada fina e' IGUAL ao que a fina usa para a
          // grossa, logo o que sai de uma entra na outra e a matriz fica
          // simetrica.
          const std::vector<int> nvz = m_->niveis_vizinhos(b, i, 2 * d + lado);
          const int nivel_self = m_->nivel(b, i);
          for (size_t k = 0; k < vz.size(); k++) {
            const int dn = nivel_self - (k < nvz.size() ? nvz[k] : nivel_self);
            const double fator = std::pow(2.0, static_cast<double>(dn));
            const double h_viz = hh[d] * fator;
            const double face = hh[t] * std::fmin(1.0, fator);
            const double dist = (hh[d] + h_viz) / 2.0;
            // O COEFICIENTE E' `1/rho` NA FACE, e e' o que faz o operador
            // deixar de ser laplaciano.  Tem de ser o MESMO `rho` que a
            // projecao usa na mesma faceta: se diferirem, o par `div`/`grad`
            // nao fecha e a divergencia sobra -- o mesmo modo de falha que o
            // gradiente transversal com duas fontes ja' produziu.
            const double coef =
                face / (dist * rho_na_face(c, d, static_cast<int>(sinal), hh[d]));

            linha[vz[k]] += coef;
            linha[id] -= coef;

            // ---- A CORRECAO TRANSVERSAL, e ela NAO e' opcional ----------
            //
            // Na transicao 2:1 os dois centros NAO estao alinhados com a
            // normal: o centro da fina esta' deslocado um quarto de celula do
            // centro da grossa.  A diferenca bruta mede a derivada ao longo da
            // LINHA DOS CENTROS, que tem componente transversal -- e para a
            // celula FINA esse erro NAO cancela, porque ela ve' uma unica face
            // contra a grossa.
            //
            // MEDIDO no Poisson deste projeto, sem ela: o grau 1, exato em
            // malha uniforme, deu 3,9e-02 na adaptada, e a trigonometrica teve
            // ordem NEGATIVA.  Nao era erro de O(h): era de O(1).
            //
            // ENTRA COMO COEFICIENTE, e nao como correcao diferida com o
            // incremento do passo anterior: diferida introduziria um erro de
            // ordem `dt` que se atribuiria ao esquema.
            const auto itv = bi_por_id_.find(vz[k]);
            if (itv == bi_por_id_.end()) continue;
            const long bv = itv->second.first, iv = itv->second.second;
            const Ponto2 cv = m_->centroide(bv, iv);
            const double delta = cv[t] - c[t];
            if (std::fabs(delta) < 1e-14 * hh[t]) continue;

            // O GRADIENTE E' AVALIADO NA CELULA GROSSA -- e' ela que precisa
            // ser deslocada para a linha do centro da sub-face.  O SINAL nao
            // depende de quem e' o grosso: refazendo a algebra, o fluxo da'
            // `coef*p_viz - coef*p_self - coef*delta*(dp/dt)` nos dois casos,
            // e o que muda e' so' ONDE a derivada e' avaliada.
            const bool self_grosso = (h_viz < hh[d]);
            const long gb = self_grosso ? b : bv;
            const long gi = self_grosso ? i : iv;

            // GRADIENTE TRANSVERSAL na grossa, como coeficientes.  Cada lado
            // contribui com a media dos vizinhos dele; parede com Neumann
            // homogeneo contribui com a PROPRIA celula, que e' o que
            // `d(dp)/dn = 0` quer dizer.
            // OS MESMOS PESOS que a projecao avalia -- uma fonte so'.  Ter
            // duas implementacoes disto foi o que fez o par `div`/`grad`
            // deixar de fechar perto da parede.
            std::map<long, double> wg;
            pesos_grad_transversal(gb, gi, t, &wg);
            for (const auto &kv : wg) linha[kv.first] -= coef * delta * kv.second;
          }
        }
      }

      if (fonte_celula_) rhs -= fonte_celula_(c) * vol;

      std::vector<int> cols;
      std::vector<double> vals;
      cols.reserve(linha.size());
      vals.reserve(linha.size());
      for (const auto &kv : linha) {
        cols.push_back(static_cast<int>(kv.first));
        vals.push_back(kv.second);
      }
      s.linha(static_cast<int>(id), cols, vals);
      s.lado_direito(static_cast<int>(id), rhs);
    }

  soma_rhs_ = soma_b;
  s.monta();
  s.resolve();
  it_pressao_ = s.iteracoes();
  conv_pressao_ = s.convergiu();

  const std::vector<double> &sol = s.solucao();
  for (long k = 0; k < n; k++) dp_.em(p0 + k) = sol[static_cast<size_t>(k)];
}

void NavierStokes::projeta()
{
  // `u^{n+1} = u* - dt grad(dp)`, e e' AQUI que a divergencia vai a zero.
  //
  // O gradiente usado e' o MESMO operador do preditor, e tem de ser: se a
  // projecao usasse um gradiente diferente do que a divergencia do Poisson
  // supos, o campo corrigido nao teria divergencia nula e nada diria -- o
  // solver convergiria igual.  E' o par `div`/`grad` que fecha, nao cada um.
  double maxd = 0.0;
  for (int d = 0; d < 2; d++)
    for (long f = 0; f < fd_->num_facetas(d); f++) {
      const Ponto2 x = fd_->centro(d, f);
      Ponto2 a = x, b = x;
      a[d] -= h_min_[d] / 4.0;
      b[d] += h_min_[d] / 4.0;
      if (m_->id_do_ponto(a) < 0 || m_->id_do_ponto(b) < 0) continue;

      // O GRADIENTE AQUI E' O MESMO OPERADOR DO PREDITOR, e tem de ser.
      //
      // Se a projecao usasse um gradiente diferente do que a divergencia do
      // Poisson supos, o campo corrigido NAO teria divergencia nula e nada
      // diria -- o solver convergiria igual.  E' o par `div`/`grad` que fecha,
      // nao cada um.  Por isso a mesma funcao, com outro campo.
      const double antes = u_.em(d, f);
      u_.em(d, f) = ustar_.em(d, f) - dt_ * grad_campo(dp_, d, f) / rho_na_faceta(d, f);
      maxd = std::max(maxd, std::fabs(u_.em(d, f) - antes));
    }
  resid_estac_ = maxd / dt_;

  for (long b = 0; b < m_->num_blocos(); b++)
    for (long i = 0; i < m_->num_folhas_do_bloco(b); i++) {
      const long id = m_->id_global(b, i);
      p_.em(id) = (proj_ == Projecao::Incremental) ? p_.em(id) + dp_.em(id)
                                                   : dp_.em(id);
    }
}

void NavierStokes::passo()
{
  if (!(dt_ > 0.0))
    throw std::runtime_error("NavierStokes::passo: o passo de tempo nao foi "
                             "instalado -- chame passo_de_tempo primeiro");

  // A TRANSICAO 2:1 NAO E' MAIS RECUSADA.  Os tres termos -- convectivo,
  // difusivo e gradiente de pressao -- e o Poisson da projecao tratam faces
  // de tamanhos diferentes pela mesma formula que trata as iguais, e os
  // casos uniformes reproduzem BIT A BIT o que davam antes: o caminho geral
  // reduz ao particular, logo eles verificam o codigo adaptado.

  // O RELOGIO DO CONTORNO: `t^n` para o preditor, `t^{n+1}` para a aplicacao
  // final.  Ver `instante_do_contorno` no cabecalho.
  t_bc_ = t_;
  if (relogio_) *relogio_ = t_bc_;
  aplica_contorno_velocidade();
  // O AVANCO TEMPORAL ESCOLHE O PREDITOR, e e' a unica coisa que muda entre os
  // dois: projecao, Poisson e correcao sao identicos.  Com `dt` pequeno os
  // dois convergem para a mesma resposta, e e' isso que o arreio exige.
  if (tempo_ == Tempo::EulerSemiImplicito)
    preditor_semi_implicito();
  else
    preditor();
  // O CONTORNO VALE PARA `u*` TAMBEM: o Poisson usa a divergencia de `u*`, e
  // se a componente normal na parede nao for a da condicao, a compatibilidade
  // do Poisson quebra -- com Neumann em toda a volta ele fica insoluvel, e com
  // Dirichlet em algum lugar ele resolve e da' a resposta errada.
  for (int d = 0; d < 2; d++)
    for (long f = 0; f < fd_->num_facetas(d); f++) {
      const Ponto2 x = fd_->centro(d, f);
      Ponto2 a = x, b = x;
      a[d] -= h_min_[d] / 4.0;
      b[d] += h_min_[d] / 4.0;
      if (m_->id_do_ponto(a) >= 0 && m_->id_do_ponto(b) >= 0) continue;
      ustar_.em(d, f) = u_.em(d, f);
    }

  resolve_pressao();
  projeta();
  t_bc_ = t_ + dt_;
  if (relogio_) *relogio_ = t_bc_;
  aplica_contorno_velocidade();

  t_ += dt_;
  n_passos_++;
}

std::vector<NavierStokes::SubFaceta> NavierStokes::sub_facetas(long b, long i,
                                                                 int d,
                                                                 int lado) const
{
  const Ponto2 c = m_->centroide(b, i);
  const Ponto2 hh = m_->lados_da_folha(b, i);
  const int t = 1 - d;                     // direcao TRANSVERSAL a' face
  const double sinal = (lado == 0) ? -1.0 : 1.0;

  // QUANTAS SUB-FACETAS, pelo NIVEL do vizinho e nao por busca.
  //
  // Com 2:1 a razao e' `2^(nivel_viz - nivel_self)`, e `niveis_vizinhos` ja' a
  // da'.  Procurar as facetas por varredura seria O(facetas) por face, dentro
  // do laco de montagem -- o mesmo quadratico que ja' custou caro.
  const std::vector<int> nvz = m_->niveis_vizinhos(b, i, 2 * d + lado);
  const int nivel_self = m_->nivel(b, i);
  int nsub = 1;
  double larg = hh[t];
  if (!nvz.empty() && nvz[0] > nivel_self) {
    nsub = 1 << (nvz[0] - nivel_self);     // 2 em 2:1
    larg = hh[t] / nsub;
  }

  std::vector<SubFaceta> out;
  out.reserve(static_cast<size_t>(nsub));
  for (int sub = 0; sub < nsub; sub++) {
    Ponto2 cf = c;
    cf[d] += sinal * hh[d] / 2.0;
    cf[t] += -hh[t] / 2.0 + (sub + 0.5) * larg;
    const long fid = faceta_em(d, cf);
    if (fid < 0)
      throw std::runtime_error(
          "NavierStokes::sub_facetas: nao ha' faceta da direcao " +
          std::to_string(d) + " em (" + std::to_string(cf[0]) + ", " +
          std::to_string(cf[1]) + "), sub " + std::to_string(sub) + " de " +
          std::to_string(nsub) + ".  Uma face sem faceta tem o fluxo DESCARTADO,"
          " e descartado em silencio -- a divergencia desta celula sairia errada"
          " sem nada acusar.");
    out.push_back({fid, larg, cf});
  }
  return out;
}

double NavierStokes::divergencia_em(long b, long i) const
{
  // A MESMA CONTA de `divergencia_maxima`, por celula.  Ela toma o maximo
  // disto.
  const Ponto2 hh = m_->lados_da_folha(b, i);
  const double vol = hh[0] * hh[1];
  double fluxo = 0.0;
  for (int d = 0; d < 2; d++)
    for (int lado = 0; lado <= 1; lado++) {
      const double sinal = (lado == 0) ? -1.0 : 1.0;
      for (const SubFaceta &sf : sub_facetas(b, i, d, lado))
        fluxo += sinal * u_.em(d, sf.id) * sf.medida;
    }
  return fluxo / vol;
}

double NavierStokes::volume_dual_da_faceta(int d, long fid) const
{
  // `medida x distancia entre os centros`.  A distancia e' a MESMA que
  // `grad_campo` usa -- a soma dos dois meios-lados --, e tem de ser: o peso
  // so' faz sentido se for o da forma bilinear que o gradiente define.
  const Ponto2 xf = fd_->centro(d, fid);
  double soma_meios = 0.0;
  int achados = 0;
  for (int lado = 0; lado < 2; lado++) {
    Ponto2 p = xf;
    p[d] += (lado == 0 ? -1.0 : 1.0) * h_min_[d] / 4.0;
    const long id = m_->id_do_ponto(p);
    const auto it = (id >= 0) ? bi_por_id_.find(id) : bi_por_id_.end();
    if (it == bi_por_id_.end()) continue;
    soma_meios += m_->lados_da_folha(it->second.first, it->second.second)[d]
                  / 2.0;
    achados++;
  }
  // FACETA DE CONTORNO devolve zero, e nao meia distancia: `grad_campo`
  // devolve zero ali ("nao se resolve ali"), entao o peso tem de anula-la na
  // forma bilinear tambem -- senao os dois lados da adjuncao somariam termos
  // diferentes.
  if (achados < 2) return 0.0;
  return fd_->medida(d, fid) * soma_meios;
}

double NavierStokes::divergencia_maxima() const
{
  // O FLUXO PELA FACE E' A SOMA DAS SUB-FACETAS, e nao o valor numa delas.
  //
  // Na transicao 2:1 a face da celula grossa e' coberta por DUAS facetas finas,
  // e o fluxo total e' `u_1 . larg + u_2 . larg`.  Tomar o valor no centro da
  // face grossa -- que nao e' o centro de nenhuma das duas -- daria -1 na busca
  // e o fluxo seria descartado.
  //
  // E A DIVERGENCIA E' DIVIDIDA PELO VOLUME, nao pelo lado: com faces de
  // tamanhos diferentes na mesma celula, `soma(u)/h` nao e' divergencia
  // nenhuma.
  double maxv = 0.0;
  for (long b = 0; b < m_->num_blocos(); b++)
    for (long i = 0; i < m_->num_folhas_do_bloco(b); i++) {
      const Ponto2 hh = m_->lados_da_folha(b, i);
      const double vol = hh[0] * hh[1];
      double fluxo = 0.0;
      for (int d = 0; d < 2; d++)
        for (int lado = 0; lado <= 1; lado++) {
          const double sinal = (lado == 0) ? -1.0 : 1.0;
          for (const SubFaceta &sf : sub_facetas(b, i, d, lado))
            fluxo += sinal * u_.em(d, sf.id) * sf.medida;
        }
      maxv = std::max(maxv, std::fabs(fluxo / vol));
    }
  return maxv;
}

// ===========================================================================
// AS PRIMITIVAS DO ESTENCIL EM MALHA ADAPTADA
// ===========================================================================

double NavierStokes::media_de_fluxo(long b, long i, int d, int lado) const
{
  double num = 0.0, den = 0.0;
  for (const SubFaceta &sf : sub_facetas(b, i, d, lado)) {
    num += u_.em(d, sf.id) * sf.medida;
    den += sf.medida;
  }
  return (den > 0.0) ? num / den : 0.0;
}

double NavierStokes::u_no_ponto(int d, const Ponto2 &x, bool *dentro) const
{
  // CAMINHO CURTO: ha' faceta exatamente aqui.  E' o unico que roda em malha
  // uniforme, e e' por isso que os casos uniformes verificam este codigo.
  const long f = faceta_em(d, x);
  if (f >= 0) {
    if (dentro) *dentro = true;
    return u_.em(d, f);
  }

  // NAO HA' FACETA: o ponto cai no MEIO de uma celula na direcao `d` -- e' o
  // que acontece com a face de uma fina contra a grossa vizinha.  O valor sai
  // da interpolacao linear entre as DUAS faces de direcao `d` dessa celula.
  const long id = m_->id_do_ponto(x);
  if (id < 0) {
    if (dentro) *dentro = false;
    return 0.0;
  }
  const auto it = bi_por_id_.find(id);
  if (it == bi_por_id_.end()) {
    if (dentro) *dentro = false;
    return 0.0;
  }
  const long b = it->second.first, i = it->second.second;
  const Ponto2 c = m_->centroide(b, i);
  const Ponto2 hh = m_->lados_da_folha(b, i);

  auto dentro_da_celula = [&](long bb, long ii) {
    const Ponto2 cc = m_->centroide(bb, ii);
    const Ponto2 h2 = m_->lados_da_folha(bb, ii);
    const double vm = media_de_fluxo(bb, ii, d, 0);
    const double vp = media_de_fluxo(bb, ii, d, 1);
    double lam = (x[d] - (cc[d] - h2[d] / 2.0)) / h2[d];
    if (lam < 0.0) lam = 0.0;
    if (lam > 1.0) lam = 1.0;
    return (1.0 - lam) * vm + lam * vp;
  };

  const double v_aqui = dentro_da_celula(b, i);

  // E A INTERPOLACAO TRANSVERSAL, que faltava e nao e' detalhe.
  //
  // As facetas de direcao `d` ficam no centro da celula na direcao `q = 1-d`.
  // Um ponto FORA desse centro -- a QUINA do volume de controle, por exemplo,
  // que e' onde a velocidade que advecta e' avaliada -- nao e' representado
  // por interpolacao so' ao longo de `d`.
  //
  // MEDIDO ao esquecer disto: a quina devolvia a faceta de UMA das duas
  // celulas vizinhas em vez da media das duas, e o erro da caixa manufaturada
  // subiu 20% em `u` e 37% em `v`.  O canal mal se mexeu -- nele `v` e' nulo e
  // o convectivo se anula --, que e' exatamente por que ele nao serve para
  // verificar isto.
  const int q = 1 - d;
  const double delta = x[q] - c[q];
  if (std::fabs(delta) < 1e-12 * hh[q]) {
    if (dentro) *dentro = true;
    return v_aqui;
  }

  Ponto2 pv = x;
  pv[q] = c[q] + (delta > 0.0 ? 1.0 : -1.0) * (hh[q] / 2.0 + h_min_[q] / 4.0);
  pv[d] = c[d];                        // dentro da faixa, para a busca nao errar
  const long idv = m_->id_do_ponto(pv);
  const auto itv = (idv >= 0) ? bi_por_id_.find(idv) : bi_por_id_.end();
  if (itv == bi_por_id_.end()) {        // sem vizinha: nao ha' o que interpolar
    if (dentro) *dentro = true;
    return v_aqui;
  }
  const long bv = itv->second.first, iv = itv->second.second;
  const double hv = m_->lados_da_folha(bv, iv)[q];
  const double v_viz = dentro_da_celula(bv, iv);
  double mu = std::fabs(delta) / ((hh[q] + hv) / 2.0);
  if (mu > 1.0) mu = 1.0;
  if (dentro) *dentro = true;
  return (1.0 - mu) * v_aqui + mu * v_viz;
}

NavierStokes::Vizinha NavierStokes::vizinha_de(int d, long fid, int e,
                                               int s) const
{
  const Dono &dn = dono_[d][static_cast<size_t>(fid)];
  const Ponto2 xf = fd_->centro(d, fid);
  const Ponto2 hC = m_->lados_da_folha(dn.b, dn.i);
  const Ponto2 cC = m_->centroide(dn.b, dn.i);
  Vizinha v;

  if (e == d) {
    // AO LONGO DA NORMAL: o vizinho e' a OUTRA FACE de uma celula, e a
    // distancia e' o lado DELA.  Isso vale sempre -- entre duas faces opostas
    // da mesma celula a distancia e' exatamente o lado, em qualquer nivel.
    const bool para_dentro = (dn.lado == 0) ? (s > 0) : (s < 0);
    if (para_dentro) {
      v.valor = media_de_fluxo(dn.b, dn.i, d, 1 - dn.lado);
      v.dist = hC[d];
      v.dentro = true;
      return v;
    }
    // Para fora da celula dona: atravessa a do outro lado da faceta.
    Ponto2 p = xf;
    p[d] += s * h_min_[d] / 4.0;          // um cutucao, so' para cair dentro
    const long id = m_->id_do_ponto(p);
    const auto it = (id >= 0) ? bi_por_id_.find(id) : bi_por_id_.end();
    if (it == bi_por_id_.end()) { v.dentro = false; return v; }
    const long bb = it->second.first, ii = it->second.second;
    // A faceta e' a face "menos" dessa celula quando se anda em +d.
    const int lado_perto = (s > 0) ? 0 : 1;
    v.valor = media_de_fluxo(bb, ii, d, 1 - lado_perto);
    v.dist = m_->lados_da_folha(bb, ii)[d];
    v.dentro = true;
    return v;
  }

  // TRANSVERSALMENTE: atravessa a face de `t` da celula dona.  A distancia e'
  // a MEIA-SOMA dos lados das duas -- a distancia entre os centros --, e nao o
  // lado de nenhuma delas.
  const int t = e;
  Ponto2 p = xf;
  // O cutucao em `d` poe o ponto DENTRO da faixa da celula dona: `xf` esta'
  // exatamente sobre a face, e sem isso a busca pode cair na celula errada.
  p[d] += (dn.lado == 0 ? 1.0 : -1.0) * h_min_[d] / 4.0;
  p[t] = cC[t] + s * (hC[t] / 2.0 + h_min_[t] / 4.0);
  const long id = m_->id_do_ponto(p);
  const auto it = (id >= 0) ? bi_por_id_.find(id) : bi_por_id_.end();
  if (it == bi_por_id_.end()) { v.dentro = false; return v; }
  const long bb = it->second.first, ii = it->second.second;
  const Ponto2 hN = m_->lados_da_folha(bb, ii);
  const Ponto2 cN = m_->centroide(bb, ii);

  // O ALVO esta' na coordenada `d` DA FACETA e na linha de centro da vizinha.
  Ponto2 alvo;
  alvo[d] = xf[d];
  alvo[t] = cN[t];
  bool ok = false;
  v.valor = u_no_ponto(d, alvo, &ok);
  v.dist = (hC[t] + hN[t]) / 2.0;
  v.dentro = ok;
  return v;
}

void NavierStokes::pesos_grad_transversal(long b, long i, int t,
                                          std::map<long, double> *w) const
{
  // OS PESOS DO GRADIENTE TRANSVERSAL, e esta e' a UNICA fonte deles.
  //
  // Eles aparecem em DOIS lugares -- como coeficientes na matriz do Poisson e
  // como avaliacao na projecao --, e tinham DUAS implementacoes.  Elas
  // concordavam no interior e discordavam NA PAREDE: uma fechava com o valor
  // da propria celula a `h/2` (que e' o que `d(dp)/dn = 0` quer dizer), a
  // outra usava a formula lateral sem esse ponto.
  //
  // O EFEITO: o par `div`/`grad` deixava de fechar, e a divergencia depois da
  // projecao ia de 1e-14 para 9,5e-06.  E so' aparecia quando o refino chegava
  // perto da parede -- com o bloco refinado no interior, as duas concordavam e
  // o defeito ficava invisivel.
  //
  // Nao e' erro de implementacao: e' erro de TER DUAS.
  w->clear();
  const Ponto2 c = m_->centroide(b, i);
  const Ponto2 hh = m_->lados_da_folha(b, i);
  const long id = m_->id_global(b, i);

  std::map<long, double> lado_v[2];
  double dist_l[2] = {0.0, 0.0};
  for (int lt = 0; lt < 2; lt++) {
    const std::vector<long> vt = m_->ids_vizinhos(b, i, 2 * t + lt);
    if (vt.empty()) {
      // PAREDE com Neumann homogeneo: o valor ali e' o da propria celula, e
      // esta' a meio lado de distancia.
      lado_v[lt][id] = 1.0;
      dist_l[lt] = hh[t] / 2.0;
      continue;
    }
    double soma_d = 0.0;
    int nn = 0;
    for (long q : vt) {
      const auto itq = bi_por_id_.find(q);
      if (itq == bi_por_id_.end()) continue;
      lado_v[lt][q] += 1.0;
      soma_d += std::fabs(
          m_->centroide(itq->second.first, itq->second.second)[t] - c[t]);
      nn++;
    }
    if (nn == 0) { lado_v[lt][id] = 1.0; dist_l[lt] = hh[t] / 2.0; continue; }
    for (auto &kv : lado_v[lt]) kv.second /= nn;
    dist_l[lt] = soma_d / nn;
  }
  const double dtot = dist_l[0] + dist_l[1];
  if (dtot <= 0.0) return;
  for (const auto &kv : lado_v[1]) (*w)[kv.first] += kv.second / dtot;
  for (const auto &kv : lado_v[0]) (*w)[kv.first] -= kv.second / dtot;
}

double NavierStokes::grad_transversal(const CampoEscalar &campo, long b,
                                      long i, int t) const
{
  // A AVALIACAO dos pesos acima.  Nao reimplementa nada -- se reimplementasse,
  // voltaria a divergir da matriz do Poisson, que e' o defeito que custou a
  // divergencia de 9,5e-06.
  std::map<long, double> w;
  pesos_grad_transversal(b, i, t, &w);
  double g = 0.0;
  for (const auto &kv : w) g += kv.second * campo.em(kv.first);
  return g;
}

// ===========================================================================
// PROPRIEDADES VARIAVEIS
// ===========================================================================

double NavierStokes::rho_em(const Ponto2 &x) const
{
  // SEM CAMPO, densidade UM.  E' o caminho monofasico, e ele tem de ser exato
  // e nao "equivalente": e' por ele que os casos monofasicos continuam
  // verificando este codigo depois da mudanca.
  if (rho_ == nullptr) return 1.0;
  const long id = m_->id_do_ponto(x);
  if (id < 0) return 1.0;
  return rho_->em(id);
}

double NavierStokes::mu_em(const Ponto2 &x) const
{
  if (mu_ == nullptr) return 1.0 / re_;
  const long id = m_->id_do_ponto(x);
  if (id < 0) return 1.0 / re_;
  return mu_->em(id);
}

double NavierStokes::rho_na_face(const Ponto2 &c, int d, int s2,
                                 double hd) const
{
  // A MEDIA ARITMETICA das duas celulas, como no HiGFlow.
  //
  // E TEM DE SER A MESMA que a projecao usa na mesma faceta.  Se o Poisson e a
  // projecao usarem densidades diferentes na mesma face, o par `div`/`grad`
  // nao fecha e a divergencia sobra -- o modo de falha que o gradiente
  // transversal com duas fontes ja' produziu aqui, e que custou 9,5e-06.
  if (rho_ == nullptr) return 1.0;
  const long a = m_->id_do_ponto(c);
  Ponto2 q = c;
  q[d] += s2 * hd;
  const long b = m_->id_do_ponto(q);
  const double ra = (a >= 0) ? rho_->em(a) : 1.0;
  const double rb = (b >= 0) ? rho_->em(b) : ra;
  return 0.5 * (ra + rb);
}

double NavierStokes::rho_na_faceta(int d, long fid) const
{
  if (rho_ == nullptr) return 1.0;
  const Ponto2 xf = fd_->centro(d, fid);
  double r[2];
  int n = 0;
  for (int lado = 0; lado < 2; lado++) {
    Ponto2 p = xf;
    p[d] += (lado == 0 ? -1.0 : 1.0) * h_min_[d] / 4.0;
    const long id = m_->id_do_ponto(p);
    r[lado] = (id >= 0) ? rho_->em(id) : 0.0;
    if (id >= 0) n++;
  }
  if (n == 2) return 0.5 * (r[0] + r[1]);
  if (n == 1) return r[0] + r[1];
  return 1.0;
}

double NavierStokes::mu_na_quina(const Ponto2 &xq) const
{
  // A MEDIA HARMONICA DAS QUATRO CELULAS que cercam a quina, como no HiGFlow.
  //
  //     mu = 4 / (1/v1 + 1/v2 + 1/v3 + 1/v4)
  //
  // Harmonica porque e' ela que preserva a continuidade da tensao cisalhante
  // atraves de uma interface plana.  Com a aritmetica, a tensao salta ali -- e
  // "ali" e' justamente onde os dois fluidos se tocam.
  //
  // Uma celula que nao exista -- quina no contorno -- simplesmente nao entra,
  // e a media e' das que existem.  Trata-la como viscosidade zero faria a
  // harmonica ir a zero e a parede virar livre de tensao.
  if (mu_ == nullptr) return 1.0 / re_;
  double soma = 0.0;
  int n = 0;
  for (int sy = -1; sy <= 1; sy += 2)
    for (int sx = -1; sx <= 1; sx += 2) {
      const Ponto2 q{xq[0] + sx * h_min_[0] / 4.0, xq[1] + sy * h_min_[1] / 4.0};
      const long id = m_->id_do_ponto(q);
      if (id < 0) continue;
      const double v = mu_->em(id);
      if (v <= 0.0) continue;
      soma += 1.0 / v;
      n++;
    }
  if (n == 0) return 1.0 / re_;
  return n / soma;
}

// ===========================================================================
// OS PESOS: a mesma dependencia linear, em forma de MATRIZ.
//
// Cada funcao aqui e' o espelho de uma funcao de VALOR, e o contrato entre as
// duas e' `valor == soma_j w_j u_j (+ cte)`.  O arreio exige isso em precisao
// de maquina por `difusao_explicita` contra `difusao_por_pesos`.
//
// POR QUE DUAS E NAO UMA.  Colapsar as duas -- fazer a avaliacao explicita
// somar os pesos -- trocaria a ORDEM DAS SOMAS de ponto flutuante, e com ela
// os ultimos bits de todo caso monofasico ja' gravado.  O criterio de
// "reproduz bit a bit" e' o que pegou duas das minhas refatoracoes mentindo,
// e perde-lo para ganhar uma funcao a menos e' troca ruim.  O preco e' a
// gemea, e a gemea e' testada.
// ===========================================================================

void NavierStokes::pesos_media_de_fluxo(long b, long i, int d, int lado,
                                        double esc, Pesos *w) const
{
  const std::vector<SubFaceta> sfs = sub_facetas(b, i, d, lado);
  double den = 0.0;
  for (const SubFaceta &sf : sfs) den += sf.medida;
  if (den <= 0.0) return;                 // face sem faceta: valor zero
  for (const SubFaceta &sf : sfs) (*w)[sf.id] += esc * sf.medida / den;
}

bool NavierStokes::pesos_u_no_ponto(int d, const Ponto2 &x, double esc,
                                    Pesos *w) const
{
  const long f = faceta_em(d, x);
  if (f >= 0) { (*w)[f] += esc; return true; }

  const long id = m_->id_do_ponto(x);
  if (id < 0) return false;
  const auto it = bi_por_id_.find(id);
  if (it == bi_por_id_.end()) return false;
  const long b = it->second.first, i = it->second.second;
  const Ponto2 c = m_->centroide(b, i);
  const Ponto2 hh = m_->lados_da_folha(b, i);

  auto pesos_da_celula = [&](long bb, long ii, double sc) {
    const Ponto2 cc = m_->centroide(bb, ii);
    const Ponto2 h2 = m_->lados_da_folha(bb, ii);
    double lam = (x[d] - (cc[d] - h2[d] / 2.0)) / h2[d];
    if (lam < 0.0) lam = 0.0;
    if (lam > 1.0) lam = 1.0;
    pesos_media_de_fluxo(bb, ii, d, 0, sc * (1.0 - lam), w);
    pesos_media_de_fluxo(bb, ii, d, 1, sc * lam, w);
  };

  const int q = 1 - d;
  const double delta = x[q] - c[q];
  if (std::fabs(delta) < 1e-12 * hh[q]) {
    pesos_da_celula(b, i, esc);
    return true;
  }

  Ponto2 pv = x;
  pv[q] = c[q] + (delta > 0.0 ? 1.0 : -1.0) * (hh[q] / 2.0 + h_min_[q] / 4.0);
  pv[d] = c[d];
  const long idv = m_->id_do_ponto(pv);
  const auto itv = (idv >= 0) ? bi_por_id_.find(idv) : bi_por_id_.end();
  if (itv == bi_por_id_.end()) {
    pesos_da_celula(b, i, esc);
    return true;
  }
  const long bv = itv->second.first, iv = itv->second.second;
  const double hv = m_->lados_da_folha(bv, iv)[q];
  double mu = std::fabs(delta) / ((hh[q] + hv) / 2.0);
  if (mu > 1.0) mu = 1.0;
  pesos_da_celula(b, i, esc * (1.0 - mu));
  pesos_da_celula(bv, iv, esc * mu);
  return true;
}

NavierStokes::VizinhaPesos NavierStokes::pesos_vizinha_de(int d, long fid,
                                                          int e, int s,
                                                          double esc,
                                                          Pesos *w) const
{
  const Dono &dn = dono_[d][static_cast<size_t>(fid)];
  const Ponto2 xf = fd_->centro(d, fid);
  const Ponto2 hC = m_->lados_da_folha(dn.b, dn.i);
  const Ponto2 cC = m_->centroide(dn.b, dn.i);
  VizinhaPesos v;

  if (e == d) {
    const bool para_dentro = (dn.lado == 0) ? (s > 0) : (s < 0);
    if (para_dentro) {
      pesos_media_de_fluxo(dn.b, dn.i, d, 1 - dn.lado, esc, w);
      v.dist = hC[d];
      v.dentro = true;
      return v;
    }
    Ponto2 p = xf;
    p[d] += s * h_min_[d] / 4.0;
    const long id = m_->id_do_ponto(p);
    const auto it = (id >= 0) ? bi_por_id_.find(id) : bi_por_id_.end();
    if (it == bi_por_id_.end()) { v.dentro = false; return v; }
    const long bb = it->second.first, ii = it->second.second;
    const int lado_perto = (s > 0) ? 0 : 1;
    pesos_media_de_fluxo(bb, ii, d, 1 - lado_perto, esc, w);
    v.dist = m_->lados_da_folha(bb, ii)[d];
    v.dentro = true;
    return v;
  }

  const int t = e;
  Ponto2 p = xf;
  p[d] += (dn.lado == 0 ? 1.0 : -1.0) * h_min_[d] / 4.0;
  p[t] = cC[t] + s * (hC[t] / 2.0 + h_min_[t] / 4.0);
  const long id = m_->id_do_ponto(p);
  const auto it = (id >= 0) ? bi_por_id_.find(id) : bi_por_id_.end();
  if (it == bi_por_id_.end()) { v.dentro = false; return v; }
  const long bb = it->second.first, ii = it->second.second;
  const Ponto2 hN = m_->lados_da_folha(bb, ii);
  const Ponto2 cN = m_->centroide(bb, ii);

  Ponto2 alvo;
  alvo[d] = xf[d];
  alvo[t] = cN[t];
  v.dentro = pesos_u_no_ponto(d, alvo, esc, w);
  v.dist = (hC[t] + hN[t]) / 2.0;
  return v;
}

void NavierStokes::contorno_afim(int d, long fid, int e, int s, double *A,
                                 double *B) const
{
  // ESPELHO EXATO de `contorno_vizinho`: fantasma = `A * u_f + B`.
  const Dono &dn = dono_[d][static_cast<size_t>(fid)];
  const Ponto2 hC = m_->lados_da_folha(dn.b, dn.i);
  const Ponto2 xf = fd_->centro(d, fid);
  *A = 1.0;
  *B = 0.0;
  if (e == d) return;                     // nao ocorre em faceta interior

  Ponto2 xp = xf;
  xp[e] += s * hC[e] / 2.0;
  const Contorno &c = bcu_[d]->em(xp);
  switch (c.tipo) {
    case TipoContorno::Dirichlet:
      *A = -1.0;
      *B = 2.0 * c.g(xp);
      return;
    case TipoContorno::Neumann:
      *A = 1.0;
      *B = s * hC[e] * c.g(xp);
      return;
    case TipoContorno::Robin:
      c.afim(xp, hC[e], A, B);
      return;
  }
}

void NavierStokes::pesos_difusao(int d, long fid, Pesos *w, double *cte) const
{
  // ESPELHO DE `tensao_viscosa_fluxo`.  Mesma escolha de `mu`, mesmas
  // distancias, mesmo divisor: o que muda e' so' que o vizinho entra como
  // DEPENDENCIA em vez de numero.
  const Ponto2 xf = fd_->centro(d, fid);
  *cte = 0.0;

  for (int e = 0; e < 2; e++) {
    // O divisor e' `(dist_- + dist_+)/2` e so' se conhece depois dos dois
    // lados -- por isso a direcao acumula num mapa proprio antes de dividir.
    Pesos wn;
    double cn = 0.0;
    double dist[2] = {1.0, 1.0};

    for (int lado = 0; lado < 2; lado++) {
      const int sgn = (lado == 0) ? -1 : +1;
      Pesos wv;
      const VizinhaPesos v = pesos_vizinha_de(d, fid, e, sgn, 1.0, &wv);
      double dd = 0.0, A = 0.0, B = 0.0;
      if (v.dentro) {
        dd = v.dist;
      } else {
        dd = dist_ate_parede(d, fid, e, sgn);
        wv.clear();
        contorno_afim(d, fid, e, sgn, &A, &B);
      }
      if (dd <= 0.0) { dist[lado] = 1.0; continue; }

      Ponto2 xm = xf;
      xm[e] += sgn * dd / 2.0;
      const double mu_m = (e == d) ? mu_em(xm) : mu_na_quina(xm);
      const double coef = mu_m / dd;

      // `fluxo = mu (u_viz - u_aqui) / dd`
      if (v.dentro)
        for (const auto &kv : wv) wn[kv.first] += coef * kv.second;
      else {
        wn[fid] += coef * A;
        cn += coef * B;
      }
      wn[fid] -= coef;
      dist[lado] = dd;
    }

    const double D = (dist[0] + dist[1]) / 2.0;
    for (const auto &kv : wn) (*w)[kv.first] += kv.second / D;
    *cte += cn / D;
  }
}

double NavierStokes::difusao_explicita(int d, long fid) const
{
  return tensao_viscosa_fluxo(d, fid);
}

double NavierStokes::difusao_por_pesos(int d, long fid) const
{
  Pesos w;
  double cte = 0.0;
  pesos_difusao(d, fid, &w, &cte);
  double v = cte;
  for (const auto &kv : w) v += kv.second * u_.em(d, kv.first);
  return v;
}

void NavierStokes::preditor_semi_implicito()
{
  // `(I - (dt/rho) L_mu) u* = u^n + dt R`, UMA COMPONENTE POR VEZ.
  //
  // O QUE ESTA' IMPLICITO E' SO' `div(mu grad u)`, como no HiGFlow
  // (`higflow_semi_implicit_euler_intermediate_velocity_multiphase`, com
  // `wr = -dt viscr/(Re h^2)/dens`).  O convectivo, a gravidade, a forca de
  // tensao, o termo fonte e o TRANSPOSTO ficam explicitos.  O transposto em
  // particular NAO pode ir para a matriz: ele depende da derivada CRUZADA
  // `d_d u_e`, e isso acoplaria as duas componentes num sistema so' -- o dobro
  // de incognitas por um termo que e' pequeno fora da interface.
  //
  // O GANHO E' O LIMITE VISCOSO DO PASSO.  No caso de Hysing o `dt` explicito
  // sai de `min(visc, capilar, convect)`: no `h = 0,00625` isso da' visc
  // 9,77e-4 contra capilar 9,34e-4 -- os dois apertados. Com o viscoso
  // implicito sobra o capilar, que e' propriedade da tensao superficial e nao
  // da difusao, e nao se remove com implicitacao.
  //
  // E' SERIAL, e isso esta' dito e nao suposto: `DominioFacetas` numera as
  // facetas LOCALMENTE, sem deslocamento global, de modo que o indice da
  // faceta e' a linha da matriz so' em `np = 1`.  O exemplo bifasico ja' recusa
  // `np > 1`; quando a troca de franja dos campos entrar, esta funcao precisa
  // de numeracao global antes de rodar em paralelo.
  ustar_ = u_;
  rebaixos_conv_ = 0;

  for (int d = 0; d < 2; d++) {
    const long nf = fd_->num_facetas(d);
    t8flow::Solver s(t8flow::Motor::Petsc, 0, static_cast<std::size_t>(nf), 9);
    s.tolerancia(1e-12);

    for (long f = 0; f < nf; f++) {
      const Ponto2 x = fd_->centro(d, f);
      Ponto2 a = x, b = x;
      a[d] -= h_min_[d] / 4.0;
      b[d] += h_min_[d] / 4.0;

      // FACETA DE CONTORNO: LINHA IDENTIDADE.  O preditor explicito nao a
      // toca -- ali a velocidade vem da condicao --, e o implicito tem de
      // fazer o mesmo.  Resolver a equacao nela sobrescreveria a condicao de
      // contorno sem que nada acusasse.
      if (m_->id_do_ponto(a) < 0 || m_->id_do_ponto(b) < 0) {
        const std::vector<int> col{static_cast<int>(f)};
        const std::vector<double> val{1.0};
        s.linha(static_cast<int>(f), col, val);
        s.lado_direito(static_cast<int>(f), u_.em(d, f));
        continue;
      }

      const double rho_f = rho_na_faceta(d, f);

      double rhs = -conveccao(d, f) + g_[d]
                   + (tensao_viscosa_transposto(d, f)
                      + (forca_campo_ ? forca_campo_->em(d, f) : 0.0))
                         / rho_f;
      if (proj_ == Projecao::Incremental) rhs -= grad_campo(p_, d, f) / rho_f;
      if (fonte_faceta_) rhs += fonte_faceta_(x, d, t_bc_);

      Pesos w;
      double cte = 0.0;
      pesos_difusao(d, f, &w, &cte);
      const double fat = dt_ / rho_f;

      std::map<long, double> linha;
      linha[f] += 1.0;
      for (const auto &kv : w) linha[kv.first] -= fat * kv.second;

      std::vector<int> cols;
      std::vector<double> vals;
      cols.reserve(linha.size());
      vals.reserve(linha.size());
      for (const auto &kv : linha) {
        cols.push_back(static_cast<int>(kv.first));
        vals.push_back(kv.second);
      }
      s.linha(static_cast<int>(f), cols, vals);
      s.lado_direito(static_cast<int>(f),
                     u_.em(d, f) + dt_ * rhs + fat * cte);
    }

    s.monta();
    s.resolve();
    if (!s.convergiu())
      throw std::runtime_error(
          "NavierStokes: o sistema viscoso da componente " +
          std::to_string(d) + " nao convergiu");
    const std::vector<double> &sol = s.solucao();
    for (long f = 0; f < nf; f++) ustar_.em(d, f) = sol[static_cast<size_t>(f)];
  }
}

double NavierStokes::tensao_viscosa(int d, long fid) const
{
  // `div(mu grad u) + (grad u^T) . grad(mu)`, em FORMA DE FLUXO e com `mu` no
  // MEIO do caminho.  Quem chama divide por `rho`.
  //
  // E' O TENSOR COMPLETO, E E' TAMBEM A FORMA DO HIGFLOW.  Aqui estava escrito
  // o contrario -- que o HiGFlow usava `div(mu grad u)` e que o transposto
  // ficava de fora por comparabilidade.  Isso estava ERRADO, e o erro nao era
  // de opiniao: era de leitura do codigo deles.  O que o HiGFlow descarta nao
  // e' o transposto inteiro.  Abrindo,
  //
  //     div(mu grad u^T)_i = d_j ( mu d_i u_j )
  //                        = (d_j mu)(d_i u_j)  +  mu d_i (d_j u_j)
  //
  // o SEGUNDO termo e' `mu grad(div u)` e some num campo solenoidal -- e' esse
  // que eles omitem, dizendo-o no comentario ("div(Du^t) has to be omitted").
  // O PRIMEIRO sobrevive, e esta' la', em `hig-flow-discret.c`:
  //
  //     ns->cc.d2vdx2[dim2] = dujdxi*dviscdxj;
  //
  // somado em `higflow_extra_diffusive_term` e dividido por `Re` e por `rho`.
  //
  // DESCARTAR O TRANSPOSTO INTEIRO JOGA FORA A TENSAO CISALHANTE ATRAVES DA
  // INTERFACE, que e' exatamente onde o caso bifasico vive: com razao 10 de
  // viscosidade, `grad mu` ali e' da ordem de `9/h` e nao e' correcao pequena.
  //
  // O QUE CUSTOU, medido no benchmark de Hysing caso 1, nivel 5: circularidade
  // minima 0,8361 contra 0,9012 +- 0,0001 das tres referencias publicadas
  // (650 vezes a dispersao entre elas), velocidade de ascensao +5%, area +2,5%.
  // E o excesso de velocidade NAO caia com refino -- nivel 6 ainda +4,3%.
  // PATAMAR SOB REFINO E' ASSINATURA DE TERMO FALTANDO, nao de malha grossa:
  // segunda ordem teria dividido o erro por quatro.
  //
  // COM `mu` UNIFORME O TERMO EXTRA E' IDENTICAMENTE NULO (`grad mu = 0`), e
  // isto volta a ser o laplaciano sobre `Re` bit a bit.  Os casos monofasicos
  // -- canal e caixa manufaturada -- continuam verificando este codigo
  // exatamente, e e' por isso que a mudanca pode entrar sem reabrir aquilo.
  return tensao_viscosa_fluxo(d, fid) + tensao_viscosa_transposto(d, fid);
}

double NavierStokes::tensao_viscosa_fluxo(int d, long fid) const
{
  // SO' A PARTE DE FLUXO.  A soma sai na MESMA ORDEM de antes da separacao, e
  // isso e' deliberado: os casos monofasicos tem de reproduzir bit a bit.
  const double aqui = u_.em(d, fid);
  const Ponto2 xf = fd_->centro(d, fid);
  double total = 0.0;

  for (int e = 0; e < 2; e++) {
    double fluxo[2], dist[2];
    for (int lado = 0; lado < 2; lado++) {
      const int sgn = (lado == 0) ? -1 : +1;
      const Vizinha v = vizinha_de(d, fid, e, sgn);
      const double uv = v.dentro ? v.valor : contorno_vizinho(d, fid, e, sgn);
      const double dd = v.dentro ? v.dist : dist_ate_parede(d, fid, e, sgn);
      if (dd <= 0.0) { fluxo[lado] = 0.0; dist[lado] = 1.0; continue; }

      // `mu` NO MEIO DO CAMINHO, e por duas regras DIFERENTES conforme a
      // direcao -- e' assim no HiGFlow, e a diferenca nao e' detalhe.
      //
      //   `e == d`: a posicao e' um CENTRO DE CELULA, e `mu` e' o valor dela.
      //   Nao ha' media a fazer: o ponto E' onde a propriedade vive;
      //
      //   `e != d`: a posicao e' uma QUINA, e `mu` e' a media HARMONICA das
      //   QUATRO celulas que a cercam.  Harmonica, e nao aritmetica: e' ela
      //   que preserva a continuidade da TENSAO CISALHANTE atraves de uma
      //   interface plana.  Com a aritmetica a tensao salta na interface, e
      //   isso nao falha -- muda a taxa de deformacao exatamente onde os dois
      //   fluidos se tocam, que e' o que o caso bifasico mede.
      //
      // Tomar o valor de UMA celula na quina -- o que a regra semiaberta de
      // `id_do_ponto` devolveria -- seria escolher arbitrariamente uma das
      // quatro, e com razao 10 de viscosidade a escolha vale uma ordem de
      // grandeza.
      Ponto2 xm = xf;
      xm[e] += sgn * dd / 2.0;
      const double mu_m = (e == d) ? mu_em(xm) : mu_na_quina(xm);
      // O FLUXO ENTRA SEM O SINAL DO LADO.  A divergencia e'
      // `(fluxo_+ - fluxo_-)/D`, e com `fluxo_s = mu (u_s - u)/d_s` os dois
      // ja' saem com o sinal certo -- o de fora aponta para fora nos dois
      // lados.  Multiplicar pelo sinal do lado soma onde devia subtrair, e o
      // laplaciano vira o oposto do que e': medido, a caixa manufaturada foi
      // de 9,97e-02 para 1,15e+01.
      fluxo[lado] = mu_m * (uv - aqui) / dd;
      dist[lado] = dd;
    }
    total += (fluxo[0] + fluxo[1]) / ((dist[0] + dist[1]) / 2.0);
  }
  return total;
}


double NavierStokes::tensao_viscosa_transposto(int d, long fid) const
{
  // O TRANSPOSTO QUE SOBREVIVE.  Recalcula as vizinhancas em vez de receber as
  // do fluxo -- custa uma varredura e evita que as duas partes fiquem
  // amarradas uma na outra, que e' o que impediria o fluxo de ir para a
  // matriz do semi-implicito.
  const Ponto2 xf = fd_->centro(d, fid);
  double u_viz[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
  double d_viz[2][2] = {{1.0, 1.0}, {1.0, 1.0}};
  double mu_viz[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
  bool dir_ok[2] = {true, true};

  for (int e = 0; e < 2; e++)
    for (int lado = 0; lado < 2; lado++) {
      const int sgn = (lado == 0) ? -1 : +1;
      const Vizinha v = vizinha_de(d, fid, e, sgn);
      const double uv = v.dentro ? v.valor : contorno_vizinho(d, fid, e, sgn);
      const double dd = v.dentro ? v.dist : dist_ate_parede(d, fid, e, sgn);
      if (dd <= 0.0) { dir_ok[e] = false; continue; }
      Ponto2 xm = xf;
      xm[e] += sgn * dd / 2.0;
      u_viz[e][lado] = uv;
      d_viz[e][lado] = dd;
      mu_viz[e][lado] = (e == d) ? mu_em(xm) : mu_na_quina(xm);
    }

  if (!dir_ok[d]) return 0.0;
  const double hd = (d_viz[d][0] + d_viz[d][1]) / 2.0;
  double extra = 0.0;
  for (int e = 0; e < 2; e++) {
    if (!dir_ok[e]) continue;
    const double dmu =
        (mu_viz[e][1] - mu_viz[e][0]) / ((d_viz[e][0] + d_viz[e][1]) / 2.0);
    if (dmu == 0.0) continue;

    double dujdxi;
    if (e == d) {
      dujdxi = (u_viz[d][1] - u_viz[d][0]) / (d_viz[d][0] + d_viz[d][1]);
    } else {
      Ponto2 xa = xf, xb = xf;
      xa[d] -= hd / 2.0;
      xb[d] += hd / 2.0;
      bool ia = false, ib = false;
      const double ua = u_no_ponto(e, xa, &ia);
      const double ub = u_no_ponto(e, xb, &ib);
      if (!ia || !ib) continue;
      dujdxi = (ub - ua) / hd;
    }
    extra += dujdxi * dmu;
  }
  return extra;
}

}  // namespace t8flow2d
