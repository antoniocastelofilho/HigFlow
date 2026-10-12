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

  // A MALHA ADAPTADA E' RECUSADA, e nao tolerada.
  //
  // A transicao 2:1 numa faceta de VELOCIDADE nao e' o mesmo problema que numa
  // faceta de fluxo escalar: a componente normal se restringe por conservacao
  // de vazao e a tangencial por interpolacao, e usar a mesma regra para as duas
  // produz divergencia que nao fecha -- um defeito que aparece como o solver de
  // pressao nao convergindo, longe da causa.
  //
  // Rodar e dar numero seria pior que recusar: o numero sairia plausivel.
  if (m.num_folhas_em_transicao() > 0)
    throw std::invalid_argument(
        "NavierStokes: a malha tem " +
        std::to_string(m.num_folhas_em_transicao()) +
        " folhas em transicao 2:1, e este solver ainda nao a trata.  "
        "A restricao da componente normal e a da tangencial sao diferentes, e "
        "aplicar a mesma as duas produz divergencia que nao fecha.");

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
  for (long b = 0; b < m.num_blocos(); b++)
    for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
      const Ponto2 hh = m.lados_da_folha(b, i);
      for (int d = 0; d < 2; d++)
        if (std::fabs(hh[d] - h_[d]) > 1e-12 * h_[d])
          throw std::invalid_argument(
              "NavierStokes: a folha (" + std::to_string(b) + "," +
              std::to_string(i) + ") tem lado " + std::to_string(hh[d]) +
              " na direcao " + std::to_string(d) + " contra " +
              std::to_string(h_[d]) + " da primeira.  Este solver exige malha "
              "uniforme -- blocos com niveis diferentes tambem nao servem.");
    }

  // O INDICE DE FACETA POR PONTO.  Construido uma vez; a alternativa seria
  // varrer as facetas a cada consulta, e a montagem do passo consulta uma vez
  // por termo por faceta -- o mesmo quadratico que ja' custou caro no Poisson.
  for (int d = 0; d < 2; d++) {
    fac_por_ponto_[d].reserve(static_cast<size_t>(fd.num_facetas(d)) * 2);
    for (long f = 0; f < fd.num_facetas(d); f++) {
      const unsigned long long ch = chave_de_ponto(fd.centro(d, f), h_);
      const auto it = fac_por_ponto_[d].find(ch);
      if (it != fac_por_ponto_[d].end())
        throw std::runtime_error(
            "NavierStokes: duas facetas da direcao " + std::to_string(d) +
            " caem na mesma chave de indice -- a quantizacao nao serve para "
            "esta malha");
      fac_por_ponto_[d][ch] = f;
    }
  }

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
  const auto it = fac_por_ponto_[d].find(chave_de_ponto(x, h_));
  return (it == fac_por_ponto_[d].end()) ? -1 : it->second;
}

double NavierStokes::u_vizinha(int d, const Ponto2 &xf, int e, int s,
                               bool *achou) const
{
  const long aqui_f = faceta_em(d, xf);
  if (aqui_f < 0)
    throw std::runtime_error(
        "NavierStokes::u_vizinha: nao ha' faceta da direcao " +
        std::to_string(d) + " em (" + std::to_string(xf[0]) + ", " +
        std::to_string(xf[1]) + ").  Esta funcao fecha o contorno supondo que "
        "a parede esta' a MEIO PASSO, e isso so' vale partindo de uma faceta "
        "que existe.");

  Ponto2 xn = xf;
  xn[e] += s * h_[e];
  const long f = faceta_em(d, xn);
  if (f >= 0) {
    if (achou) *achou = true;
    return u_.em(d, f);
  }
  if (achou) *achou = false;

  // FORA DO DOMINIO: a parede esta' no MEIO DO CAMINHO, e e' la' que a
  // condicao vive.  Tomar a condicao no ponto fantasma em vez de na parede
  // seria perguntar a uma regiao que nao o cobre -- e, com `padrao` instalado,
  // obter resposta errada em silencio.
  Ponto2 xp = xf;
  xp[e] += s * h_[e] / 2.0;
  const Contorno &c = bcu_[d]->em(xp);
  const double dentro = u_.em(d, aqui_f);
  switch (c.tipo) {
    case TipoContorno::Dirichlet:
      // REFLEXAO: a media entre o de dentro e o fantasma e' o valor na parede.
      return 2.0 * c.g(xp) - dentro;
    case TipoContorno::Neumann:
      return dentro + s * h_[e] * c.g(xp);
    case TipoContorno::Robin: {
      double A = 0.0, B = 0.0;
      c.afim(xp, h_[e], &A, &B);
      return A * dentro + B;
    }
  }
  return dentro;
}

double NavierStokes::velocidade_que_advecta(int d, const Ponto2 &xf, int e,
                                            int s) const
{
  // A VELOCIDADE NA POSICAO DO FLUXO, que fica meio passo adiante na direcao
  // `e`.  Para `e == d` essa posicao e' um CENTRO DE CELULA e a velocidade ali
  // e' a media das duas facetas vizinhas; para `e != d` e' uma QUINA, e a
  // velocidade que atravessa e' a componente `e`, media das duas facetas de
  // direcao `e` que a cercam.
  //
  // Confundir as duas e' o erro classico do arranjo deslocado, e ele nao falha:
  // produz um termo convectivo plausivel e errado.
  if (e == d) {
    bool achou = false;
    const double adiante = u_vizinha(d, xf, d, s, &achou);
    return 0.5 * (u_.em(d, faceta_em(d, xf)) + adiante);
  }

  // Quina: as duas facetas de direcao `e` adjacentes, uma de cada lado em `d`.
  Ponto2 a = xf, b = xf;
  a[e] += s * h_[e] / 2.0;
  a[d] -= h_[d] / 2.0;
  b[e] += s * h_[e] / 2.0;
  b[d] += h_[d] / 2.0;
  const long fa = faceta_em(e, a), fb = faceta_em(e, b);
  if (fa >= 0 && fb >= 0) return 0.5 * (u_.em(e, fa) + u_.em(e, fb));
  if (fa >= 0) return u_.em(e, fa);
  if (fb >= 0) return u_.em(e, fb);
  return 0.0;
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

double NavierStokes::valor_na_face(int d, const Ponto2 &xf, int e, int s,
                                   double V) const
{
  // QUEM E' `C`, `D` E `U` DEPENDE DO SINAL DE `V`, e nao da geometria: o
  // esquema olha para montante, e montante troca de lado quando o escoamento
  // troca.  Fixar `C` como "o daqui" e' o erro que faz o upwind virar
  // downwind na metade do dominio, e ele nao falha -- difunde.
  bool tem_d = false, tem_u = false;
  const double aqui = u_.em(d, faceta_em(d, xf));
  const double alem = u_vizinha(d, xf, e, s, &tem_d);

  double phiC, phiD, phiU;
  if (V >= 0.0) {
    phiC = aqui;
    phiD = alem;
    phiU = u_vizinha(d, xf, e, -s, &tem_u);
  } else {
    phiC = alem;
    phiD = aqui;
    // O ponto a montante agora esta' DOIS passos adiante, e aqui o fechamento
    // de contorno NAO vale.
    //
    // `u_vizinha` supoe que a parede esta' a meio passo -- o que e' verdade
    // para o PRIMEIRO passo para fora, e falso para o segundo: dali a parede
    // ja' ficou para tras, e o ponto medio cai fora do dominio.  Chamando-a
    // mesmo assim, o contorno e' consultado num ponto que regiao nenhuma
    // cobre; com `padrao` instalado isso nao lanca, devolve a condicao errada
    // em silencio.
    //
    // Entao a busca aqui e' ESTRITA: ou a faceta existe, ou o esquema de
    // segunda ordem nao tem como ser montado e cai para upwind, contado.
    Ponto2 x2 = xf, x3 = xf;
    x2[e] += s * h_[e];
    x3[e] += 2.0 * s * h_[e];
    const long f3 = faceta_em(d, x3);
    if (f3 >= 0) {
      phiU = u_.em(d, f3);
      tem_u = true;
    } else {
      tem_u = false;
      phiU = phiC;
    }
  }

  if (conv_ == Convectivo::Upwind1) return phiC;
  if (conv_ == Convectivo::Central) return 0.5 * (phiC + phiD);

  // SEM O PONTO A MONTANTE o esquema de segunda ordem nao existe, e a queda
  // para upwind e' CONTADA.  Perto da parede isso e' inevitavel; no meio do
  // dominio seria sintoma.
  if (!tem_u) {
    rebaixos_conv_++;
    return phiC;
  }

  const double den = phiD - phiU;
  if (std::fabs(den) < EPS) return phiC;        // campo localmente plano
  const double fc = (phiC - phiU) / den;
  return phiU + nvd(conv_, fc) * den;
}

double NavierStokes::conveccao(int d, const Ponto2 &xf) const
{
  // FORMA DIVERGENTE, `sum_e d(u_e u_d)/dx_e`, e nao a advectiva `u_e du_d/dx_e`.
  //
  // As duas coincidem no continuo porque a divergencia e' nula, e NAO coincidem
  // no discreto: a divergente conserva quantidade de movimento por construcao,
  // e a advectiva so' a conserva quando a divergencia discreta ja' e' nula --
  // o que no meio do passo ela nao e'.
  double soma = 0.0;
  for (int e = 0; e < 2; e++) {
    const double Vp = velocidade_que_advecta(d, xf, e, +1);
    const double Vm = velocidade_que_advecta(d, xf, e, -1);
    const double fp = valor_na_face(d, xf, e, +1, Vp);
    const double fm = valor_na_face(d, xf, e, -1, Vm);
    soma += (Vp * fp - Vm * fm) / h_[e];
  }
  return soma;
}

double NavierStokes::difusao(int d, const Ponto2 &xf) const
{
  const double aqui = u_.em(d, faceta_em(d, xf));
  double lap = 0.0;
  for (int e = 0; e < 2; e++) {
    const double mais = u_vizinha(d, xf, e, +1, nullptr);
    const double menos = u_vizinha(d, xf, e, -1, nullptr);
    lap += (mais - 2.0 * aqui + menos) / (h_[e] * h_[e]);
  }
  return lap / re_;
}

double NavierStokes::grad_pressao(int d, const Ponto2 &xf) const
{
  // A PRESSAO VIVE NOS CENTROS e a faceta esta' EXATAMENTE entre dois deles --
  // e' isso que o arranjo deslocado compra: o gradiente na faceta e' uma
  // diferenca de dois pontos, de segunda ordem, sem interpolacao e sem o modo
  // de tabuleiro que o arranjo colocado admite.
  Ponto2 a = xf, b = xf;
  a[d] -= h_[d] / 2.0;
  b[d] += h_[d] / 2.0;
  const long ca = m_->id_do_ponto(a), cb = m_->id_do_ponto(b);
  if (ca < 0 || cb < 0) return 0.0;   // faceta de contorno: nao se resolve ali
  return (p_.em(cb) - p_.em(ca)) / h_[d];
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
      Ponto2 a = x, b = x;
      a[d] -= h_[d] / 2.0;
      b[d] += h_[d] / 2.0;
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
        Ponto2 xi = x;
        xi[d] += s * h_[d];
        const long fi = faceta_em(d, xi);
        if (fi < 0) continue;                   // dominio de uma celula so'
        u_.em(d, f) = u_.em(d, fi) - s * h_[d] * c.g(x);
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
      a[d] -= h_[d] / 2.0;
      b[d] += h_[d] / 2.0;
      if (m_->id_do_ponto(a) < 0 || m_->id_do_ponto(b) < 0) continue;

      double rhs = -conveccao(d, x) + difusao(d, x);

      // A PRESSAO ANTIGA ENTRA SO' NA INCREMENTAL, e e' nisso que as duas
      // diferem.  Sem ela o preditor ignora o gradiente que ja' se conhece, e
      // a projecao tem de reconstrui-lo inteiro a cada passo -- o que custa
      // uma ordem na pressao.
      if (proj_ == Projecao::Incremental) rhs -= grad_pressao(d, x);

      if (fonte_faceta_) rhs += fonte_faceta_(x, d);

      ustar_.em(d, f) = u_.em(d, f) + dt_ * rhs;
    }
}

void NavierStokes::resolve_pressao()
{
  // O POISSON DA PROJECAO: `lap(dp) = div(u*) / dt`.
  //
  // A CONDICAO DELE NAO E' A DA PRESSAO, e' DERIVADA dela, e confundir as duas
  // e' erro de formulacao e nao de implementacao:
  //
  //   onde a pressao tem NEUMANN, o incremento tem Neumann HOMOGENEO -- a
  //   componente normal da velocidade ja' satisfaz a condicao depois do
  //   preditor, e a projecao nao pode altera-la;
  //
  //   onde a pressao tem DIRICHLET, o incremento tem Dirichlet do VALOR na
  //   nao incremental (ali `p = dp`) e HOMOGENEO na incremental (ali o valor
  //   ja' esta' em `p^n`).  Na incremental isso exige que `p^0` satisfaca a
  //   condicao -- com `p^0 = 0` e saida em zero, satisfaz.
  const long n = m_->num_folhas();
  const long p0 = m_->primeiro_id();
  t8flow::Solver s(t8flow::Motor::Petsc, static_cast<int>(p0), static_cast<std::size_t>(n), 5);
  s.tolerancia(1e-12);

  for (long b = 0; b < m_->num_blocos(); b++)
    for (long i = 0; i < m_->num_folhas_do_bloco(b); i++) {
      const long id = m_->id_global(b, i);
      const Ponto2 c = m_->centroide(b, i);
      const Ponto2 hh = m_->lados_da_folha(b, i);
      const double vol = hh[0] * hh[1];

      std::map<long, double> linha;
      double diag = 0.0, rhs = 0.0;

      for (int d = 0; d < 2; d++) {
        const double area = hh[(d + 1) % 2];
        for (int s2 = -1; s2 <= 1; s2 += 2) {
          Ponto2 xn = c;
          xn[d] += s2 * hh[d];
          const long viz = m_->id_do_ponto(xn);

          Ponto2 xp = c;
          xp[d] += s2 * hh[d] / 2.0;          // a parede/faceta

          // A DIVERGENCIA DE `u*`, montada na MESMA varredura: a faceta e' a
          // mesma, e percorre-la duas vezes abriria espaco para as duas
          // contas discordarem na borda.
          const long fid = faceta_em(d, xp);
          if (fid >= 0) rhs += s2 * ustar_.em(d, fid) * area / dt_;

          if (viz >= 0) {
            linha[viz] += area / hh[d];
            diag -= area / hh[d];
            continue;
          }

          const Contorno &cp = bcp_->em(xp);
          if (cp.tipo == TipoContorno::Dirichlet) {
            // O VALOR DO INCREMENTO NA PAREDE NAO E' ZERO NA INCREMENTAL, e
            // supor que fosse custou um nivel de pressao inteiro.
            //
            // A condicao e' `p^{n+1} = g`, logo `dp_parede = g - p^n_parede`.
            // Com Dirichlet HOMOGENEO no incremento, `dp = 0` resolve a
            // equacao do ponto fixo seja qual for o nivel de `p` -- a
            // condicao nao prende `p`, prende so' a VARIACAO dele.  Medido no
            // canal: o primeiro passo, de partida impulsiva, deixava `p` em
            // 62,99 perto da saida e o nivel ficava la' para sempre, com o
            // GRADIENTE certo (22,74 contra 23,25 de faixa) e o campo inteiro
            // deslocado.  A nao incremental nao sofria disso porque recalcula
            // `p` do zero a cada passo.
            //
            // `p^n` na parede vem de EXTRAPOLACAO LINEAR das duas celulas mais
            // proximas -- a parede esta' a `h/2` da primeira e a `3h/2` da
            // segunda, logo `p_parede = (3 p_c - p_cc)/2`.  Extrapolar a ordem
            // zero (`p_parede = p_c`) forcaria `p_c -> g` e introduziria um
            // erro de `h/2` vezes o gradiente -- exatamente o tipo de erro que
            // parece pequeno e e' de primeira ordem.
            double val = cp.g(xp);
            if (proj_ == Projecao::Incremental) {
              Ponto2 xcc = c;
              xcc[d] -= s2 * hh[d];
              const long cc = m_->id_do_ponto(xcc);
              const double p_parede = (cc >= 0)
                                          ? 1.5 * p_.em(id) - 0.5 * p_.em(cc)
                                          : p_.em(id);
              val -= p_parede;
            }
            diag -= 2.0 * area / hh[d];
            rhs -= 2.0 * area * val / hh[d];
          }
          // NEUMANN HOMOGENEO: fluxo nulo, nada a somar -- e e' por isso que
          // nao ha' ramo.  Escrever `+= 0` convidaria alguem a "consertar".
        }
      }

      if (fonte_celula_) rhs -= fonte_celula_(c) * vol;

      linha[id] += diag;
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
      a[d] -= h_[d] / 2.0;
      b[d] += h_[d] / 2.0;
      const long ca = m_->id_do_ponto(a), cb = m_->id_do_ponto(b);
      if (ca < 0 || cb < 0) continue;

      const double antes = u_.em(d, f);
      u_.em(d, f) = ustar_.em(d, f) - dt_ * (dp_.em(cb) - dp_.em(ca)) / h_[d];
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

  aplica_contorno_velocidade();
  preditor();
  // O CONTORNO VALE PARA `u*` TAMBEM: o Poisson usa a divergencia de `u*`, e
  // se a componente normal na parede nao for a da condicao, a compatibilidade
  // do Poisson quebra -- com Neumann em toda a volta ele fica insoluvel, e com
  // Dirichlet em algum lugar ele resolve e da' a resposta errada.
  for (int d = 0; d < 2; d++)
    for (long f = 0; f < fd_->num_facetas(d); f++) {
      const Ponto2 x = fd_->centro(d, f);
      Ponto2 a = x, b = x;
      a[d] -= h_[d] / 2.0;
      b[d] += h_[d] / 2.0;
      if (m_->id_do_ponto(a) >= 0 && m_->id_do_ponto(b) >= 0) continue;
      ustar_.em(d, f) = u_.em(d, f);
    }

  resolve_pressao();
  projeta();
  aplica_contorno_velocidade();

  t_ += dt_;
  n_passos_++;
}

double NavierStokes::divergencia_maxima() const
{
  double maxv = 0.0;
  for (long b = 0; b < m_->num_blocos(); b++)
    for (long i = 0; i < m_->num_folhas_do_bloco(b); i++) {
      const Ponto2 c = m_->centroide(b, i);
      const Ponto2 hh = m_->lados_da_folha(b, i);
      double div = 0.0;
      for (int d = 0; d < 2; d++)
        for (int s2 = -1; s2 <= 1; s2 += 2) {
          Ponto2 xp = c;
          xp[d] += s2 * hh[d] / 2.0;
          const long fid = faceta_em(d, xp);
          if (fid >= 0) div += s2 * u_.em(d, fid) / hh[d];
        }
      maxv = std::max(maxv, std::fabs(div));
    }
  return maxv;
}

}  // namespace t8flow2d
