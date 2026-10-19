#include "navierstokes3d.hxx"

#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "solverlin.hxx"

namespace t8flow3d {

namespace {

constexpr double EPS = 1e-14;

// A FUNCAO NVD, em variavel normalizada.  `fc` e'
// `(phi_C - phi_U)/(phi_D - phi_U)` e devolve `phi_f` normalizado.
//
// FORA DE [0,1] TODOS CAEM PARA UPWIND, e isso nao e' detalhe: e' o que torna
// o esquema limitado.  Um "segunda ordem" que nao cai e' central com outro
// nome, e oscila onde o gradiente e' forte.
//
// IDENTICA A' DO LADO 2D, e de proposito: o esquema NVD nao tem nada de
// dimensional -- ele age sobre tres valores ao longo de UMA linha.  O que
// muda em tres dimensoes e' quantas linhas existem por faceta, nao a funcao.
double nvd(Convectivo c, double fc)
{
  if (fc <= 0.0 || fc >= 1.0) return fc;        // upwind: limitado
  switch (c) {
    case Convectivo::Upwind1:
      return fc;
    case Convectivo::Central:
      return 0.5 * (fc + 1.0);
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

Convectivo convectivo_do_nome(const std::string &s)
{
  if (s == "central") return Convectivo::Central;
  if (s == "first_order" || s == "upwind1") return Convectivo::Upwind1;
  if (s == "cubista") return Convectivo::Cubista;
  if (s == "quick") return Convectivo::Quick;
  throw std::invalid_argument(
      "convectivo_do_nome: \"" + s + "\" nao e' valido.  Validos: central, "
      "first_order (upwind1), cubista, quick");
}

NavierStokes::NavierStokes(const Malha3 &m, const DominioFacetas &fd,
                           const ContornoPorRegiao &bc_u0,
                           const ContornoPorRegiao &bc_u1,
                           const ContornoPorRegiao &bc_u2,
                           const ContornoPorRegiao &bc_p)
    : m_(&m), fd_(&fd), u_(fd), p_(m), ustar_(fd), dp_(m)
{
  bcu_[0] = &bc_u0;
  bcu_[1] = &bc_u1;
  bcu_[2] = &bc_u2;
  bcp_ = &bc_p;

  // AS DUAS RECUSAS, e nenhuma produz numero.
  //
  // Recusar na CONSTRUCAO, e nao no primeiro passo: quem monta uma malha
  // adaptada e so' descobre no passo mil que ela nao serve ja' gastou a
  // montagem, e pior, pode ter visto numeros no caminho.
  if (m.num_processos() > 1)
    throw std::invalid_argument(
        "NavierStokes3D: este solver e' SERIAL (" +
        std::to_string(m.num_processos()) +
        " processos).  Faltam a troca de franja dos campos -- o estencil le' "
        "celula e faceta do vizinho, e os campos cobrem so' os ids locais -- e "
        "o percurso por blocos, que supoe todas as arvores locais.  Com a "
        "franja ausente o estencil leria lixo perto da fronteira de particao e "
        "a corrida daria numero");

  const long trans = m.num_folhas_em_transicao();
  if (trans > 0)
    throw std::invalid_argument(
        "NavierStokes3D: a malha tem " + std::to_string(trans) +
        " folhas em transicao 2:1, e este solver ainda nao a trata.  Em tres "
        "dimensoes a face da celula grossa e' coberta por QUATRO facetas "
        "finas, nao duas, e cada faceta tem DUAS direcoes transversais, nao "
        "uma -- as regras de 2D aplicadas aqui produzem divergencia que nao "
        "fecha, e o sintoma aparece como o solver de pressao nao convergindo, "
        "longe da causa");

  // O MAPA de id global para `(bloco, folha)`.  O balanco de fluxos precisa
  // dele porque a faceta de um lado pode pertencer ao VIZINHO, e a' celula so'
  // sobra pedi-la a ele -- ver `divergencia_em`.
  for (long b = 0; b < m.num_blocos(); b++)
    for (long i = 0; i < m.num_folhas_do_bloco(b); i++)
      bi_por_id_[m.id_global(b, i)] = {b, i};

  // AS DUAS CELULAS DE CADA FACETA, construidas uma vez.
  for (int d = 0; d < 3; d++)
    celulas_da_faceta_[d].assign(
        static_cast<size_t>(fd.num_facetas(d)), {{{-1, -1}, {-1, -1}}});
  for (long b = 0; b < m.num_blocos(); b++)
    for (long i = 0; i < m.num_folhas_do_bloco(b); i++)
      for (int d = 0; d < 3; d++)
        for (int lado = 0; lado < 2; lado++) {
          const long f = faceta_da_face(b, i, d, lado);
          if (f < 0) continue;
          // A celula esta' do lado OPOSTO ao `lado` que olha para a faceta:
          // se a faceta e' a face "mais" da celula, a celula esta' do lado
          // "menos" DELA.
          celulas_da_faceta_[d][static_cast<size_t>(f)][1 - lado] = {b, i};
        }

  // A PRESSAO TODA-NEUMANN E' RECUSADA, pela mesma razao do lado 2D: a
  // constante fica no nucleo e o PETSc converge para um elemento qualquer do
  // subespaco -- isso DA' NUMERO, e o numero e' plausivel.
  for (long b = 0; b < m.num_blocos() && !p_tem_dirichlet_; b++)
    for (long i = 0; i < m.num_folhas_do_bloco(b) && !p_tem_dirichlet_; i++) {
      const Ponto3 c = m.centroide(b, i);
      const Ponto3 hh = m.lados_da_folha(b, i);
      for (int d = 0; d < 3 && !p_tem_dirichlet_; d++)
        for (int lado = 0; lado < 2; lado++) {
          if (!m.ids_vizinhos(b, i, 2 * d + lado).empty()) continue;
          Ponto3 xp = c;
          xp[d] += (lado == 0 ? -1.0 : 1.0) * hh[d] / 2.0;
          if (bc_p.em(xp).tipo == TipoContorno::Dirichlet) {
            p_tem_dirichlet_ = true;
            break;
          }
        }
    }
  if (!p_tem_dirichlet_)
    throw std::invalid_argument(
        "NavierStokes3D: a pressao tem Neumann em TODA a parede, e o Poisson "
        "da projecao fica singular -- a constante esta' no nucleo.  Imponha "
        "Dirichlet em pelo menos um ponto de parede; numa solucao manufaturada "
        "o valor ali e' conhecido");

  // O MAPA DE POSICAO PARA FACETA.  O quantum e' um quarto do menor lado de
  // celula: fino o bastante para nao colidir centros distintos, grosso o
  // bastante para absorver o arredondamento de somar meios-lados.
  double hmin = 1e300;
  for (long b = 0; b < m.num_blocos(); b++)
    for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
      const Ponto3 hh = m.lados_da_folha(b, i);
      for (int d = 0; d < 3; d++) hmin = std::min(hmin, hh[d]);
    }
  quantum_ = hmin / 4.0;
  for (int d = 0; d < 3; d++)
    for (long f = 0; f < fd.num_facetas(d); f++)
      faceta_por_ponto_[d][chave(fd.centro(d, f))] = f;

  u_.preenche(0.0);
  p_.preenche(0.0);
  dp_.preenche(0.0);
}

void NavierStokes::reynolds(double re)
{
  if (!(re > 0.0))
    throw std::invalid_argument("NavierStokes3D: Reynolds nao positivo");
  re_ = re;
}

long long NavierStokes::chave(const Ponto3 &x) const
{
  // Tres coordenadas quantizadas num inteiro so'.  O deslocamento de 1<<20
  // cobre coordenadas negativas e dominios de ate' ~2^20 quanta por direcao.
  long long k = 0;
  for (int d = 0; d < 3; d++) {
    const long long q =
        static_cast<long long>(std::llround(x[d] / quantum_)) + (1LL << 20);
    k = k * (1LL << 21) + q;
  }
  return k;
}

long NavierStokes::faceta_em(int d, const Ponto3 &x) const
{
  const auto it = faceta_por_ponto_[d].find(chave(x));
  return (it == faceta_por_ponto_[d].end()) ? -1 : it->second;
}

double NavierStokes::mu_em(const Ponto3 &x) const
{
  if (mu_ == nullptr) return 1.0 / re_;
  const long id = m_->id_do_ponto(x);
  if (id < 0) return 1.0 / re_;
  return mu_->em(id);
}

double NavierStokes::mu_na_aresta(const Ponto3 &xa, int d, int e) const
{
  // AS QUATRO CELULAS QUE TOCAM A ARESTA: o ponto esta' sobre os planos de
  // face de `d` e de `e`, e num centro de celula na terceira direcao.  As
  // quatro saem de cutucar meio quantum para cada lado em `d` e em `e`.
  if (mu_ == nullptr) return 1.0 / re_;
  const double cut = quantum_;
  double soma_inv = 0.0;
  int n = 0;
  for (int sd = -1; sd <= 1; sd += 2)
    for (int se = -1; se <= 1; se += 2) {
      Ponto3 p = xa;
      p[d] += sd * cut;
      p[e] += se * cut;
      const long id = m_->id_do_ponto(p);
      if (id < 0) continue;              // fora do dominio: a aresta e' de parede
      const double v = mu_->em(id);
      if (!(v > 0.0)) continue;
      soma_inv += 1.0 / v;
      n++;
    }
  if (n == 0) return 1.0 / re_;
  return n / soma_inv;
}

double NavierStokes::u_no_ponto(int d, const Ponto3 &x, bool *dentro) const
{
  // CAMINHO CURTO: ha' faceta exatamente aqui.
  const long f = faceta_em(d, x);
  if (f >= 0) {
    if (dentro) *dentro = true;
    return u_.em(d, f);
  }

  // INTERPOLACAO TRILINEAR SOBRE A REDE DE FACETAS DA DIRECAO `d`.
  //
  // A PRIMEIRA VERSAO INTERPOLAVA SO' AO LONGO DE `d`, e isso custou PRIMEIRA
  // ORDEM no termo convectivo -- medido: 0,944 e 0,988 entre os niveis 3-4-5,
  // com diferencas centrais, que deveriam dar segunda.  O viscoso, que nao
  // usa este caminho, dava 1,88 e 2,00 nas mesmas malhas, e foi esse contraste
  // que localizou o defeito.
  //
  // A RAZAO: a velocidade que advecta atraves da face transversal do volume de
  // controle e' avaliada numa ARESTA -- deslocada meia celula em `d` em
  // relacao a's facetas de direcao `e` que a carregam.  Interpolar so' em `e`
  // devolve o valor no CENTRO da celula, nao na aresta, e o erro e' `O(h)`.
  //
  // E' A MESMA CLASSE DE DEFEITO QUE NO LADO 2D, onde a quina devolvia a
  // faceta de uma das duas celulas em vez da media, e o erro da caixa
  // manufaturada subiu 20% em `u` e 37% em `v`.  Em tres dimensoes ha' DUAS
  // direcoes transversais, logo duas interpolacoes a esquecer.
  //
  // A REDE: na direcao `d` as facetas ficam nos multiplos de `h` (as faces);
  // nas outras duas, nos meios (os centros de celula).
  const double h = quantum_ * 4.0;              // malha uniforme: so' ha' um
  int k0[3];
  double w[3];
  for (int j = 0; j < 3; j++) {
    const double t = (j == d) ? (x[j] / h) : (x[j] / h - 0.5);
    k0[j] = static_cast<int>(std::floor(t + 1e-9));
    w[j] = t - k0[j];
    if (w[j] < 1e-9) w[j] = 0.0;
    if (w[j] > 1.0 - 1e-9) { w[j] = 0.0; k0[j] += 1; }
  }

  // OS CANTOS QUE EXISTEM, com os pesos RENORMALIZADOS.  Perto da parede
  // algum canto cai fora; usar o fantasma ali seria outro esquema, e
  // renormalizar degrada suavemente -- que e' o que o convectivo ja' faz
  // naquela regiao, contando o rebaixamento.
  double soma = 0.0, peso = 0.0;
  for (int a = 0; a < 2; a++)
    for (int b = 0; b < 2; b++)
      for (int c = 0; c < 2; c++) {
        const double wt = (a ? w[0] : 1.0 - w[0]) * (b ? w[1] : 1.0 - w[1]) *
                          (c ? w[2] : 1.0 - w[2]);
        if (wt <= 0.0) continue;
        Ponto3 q;
        q[0] = ((d == 0) ? (k0[0] + a) : (k0[0] + a + 0.5)) * h;
        q[1] = ((d == 1) ? (k0[1] + b) : (k0[1] + b + 0.5)) * h;
        q[2] = ((d == 2) ? (k0[2] + c) : (k0[2] + c + 0.5)) * h;
        const long fq = faceta_em(d, q);
        if (fq < 0) continue;
        soma += wt * u_.em(d, fq);
        peso += wt;
      }

  if (peso <= 0.0) {
    if (dentro) *dentro = false;
    return 0.0;
  }
  if (dentro) *dentro = true;
  return soma / peso;
}

NavierStokes::Vizinha NavierStokes::vizinha_de(int d, long fid, int e,
                                               int s) const
{
  Vizinha v;
  const Ponto3 xf = fd_->centro(d, fid);
  const auto &cel = celulas_da_faceta_[d][static_cast<size_t>(fid)];

  if (e == d) {
    // AO LONGO DA NORMAL: a vizinha e' a outra face da celula daquele lado, e
    // a distancia e' o lado DELA.
    const auto &c = (s > 0) ? cel[1] : cel[0];
    if (c.first < 0) return v;
    const Ponto3 hc = m_->lados_da_folha(c.first, c.second);
    const long f = faceta_da_face(c.first, c.second, d, (s > 0) ? 1 : 0);
    if (f < 0) return v;
    v.valor = u_.em(d, f);
    v.dist = hc[d];
    v.dentro = true;
    return v;
  }

  // TRANSVERSALMENTE: a faceta de mesma direcao, deslocada em `e`.  Em malha
  // uniforme a distancia e' o lado da celula em `e`.
  const auto &c = (cel[0].first >= 0) ? cel[0] : cel[1];
  if (c.first < 0) return v;
  const Ponto3 hc = m_->lados_da_folha(c.first, c.second);
  Ponto3 alvo = xf;
  alvo[e] += s * hc[e];
  const long f = faceta_em(d, alvo);
  if (f < 0) return v;
  v.valor = u_.em(d, f);
  v.dist = hc[e];
  v.dentro = true;
  return v;
}

double NavierStokes::contorno_vizinho(int d, long fid, int e, int s) const
{
  // O FECHAMENTO, com a parede a MEIO LADO da celula.  Afim na propria
  // incognita, como no lado 2D.
  const Ponto3 xf = fd_->centro(d, fid);
  const auto &cel = celulas_da_faceta_[d][static_cast<size_t>(fid)];
  const auto &c = (cel[0].first >= 0) ? cel[0] : cel[1];
  const double aqui = u_.em(d, fid);
  if (c.first < 0 || e == d) return aqui;
  const Ponto3 hc = m_->lados_da_folha(c.first, c.second);

  Ponto3 xp = xf;
  xp[e] += s * hc[e] / 2.0;
  const Contorno &cc = bcu_[d]->em(xp);
  switch (cc.tipo) {
    case TipoContorno::Dirichlet: return 2.0 * cc.g(xp) - aqui;
    case TipoContorno::Neumann:   return aqui + s * hc[e] * cc.g(xp);
    case TipoContorno::Robin: {
      double A = 0.0, B = 0.0;
      cc.afim(xp, hc[e], &A, &B);
      return A * aqui + B;
    }
  }
  return aqui;
}

double NavierStokes::comprimento_do_volume(int d, long fid, int e) const
{
  // O VOLUME DE CONTROLE DA FACETA: meia celula de cada lado ao longo da
  // normal, e a celula inteira nas DUAS transversais.
  const auto &cel = celulas_da_faceta_[d][static_cast<size_t>(fid)];
  const auto &c = (cel[0].first >= 0) ? cel[0] : cel[1];
  if (c.first < 0) return 0.0;
  if (e != d) return m_->lados_da_folha(c.first, c.second)[e];

  double soma = 0.0;
  int n = 0;
  for (int lado = 0; lado < 2; lado++) {
    if (cel[lado].first < 0) continue;
    soma += m_->lados_da_folha(cel[lado].first, cel[lado].second)[d] / 2.0;
    n++;
  }
  // Na faceta de contorno so' ha' meia celula: o volume e' meio, e usar o
  // inteiro daria um fluxo diluido pela metade ali.
  return (n > 0) ? soma : 0.0;
}

double NavierStokes::velocidade_que_advecta(int d, long fid, int e,
                                            int s) const
{
  // A VELOCIDADE NA POSICAO DO FLUXO.  Para `e == d` essa posicao e' um CENTRO
  // DE CELULA e a componente que atravessa e' a propria `d`; para `e != d` e'
  // uma ARESTA, e a que atravessa e' a componente `e`.
  //
  // Confundir as duas e' o erro classico do arranjo deslocado, e ele nao
  // falha: produz um termo convectivo plausivel e errado.  Em tres dimensoes
  // ha' DUAS direcoes transversais, logo duas chances de confundir.
  const auto &cel = celulas_da_faceta_[d][static_cast<size_t>(fid)];
  const Ponto3 xf = fd_->centro(d, fid);
  bool ok = false;

  if (e == d) {
    const auto &c = (s > 0) ? cel[1] : cel[0];
    if (c.first < 0) return u_.em(d, fid);
    const double v = u_no_ponto(d, m_->centroide(c.first, c.second), &ok);
    return ok ? v : u_.em(d, fid);
  }

  const auto &c = (cel[0].first >= 0) ? cel[0] : cel[1];
  if (c.first < 0) return 0.0;
  const Ponto3 hc = m_->lados_da_folha(c.first, c.second);
  const Ponto3 cc = m_->centroide(c.first, c.second);
  Ponto3 aresta = xf;
  aresta[e] = cc[e] + s * hc[e] / 2.0;
  const double v = u_no_ponto(e, aresta, &ok);
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
    // fantasma do contorno, e usa-lo seria legitimo -- mas e' OUTRO esquema.
    // Fora do dominio, rebaixa e CONTA.
    const Vizinha vu = vizinha_de(d, fid, e, -s);
    tem_u = vu.dentro;
    phiU = vu.dentro ? vu.valor : phiC;
  } else {
    phiC = alem;
    phiD = aqui;
    // O PONTO A MONTANTE ESTA' DOIS PASSOS ADIANTE, e ali o fechamento de
    // contorno NAO vale: ele supoe a parede a meio caminho, e do segundo
    // passo em diante ela ja' ficou para tras.  A busca e' ESTRITA.
    const auto &cel = celulas_da_faceta_[d][static_cast<size_t>(fid)];
    const auto &c = (cel[0].first >= 0) ? cel[0] : cel[1];
    if (c.first < 0) { rebaixos_conv_++; return phiC; }
    Ponto3 x2 = fd_->centro(d, fid);
    x2[e] += 2.0 * s * m_->lados_da_folha(c.first, c.second)[e];
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
  //
  // TRES direcoes, e nao duas.
  double soma = 0.0;
  for (int e = 0; e < 3; e++) {
    const double Vp = velocidade_que_advecta(d, fid, e, +1);
    const double Vm = velocidade_que_advecta(d, fid, e, -1);
    const double fp = valor_na_face(d, fid, e, +1, Vp);
    const double fm = valor_na_face(d, fid, e, -1, Vm);
    const double L = comprimento_do_volume(d, fid, e);
    if (L > 0.0) soma += (Vp * fp - Vm * fm) / L;
  }
  return soma;
}

double NavierStokes::tensao_viscosa(int d, long fid) const
{
  return tensao_viscosa_fluxo(d, fid) + tensao_viscosa_transposto(d, fid);
}

double NavierStokes::tensao_viscosa_fluxo(int d, long fid) const
{
  // `div(mu grad u)` EM FORMA DE FLUXO, somando as TRES direcoes.  Em duas
  // dimensoes eram duas; a direcao a mais nao e' so' um termo a mais, e' uma
  // segunda direcao TRANSVERSAL -- e e' nas transversais que `mu` vem da
  // media harmonica sobre a aresta.
  const double aqui = u_.em(d, fid);
  const Ponto3 xf = fd_->centro(d, fid);
  double total = 0.0;

  for (int e = 0; e < 3; e++) {
    double fluxo[2] = {0.0, 0.0}, dist[2] = {1.0, 1.0};
    for (int lado = 0; lado < 2; lado++) {
      const int sgn = (lado == 0) ? -1 : +1;
      const Vizinha v = vizinha_de(d, fid, e, sgn);
      double uv, dd;
      if (v.dentro) {
        uv = v.valor;
        dd = v.dist;
      } else {
        uv = contorno_vizinho(d, fid, e, sgn);
        const auto &cel = celulas_da_faceta_[d][static_cast<size_t>(fid)];
        const auto &c = (cel[0].first >= 0) ? cel[0] : cel[1];
        if (c.first < 0) continue;
        dd = m_->lados_da_folha(c.first, c.second)[e];
      }
      if (dd <= 0.0) continue;

      Ponto3 xm = xf;
      xm[e] += sgn * dd / 2.0;
      // `e == d` cai num CENTRO de celula -- `mu` e' o valor dela.
      // `e != d` cai numa ARESTA -- media harmonica das quatro.
      const double mu_m = (e == d) ? mu_em(xm) : mu_na_aresta(xm, d, e);
      // O FLUXO ENTRA SEM O SINAL DO LADO: com `fluxo_s = mu (u_s - u)/d_s` os
      // dois ja' saem com o sinal certo, o de fora apontando para fora nos
      // dois lados.
      fluxo[lado] = mu_m * (uv - aqui) / dd;
      dist[lado] = dd;
    }
    total += (fluxo[0] + fluxo[1]) / ((dist[0] + dist[1]) / 2.0);
  }
  return total;
}

double NavierStokes::tensao_viscosa_transposto(int d, long fid) const
{
  // `(d_e mu)(d_d u_e)`, somado sobre as TRES direcoes.
  //
  // IDENTICAMENTE NULO COM `mu` UNIFORME, e e' isso que faz o caso monofasico
  // reproduzir o laplaciano bit a bit.
  if (mu_ == nullptr) return 0.0;
  const Ponto3 xf = fd_->centro(d, fid);

  // A distancia em `d`, para a derivada cruzada.
  const auto &cel = celulas_da_faceta_[d][static_cast<size_t>(fid)];
  if (cel[0].first < 0 || cel[1].first < 0) return 0.0;
  const double hd = (m_->lados_da_folha(cel[0].first, cel[0].second)[d] +
                     m_->lados_da_folha(cel[1].first, cel[1].second)[d]) / 2.0;

  double extra = 0.0;
  for (int e = 0; e < 3; e++) {
    double mu_lado[2] = {0.0, 0.0}, dist[2] = {1.0, 1.0};
    bool ok = true;
    for (int lado = 0; lado < 2; lado++) {
      const int sgn = (lado == 0) ? -1 : +1;
      const Vizinha v = vizinha_de(d, fid, e, sgn);
      const auto &c = (cel[0].first >= 0) ? cel[0] : cel[1];
      const double dd = v.dentro ? v.dist
                                 : m_->lados_da_folha(c.first, c.second)[e];
      if (dd <= 0.0) { ok = false; break; }
      Ponto3 xm = xf;
      xm[e] += sgn * dd / 2.0;
      mu_lado[lado] = (e == d) ? mu_em(xm) : mu_na_aresta(xm, d, e);
      dist[lado] = dd;
    }
    if (!ok) continue;
    const double dmu = (mu_lado[1] - mu_lado[0]) / ((dist[0] + dist[1]) / 2.0);
    if (dmu == 0.0) continue;

    double dudx;
    if (e == d) {
      const Vizinha vm = vizinha_de(d, fid, d, -1);
      const Vizinha vp = vizinha_de(d, fid, d, +1);
      if (!vm.dentro || !vp.dentro) continue;
      dudx = (vp.valor - vm.valor) / (vm.dist + vp.dist);
    } else {
      // `d_d u_e`: a componente CRUZADA, nos dois centros de celula que a
      // faceta separa.
      Ponto3 xa = xf, xb = xf;
      xa[d] -= hd / 2.0;
      xb[d] += hd / 2.0;
      bool ia = false, ib = false;
      const double ua = u_no_ponto(e, xa, &ia);
      const double ub = u_no_ponto(e, xb, &ib);
      if (!ia || !ib) continue;
      dudx = (ub - ua) / hd;
    }
    extra += dudx * dmu;
  }
  return extra;
}

void NavierStokes::passo_de_tempo(double dt)
{
  if (!(dt > 0.0))
    throw std::invalid_argument("NavierStokes3D: passo de tempo nao positivo");
  dt_ = dt;
}

double NavierStokes::gradiente_em(const CampoEscalar &campo, int d,
                                  long fid) const
{
  // A PRESSAO VIVE NOS CENTROS, e a faceta esta' entre dois deles.  Sem
  // transicao 2:1 -- recusada na construcao -- os dois centros estao alinhados
  // com a normal e a distancia e' a soma dos dois meios-lados.
  //
  // NA FACETA DE CONTORNO DEVOLVE ZERO, e nao uma extrapolacao: ali a
  // velocidade nao e' incognita da projecao, vem da condicao.  Devolver algo
  // faria a correcao mexer num valor que a condicao ja' fixou.
  const auto &cel = celulas_da_faceta_[d][static_cast<size_t>(fid)];
  if (cel[0].first < 0 || cel[1].first < 0) return 0.0;
  const double h0 = m_->lados_da_folha(cel[0].first, cel[0].second)[d];
  const double h1 = m_->lados_da_folha(cel[1].first, cel[1].second)[d];
  const double v0 = campo.em(m_->id_global(cel[0].first, cel[0].second));
  const double v1 = campo.em(m_->id_global(cel[1].first, cel[1].second));
  return (v1 - v0) / ((h0 + h1) / 2.0);
}

double NavierStokes::volume_dual_da_faceta(int d, long fid) const
{
  const auto &cel = celulas_da_faceta_[d][static_cast<size_t>(fid)];
  if (cel[0].first < 0 || cel[1].first < 0) return 0.0;
  const double h0 = m_->lados_da_folha(cel[0].first, cel[0].second)[d];
  const double h1 = m_->lados_da_folha(cel[1].first, cel[1].second)[d];
  return fd_->medida(d, fid) * (h0 + h1) / 2.0;
}

bool NavierStokes::faceta_interior(int d, long fid) const
{
  const auto &cel = celulas_da_faceta_[d][static_cast<size_t>(fid)];
  return cel[0].first >= 0 && cel[1].first >= 0;
}

void NavierStokes::aplica_contorno_velocidade()
{
  // UMA FACETA E' DE PAREDE quando uma das duas celulas que ela separa nao
  // existe.  Ali a velocidade NAO se resolve: e' dada (Dirichlet) ou
  // extrapolada (Neumann).
  for (int d = 0; d < 3; d++)
    for (long f = 0; f < fd_->num_facetas(d); f++) {
      if (faceta_interior(d, f)) continue;
      const Ponto3 x = fd_->centro(d, f);
      const auto &cel = celulas_da_faceta_[d][static_cast<size_t>(f)];
      const auto &c = (cel[0].first >= 0) ? cel[0] : cel[1];
      if (c.first < 0) continue;
      const Ponto3 hc = m_->lados_da_folha(c.first, c.second);
      const Contorno &cc = bcu_[d]->em(x);
      if (cc.tipo == TipoContorno::Dirichlet) {
        u_.em(d, f) = cc.g(x);
        continue;
      }
      // NEUMANN NA COMPONENTE NORMAL: `du/dn = g`, e o unico valor interior
      // disponivel esta' uma faceta adentro.  O passo e' o lado da celula
      // DONA, nao um `h` global.
      const int sgn = (cel[0].first >= 0) ? -1 : +1;   // para dentro
      Ponto3 xi = x;
      xi[d] += sgn * hc[d];
      const long fi = faceta_em(d, xi);
      if (fi < 0) continue;                             // dominio de uma celula
      u_.em(d, f) = u_.em(d, fi) - sgn * hc[d] * cc.g(x);
    }
}

void NavierStokes::preditor()
{
  // O PREDITOR NAO TOCA AS FACETAS DE PAREDE: ali a velocidade vem da
  // condicao, e sobrescreve-la com a equacao e' como a condicao deixa de
  // valer sem que nada acuse.
  ustar_ = u_;
  rebaixos_conv_ = 0;

  for (int d = 0; d < 3; d++)
    for (long f = 0; f < fd_->num_facetas(d); f++) {
      if (!faceta_interior(d, f)) continue;
      const Ponto3 x = fd_->centro(d, f);

      // A DIVISAO POR `RHO` E' SO' DOS TERMOS DE FORCA, e nao do convectivo:
      // `div(uu)` ja' e' aceleracao.  Sem campo de densidade, `rho = 1`.
      const double rho_f = 1.0;

      double rhs = -conveccao(d, f) + g_[d] + tensao_viscosa(d, f) / rho_f;

      // A PRESSAO ANTIGA ENTRA SO' NA INCREMENTAL, e e' nisso que as duas
      // diferem.  Sem ela o preditor ignora o gradiente que ja' se conhece, e
      // a projecao tem de reconstrui-lo inteiro a cada passo -- o que custa
      // uma ordem na pressao.
      if (proj_ == Projecao::Incremental)
        rhs -= gradiente_em(p_, d, f) / rho_f;

      if (fonte_faceta_) rhs += fonte_faceta_(x, d, t_);

      ustar_.em(d, f) = u_.em(d, f) + dt_ * rhs;
    }
}

void NavierStokes::passo()
{
  if (!(dt_ > 0.0))
    throw std::runtime_error(
        "NavierStokes3D::passo: o passo de tempo nao foi instalado -- chame "
        "`passo_de_tempo` primeiro");

  // O RELOGIO DO CONTORNO: `t^n` para o preditor, `t^{n+1}` para a aplicacao
  // final.  Com condicao dependente do tempo, aplicar as duas no mesmo
  // instante deixa a parede meio passo atrasada -- e o erro ali e' `O(dt)` e
  // NAO cai com `h`, de modo que ele acaba dominando a norma do maximo nas
  // malhas finas e a ordem desaba sem que a causa esteja na malha.
  if (relogio_) *relogio_ = t_;
  aplica_contorno_velocidade();
  preditor();

  // O CONTORNO VALE PARA `u*` TAMBEM: o Poisson usa a divergencia de `u*`, e
  // se a componente normal na parede nao for a da condicao, a compatibilidade
  // quebra -- com Neumann em toda a volta ele fica insoluvel, e com Dirichlet
  // em algum lugar ele RESOLVE e da' a resposta errada.
  for (int d = 0; d < 3; d++)
    for (long f = 0; f < fd_->num_facetas(d); f++)
      if (!faceta_interior(d, f)) ustar_.em(d, f) = u_.em(d, f);

  // A PROJECAO OPERA SOBRE `u*`: troca-se o campo, projeta-se, e o resultado
  // ja' e' `u^{n+1}`.
  u_ = ustar_;
  resolve_pressao();
  if (!conv_pressao_)
    throw std::runtime_error(
        "NavierStokes3D::passo: o Poisson da projecao nao convergiu");
  projeta();

  // A PRESSAO: na incremental acumula o incremento; na nao incremental o
  // incremento E' a pressao.
  for (long b = 0; b < m_->num_blocos(); b++)
    for (long i = 0; i < m_->num_folhas_do_bloco(b); i++) {
      const long id = m_->id_global(b, i);
      p_.em(id) = (proj_ == Projecao::Incremental) ? p_.em(id) + dp_.em(id)
                                                   : dp_.em(id);
    }

  t_ += dt_;
  n_passos_++;
  if (relogio_) *relogio_ = t_;
  aplica_contorno_velocidade();
}

void NavierStokes::resolve_pressao()
{
  // O POISSON DA PROJECAO: `lap(dp) = div(u)/dt`, MULTIPLICADO PELO VOLUME.
  //
  // Montar `area/dist` dos dois lados -- e nao dividir pelo volume -- deixa a
  // matriz SIMETRICA, que e' o que o par `div`/`grad` sustenta em malha
  // uniforme (clausula O1 no lado 2D).  Dividir pelo volume daria a mesma
  // solucao e uma matriz nao simetrica, de graca.
  const long n = m_->num_folhas();
  const long p0 = m_->primeiro_id();
  t8flow::Solver s(t8flow::Motor::Petsc, static_cast<int>(p0),
                   static_cast<std::size_t>(n), 7);
  s.tolerancia(1e-12);

  for (long b = 0; b < m_->num_blocos(); b++)
    for (long i = 0; i < m_->num_folhas_do_bloco(b); i++) {
      const long id = m_->id_global(b, i);
      const Ponto3 c = m_->centroide(b, i);
      const Ponto3 hh = m_->lados_da_folha(b, i);

      std::map<long, double> linha;
      double rhs = 0.0;

      for (int d = 0; d < 3; d++)
        for (int lado = 0; lado < 2; lado++) {
          const double sinal = (lado == 0) ? -1.0 : 1.0;
          const long f = faceta_da_face(b, i, d, lado);
          if (f >= 0) rhs += sinal * u_.em(d, f) * fd_->medida(d, f) / dt_;

          const std::vector<long> vz = m_->ids_vizinhos(b, i, 2 * d + lado);
          if (vz.empty()) {
            // ---- parede ------------------------------------------------
            Ponto3 xp = c;
            xp[d] += sinal * hh[d] / 2.0;
            const Contorno &cp = bcp_->em(xp);
            if (cp.tipo != TipoContorno::Dirichlet) continue;  // fluxo nulo

            // A CONDICAO DO INCREMENTO NAO E' A DA PRESSAO, E' DERIVADA DELA.
            //
            // Queremos `p^{n+1} = g` na parede.  Na NAO INCREMENTAL `p = dp`,
            // entao o valor do incremento ali E' `g`.  Na INCREMENTAL
            // `p^{n+1} = p^n + dp`, logo `dp_parede = g - p^n_parede`.
            //
            // DIRICHLET HOMOGENEO -- ou `g` cru -- PRENDE A VARIACAO, E NAO A
            // PRESSAO: no ponto fixo `dp = 0` resolve a equacao qualquer que
            // seja o nivel de `p`, e o campo sai deslocado por uma constante.
            // O sintoma e' gradiente certo com campo errado, e diagnostica-se
            // pela FAIXA de `p`, nao pelo erro.
            //
            // ESTE DEFEITO JA' CUSTOU UM DIAGNOSTICO NO LADO 2D, e eu o
            // reintroduzi aqui ao portar.  Esta' registrado em
            // `doc/relatorio-ns-fv-2d.md`, secao 4.
            //
            // A EXTRAPOLACAO DE `p^n` ATE' A PAREDE E' LINEAR, `(3 p_c -
            // p_cc)/2`, e nao de ordem zero: usar `p_c` seria primeira ordem
            // disfarcada de segunda, e apareceria como a pressao convergindo
            // em ordem 1 enquanto a velocidade convergisse em 2.
            double val = cp.g(xp);
            if (proj_ == Projecao::Incremental) {
              Ponto3 xcc = c;
              xcc[d] -= sinal * hh[d];        // o segundo centro, adentro
              const long cc = m_->id_do_ponto(xcc);
              const double p_parede =
                  (cc >= 0) ? 1.5 * p_.em(id) - 0.5 * p_.em(cc) : p_.em(id);
              val -= p_parede;
            }

            const double area = (f >= 0) ? fd_->medida(d, f)
                                         : hh[(d + 1) % 3] * hh[(d + 2) % 3];
            const double coef = area / (hh[d] / 2.0);
            linha[id] -= coef;
            rhs -= coef * val;
            continue;
          }
          const auto it = bi_por_id_.find(vz[0]);
          if (it == bi_por_id_.end()) continue;
          const Ponto3 hv =
              m_->lados_da_folha(it->second.first, it->second.second);
          const double area = (f >= 0) ? fd_->medida(d, f)
                                       : hh[(d + 1) % 3] * hh[(d + 2) % 3];
          const double coef = area / ((hh[d] + hv[d]) / 2.0);
          linha[id] -= coef;
          linha[vz[0]] += coef;
        }

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
  // `u <- u - dt grad(dp)`, e e' AQUI que a divergencia vai a zero.
  //
  // O gradiente e' o MESMO operador que montou a matriz.  Se a projecao usasse
  // um gradiente diferente do que a divergencia do Poisson supos, o campo
  // corrigido nao teria divergencia nula e NADA DIRIA -- o solver convergiria
  // igual.  E' o par que fecha, nao cada um.
  for (int d = 0; d < 3; d++)
    for (long f = 0; f < fd_->num_facetas(d); f++)
      u_.em(d, f) -= dt_ * gradiente_em(dp_, d, f);
}

void NavierStokes::projeta_campo()
{
  if (!(dt_ > 0.0))
    throw std::invalid_argument(
        "NavierStokes3D: projetar exige passo de tempo -- chame "
        "`passo_de_tempo` antes");
  resolve_pressao();
  if (!conv_pressao_)
    throw std::runtime_error(
        "NavierStokes3D: o Poisson da projecao nao convergiu");
  projeta();
  for (long b = 0; b < m_->num_blocos(); b++)
    for (long i = 0; i < m_->num_folhas_do_bloco(b); i++) {
      const long id = m_->id_global(b, i);
      p_.em(id) += dp_.em(id);
    }
}

long NavierStokes::faceta_da_face(long b, long i, int d, int lado) const
{
  // A FACETA DE UMA FACE, VENHA ELA DE QUEM VIER.
  //
  // `DominioFacetas::id_faceta` devolve -1 quando ESTA folha nao e' dona da
  // faceta -- o que e' informacao e nao falha: e' assim que o chamador sabe
  // que ela pertence ao vizinho.  Mas o balanco de fluxos precisa das SEIS
  // faces da celula, donas ou nao: descartar as que nao sao dela deixa metade
  // do balanco de fora.
  //
  // MEDIDO AO ESQUECER DISTO: campo de divergencia 3 devolvia 18, e campo
  // solenoidal devolvia 10 em vez de zero.  Nao e' erro pequeno -- e' o
  // operador somando um subconjunto arbitrario de faces --, mas num caso de
  // ordem seria facil confundir com malha grossa.
  const long f = fd_->id_faceta(b, i, d, lado);
  if (f >= 0) return f;

  // Nao e' dona: a faceta e' a do OUTRO lado da celula vizinha por esta face.
  // A convencao de face e' a do t8code -- 0/1 em x, 2/3 em y, 4/5 em z.
  const std::vector<long> vz = m_->ids_vizinhos(b, i, 2 * d + lado);
  if (vz.size() != 1) return -1;      // contorno, ou transicao (recusada)
  const auto it = bi_por_id_.find(vz[0]);
  if (it == bi_por_id_.end()) return -1;
  return fd_->id_faceta(it->second.first, it->second.second, d, 1 - lado);
}

double NavierStokes::divergencia_em(long b, long i) const
{
  // O FLUXO LIQUIDO PELAS SEIS FACES, DIVIDIDO PELO VOLUME.
  //
  // A MEDIDA DA FACETA E' UMA AREA em tres dimensoes, e o volume e' o produto
  // dos tres lados.  Em malha uniforme `area/volume` e' `1/h` tanto quanto
  // seria com medida de comprimento e "volume" de area -- isto e', trocar as
  // duas por engano daria o MESMO numero aqui, e so' uma malha de celulas
  // nao cubicas separa as duas contas.  Por isso o caso que verifica isto usa
  // blocos de lados diferentes.
  const Ponto3 hh = m_->lados_da_folha(b, i);
  const double vol = hh[0] * hh[1] * hh[2];

  double fluxo = 0.0;
  for (int d = 0; d < 3; d++)
    for (int lado = 0; lado <= 1; lado++) {
      // O SINAL E' O DA NORMAL EXTERIOR: a face "menos" aponta para fora no
      // sentido negativo de `d`.
      const double sinal = (lado == 0) ? -1.0 : 1.0;
      const long f = faceta_da_face(b, i, d, lado);
      if (f < 0) continue;            // face de contorno sem faceta propria
      fluxo += sinal * u_.em(d, f) * fd_->medida(d, f);
    }
  return fluxo / vol;
}

double NavierStokes::divergencia_maxima() const
{
  double maxv = 0.0;
  for (long b = 0; b < m_->num_blocos(); b++)
    for (long i = 0; i < m_->num_folhas_do_bloco(b); i++)
      maxv = std::max(maxv, std::fabs(divergencia_em(b, i)));
  return maxv;
}

}  // namespace t8flow3d
