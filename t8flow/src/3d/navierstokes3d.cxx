#include "navierstokes3d.hxx"

#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "solverlin.hxx"

namespace t8flow3d {

NavierStokes::NavierStokes(const Malha3 &m, const DominioFacetas &fd,
                           const ContornoPorRegiao &bc_u0,
                           const ContornoPorRegiao &bc_u1,
                           const ContornoPorRegiao &bc_u2,
                           const ContornoPorRegiao &bc_p)
    : m_(&m), fd_(&fd), u_(fd), p_(m), dp_(m)
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
        "finas, nao duas, e a quina onde a viscosidade se interpola toca OITO "
        "celulas, nao quatro -- as regras de 2D aplicadas aqui produzem "
        "divergencia que nao fecha, e o sintoma aparece como o solver de "
        "pressao nao convergindo, longe da causa");

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

  u_.preenche(0.0);
  p_.preenche(0.0);
  dp_.preenche(0.0);
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
            // Dirichlet: o valor esta' na parede, a meio lado do centro.
            const double area = (f >= 0) ? fd_->medida(d, f)
                                         : hh[(d + 1) % 3] * hh[(d + 2) % 3];
            const double coef = area / (hh[d] / 2.0);
            linha[id] -= coef;
            rhs -= coef * cp.g(xp);
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
