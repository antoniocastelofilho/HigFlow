#include "poisson3d.hxx"

#include <cmath>
#include <memory>
#include <stdexcept>

#include "estencil3d.hxx"

namespace t8flow3d {

Poisson::Poisson(const Malha3 &m, ContornoPorRegiao bc)
    : m_(&m), contorno_(std::move(bc))
{
  if (contorno_.vazio())
    throw std::invalid_argument(
        "Poisson: contorno vazio -- sem condicao em parede nenhuma o sistema e' "
        "singular, e o solver gastaria as iteracoes todas para dizer isso mal");
  primeiro_ = m.primeiro_id();
  n_ = m.num_folhas();
}

void Poisson::discretizacao(const Discretizacao &d)
{
  if (d.familia() != Familia::VolumesFinitos)
    throw std::invalid_argument(
        "Poisson3D: \"" + d.nome() + "\" nao e' de volumes finitos.  O lado 3D "
        "traz so' volumes finitos, por decisao de escopo: os compactos recusam "
        "malha adaptada e o laplaciano pontual e' instavel nela");
  disc_ = d;
}

// ---------------------------------------------------------------------------
// BALANCO DE FLUXOS DE SEGUNDA ORDEM
// ---------------------------------------------------------------------------

void Poisson::monta_linha_fluxo(long b, long i, std::map<long, double> *linha_out,
                                double *rhs_out) const
{
  const Ponto3 c = m_->centroide(b, i);
  const Ponto3 h = m_->lados_da_folha(b, i);
  const double vol = h[0] * h[1] * h[2];
  const int nivel_self = m_->nivel(b, i);

  std::map<long, double> linha;
  // A FONTE E' INTEGRADA sobre a celula -- volumes finitos.  Regra do ponto
  // medio, que e' exata para `f` linear e casa com a ordem do esquema.
  double rhs = f_(c) * vol;

  for (int fa = 0; fa < 6; fa++) {
    const int d = fa / 2;
    const double sinal = (fa % 2 == 0) ? -1.0 : 1.0;
    // A AREA DA FACETA e' o produto dos DOIS lados transversais.  Em duas
    // dimensoes era um lado so', e traduzir aquilo daria um fluxo com a
    // dimensao errada -- que ainda converge, para a equacao errada.
    const double area_face = h[(d + 1) % 3] * h[(d + 2) % 3];

    const std::vector<long> vz = m_->ids_vizinhos(b, i, fa);
    const std::vector<int> nvz = m_->niveis_vizinhos(b, i, fa);

    if (vz.empty()) {
      // CONTORNO.  A relacao afim `u_fantasma = A u + B` cobre os tres tipos
      // sem ramo proprio, e o fluxo pela parede vira
      // `area/h * (u_fantasma - u) = area/h * ((A-1) u + B)`.
      Ponto3 parede = c;
      parede[d] += sinal * h[d] / 2.0;
      const Contorno &cond = contorno_.em(parede);

      double A = 0.0, B = 0.0;
      cond.afim(parede, h[d], &A, &B);
      const double coef = area_face / h[d];
      linha[m_->id_global(b, i)] += coef * (A - 1.0);
      rhs -= coef * B;
      continue;
    }

    for (size_t k = 0; k < vz.size(); k++) {
      const int nv = (k < nvz.size()) ? nvz[k] : nivel_self;

      // NA TRANSICAO 2:1 UMA FACE TEM ATE' QUATRO VIZINHOS, e a razao de AREAS
      // e' `1/4` e nao `1/2` como em duas dimensoes.  Os tres fatores mudam
      // juntos, em termos de `fator = 2^(nivel_self - nivel_viz)`:
      //
      //     h_viz = fator * h        (o lado da celula vizinha)
      //     |face| = area * min(1, fator^2)   <- AREA, logo o fator ao QUADRADO
      //     dist  = (h + h_viz) / 2
      const double fator = std::pow(2.0, static_cast<double>(nivel_self - nv));
      const double h_viz = fator * h[d];
      const double face = area_face * std::min(1.0, fator * fator);
      const double dist = (h[d] + h_viz) / 2.0;
      const double coef = face / dist;

      linha[vz[k]] += coef;
      linha[m_->id_global(b, i)] -= coef;
    }
  }

  *linha_out = linha;
  *rhs_out = rhs;
}

// ---------------------------------------------------------------------------
// RECONSTRUCAO MLS COM QUADRATURA DE GAUSS
// ---------------------------------------------------------------------------

void Poisson::monta_linha_mls(long b, long i, const void *construtor,
                              std::map<long, double> *linha_out,
                              double *rhs_out) const
{
  const ConstrutorEstencil &ce = *static_cast<const ConstrutorEstencil *>(construtor);

  const Ponto3 c = m_->centroide(b, i);
  const Ponto3 h = m_->lados_da_folha(b, i);
  const double vol = h[0] * h[1] * h[2];

  std::map<long, double> linha;

  // GAUSS DE DOIS PONTOS por direcao: nos em `+-1/sqrt(3)`, pesos iguais.
  const double no = 1.0 / std::sqrt(3.0);
  const double wg = 0.5;

  // A FONTE POR GAUSS 2x2x2, e nao pela regra do ponto medio.
  //
  // EM DUAS DIMENSOES ESTA FOI A DESCOBERTA QUE DESTRAVOU O ESQUEMA: com o fluxo
  // da faceta ja' exato e a fonte no ponto medio, ele TRAVAVA em segunda ordem
  // (9,58e-05, ordem 1,99, contra 2,52e-14 exato).  A regra do ponto medio e' um
  // erro `O(h^2)` RELATIVO em `integral f dV`, e no esquema de segunda ordem ele
  // se escondia atras do erro do proprio fluxo.
  double rhs = 0.0;
  for (int qx = 0; qx < 2; qx++)
    for (int qy = 0; qy < 2; qy++)
      for (int qz = 0; qz < 2; qz++) {
        const Ponto3 xg{c[0] + (qx == 0 ? -no : no) * h[0] / 2.0,
                        c[1] + (qy == 0 ? -no : no) * h[1] / 2.0,
                        c[2] + (qz == 0 ? -no : no) * h[2] / 2.0};
        rhs += wg * wg * wg * f_(xg);
      }
  rhs *= vol;

  for (int fa = 0; fa < 6; fa++) {
    const int d = fa / 2;
    const int t1 = (d + 1) % 3;
    const int t2 = (d + 2) % 3;
    const double sinal = (fa % 2 == 0) ? -1.0 : 1.0;

    // AS SUB-FACETAS.  Em transicao 2:1 a face da celula grossa e' coberta por
    // QUATRO sub-facetas em tres dimensoes -- duas por direcao transversal --,
    // e o fluxo total e' a soma.  Integrar a face inteira de uma vez daria o
    // fluxo medio sobre ela, que nao e' a soma dos quatro fluxos finos quando a
    // derivada varia ao longo dela.
    const std::vector<int> nvz = m_->niveis_vizinhos(b, i, fa);
    const int nivel_self = m_->nivel(b, i);
    int nsub = 1;
    if (!nvz.empty() && nvz[0] > nivel_self) nsub = 1 << (nvz[0] - nivel_self);

    const double larg1 = h[t1] / nsub;
    const double larg2 = h[t2] / nsub;

    for (int s1 = 0; s1 < nsub; s1++)
      for (int s2 = 0; s2 < nsub; s2++) {
        Ponto3 cf = c;
        cf[d] += sinal * h[d] / 2.0;
        cf[t1] += -h[t1] / 2.0 + (s1 + 0.5) * larg1;
        cf[t2] += -h[t2] / 2.0 + (s2 + 0.5) * larg2;

        // QUADRATURA 2x2 NA FACETA -- em duas dimensoes eram dois pontos num
        // segmento; aqui sao quatro num retangulo.
        for (int q1 = 0; q1 < 2; q1++)
          for (int q2 = 0; q2 < 2; q2++) {
            Ponto3 xg = cf;
            xg[t1] += (q1 == 0 ? -no : no) * larg1 / 2.0;
            xg[t2] += (q2 == 0 ? -no : no) * larg2 / 2.0;

            // O ESTENCIL DA DERIVADA NO PONTO DE QUADRATURA.  O suporte depende
            // so' do PONTO, logo a celula vizinha obtem o mesmo estencil para o
            // mesmo ponto -- e' daqui que vem a conservacao.
            const Estencil e = ce.derivada_em(xg, d);

            const double peso = sinal * wg * wg * larg1 * larg2;
            for (size_t k = 0; k < e.ids.size(); k++)
              linha[e.ids[k]] += peso * e.pesos[k];
            rhs -= peso * e.rhs;
          }
      }
  }

  *linha_out = linha;
  *rhs_out = rhs;
}

void Poisson::monta_linha(long b, long i, const void *construtor,
                          std::map<long, double> *linha, double *rhs) const
{
  if (disc_.variante() == Variante::MlsGauss)
    monta_linha_mls(b, i, construtor, linha, rhs);
  else
    monta_linha_fluxo(b, i, linha, rhs);
}

// ---------------------------------------------------------------------------
// RESOLUCAO
// ---------------------------------------------------------------------------

void Poisson::resolve(Motor motor, double tol, int max_it)
{
  if (!f_) throw std::invalid_argument("Poisson::resolve: fonte nao instalada");

  const bool mlsg = (disc_.variante() == Variante::MlsGauss);

  std::unique_ptr<ConstrutorEstencil> ce;
  if (mlsg) {
    ce.reset(new ConstrutorEstencil(*m_, disc_.ordem_mls(), mls_maxpts_, 6));
    ce->contorno(contorno_);
    ce->mls().zera_contadores();
    ce->zera_contadores_de_franja();
    // AS CAMADAS DE FRANJA QUE O SUPORTE PEDE.  O numero sai da algebra, nao de
    // escolha: `camadas_necessarias()` e' `raio + 1`, e o raio e' `ceil(ordem/2)`
    // porque uma caixa de raio `r` tem `2r+1` valores distintos por direcao.
    const_cast<Malha3 *>(m_)->exige_camadas(ce->camadas_necessarias());
  } else {
    // O FLUXO PRECISA DE UMA CAMADA so': a linha dele toca apenas os vizinhos de
    // face.  Em tres dimensoes nao ha' correcao transversal a montar -- ela
    // existia no lado 2D e aqui o termo e' absorvido pela quadratura de
    // sub-faceta, que o `mlsgauss` faz e o `fluxo2` dispensa por ser de segunda
    // ordem.
    const_cast<Malha3 *>(m_)->exige_camadas(1);
  }

  // O MLS DA' LINHAS MUITO MAIS CHEIAS: ate' 125 pontos de suporte por ponto de
  // quadratura, 24 pontos de quadratura por celula.  A pre-alocacao tem de
  // acompanhar, senao o PETSc realoca a cada insercao.
  const int nnz = mlsg ? 400 : 10;
  Solver s(motor, static_cast<int>(primeiro_), static_cast<std::size_t>(n_), nnz);
  s.tolerancia(tol);
  s.max_iteracoes(max_it);

  for (long b = 0; b < m_->num_blocos(); b++)
    for (long i = 0; i < m_->num_folhas_do_bloco(b); i++) {
      std::map<long, double> linha;
      double rhs = 0.0;
      monta_linha(b, i, ce.get(), &linha, &rhs);

      std::vector<int> cols;
      std::vector<double> vals;
      cols.reserve(linha.size());
      vals.reserve(linha.size());
      for (const auto &kv : linha) {
        cols.push_back(static_cast<int>(kv.first));
        vals.push_back(kv.second);
      }
      const long id = m_->id_global(b, i);
      s.linha(static_cast<int>(id), cols, vals);
      s.lado_direito(static_cast<int>(id), rhs);
    }

  if (mlsg) {
    rebaixos_ = ce->mls().rebaixamentos();
    truncados_ = static_cast<long>(reduz_soma(
        static_cast<double>(ce->truncados_pela_franja())));
  }

  s.monta();
  s.resolve();
  convergiu_ = s.convergiu();
  iteracoes_ = s.iteracoes();
  residuo_ = s.residuo();
  sol_ = s.solucao();
}

double Poisson::u(long id) const
{
  const long k = id - primeiro_;
  if (k < 0 || k >= static_cast<long>(sol_.size()))
    throw std::out_of_range("Poisson::u: id " + std::to_string(id) +
                            " fora da faixa local");
  return sol_[static_cast<size_t>(k)];
}

double Poisson::reduz_max(double v) const
{
  if (m_->num_processos() <= 1) return v;
  double r = v;
  sc_MPI_Allreduce(&v, &r, 1, sc_MPI_DOUBLE, sc_MPI_MAX, m_->comunicador());
  return r;
}

double Poisson::reduz_soma(double v) const
{
  if (m_->num_processos() <= 1) return v;
  double r = v;
  sc_MPI_Allreduce(&v, &r, 1, sc_MPI_DOUBLE, sc_MPI_SUM, m_->comunicador());
  return r;
}

double Poisson::erro_inf(const std::function<double(const Ponto3 &)> &exata) const
{
  double pior = 0.0;
  for (long b = 0; b < m_->num_blocos(); b++)
    for (long i = 0; i < m_->num_folhas_do_bloco(b); i++)
      pior = std::fmax(pior, std::fabs(u(m_->id_global(b, i)) -
                                       exata(m_->centroide(b, i))));
  // A REDUCAO E' PARTE DA NORMA: sem ela cada processo devolve o pior erro da
  // SUA particao, e o numero sai MENOR que o verdadeiro -- num relatorio de
  // verificacao, pior que um erro grande.
  return reduz_max(pior);
}

double Poisson::erro_l2(const std::function<double(const Ponto3 &)> &exata) const
{
  double soma = 0.0, vol_total = 0.0;
  for (long b = 0; b < m_->num_blocos(); b++)
    for (long i = 0; i < m_->num_folhas_do_bloco(b); i++) {
      const double e = u(m_->id_global(b, i)) - exata(m_->centroide(b, i));
      const double v = m_->volume(b, i);
      soma += e * e * v;
      vol_total += v;
    }
  // AS DUAS SOMAS REDUZEM, e o quociente vem DEPOIS: reduzir o quociente daria
  // a media das normas por processo, que nao e' a norma.
  soma = reduz_soma(soma);
  vol_total = reduz_soma(vol_total);
  return (vol_total > 0.0) ? std::sqrt(soma / vol_total) : 0.0;
}

double Poisson::residuo_da_exata(
    const std::function<double(const Ponto3 &)> &exata) const
{
  if (!f_) throw std::invalid_argument("Poisson::residuo_da_exata: fonte nao instalada");

  const bool mlsg = (disc_.variante() == Variante::MlsGauss);
  std::unique_ptr<ConstrutorEstencil> ce;
  if (mlsg) {
    ce.reset(new ConstrutorEstencil(*m_, disc_.ordem_mls(), mls_maxpts_, 6));
    ce->contorno(contorno_);
    const_cast<Malha3 *>(m_)->exige_camadas(ce->camadas_necessarias());
  }

  double pior = 0.0;
  for (long b = 0; b < m_->num_blocos(); b++)
    for (long i = 0; i < m_->num_folhas_do_bloco(b); i++) {
      std::map<long, double> linha;
      double rhs = 0.0;
      monta_linha(b, i, ce.get(), &linha, &rhs);

      // A SOMA PERCORRE A LINHA INTEIRA, inclusive as colunas de FRANJA.
      // Descartar as de fora da faixa local faria o residuo sair MENOR que o
      // verdadeiro -- um operador inconsistente se apresentando como
      // consistente, que e' o oraculo mentindo a favor do codigo.
      double soma = 0.0;
      for (const auto &kv : linha) {
        const long vx = m_->indice_do_id(kv.first);
        if (vx < 0)
          throw std::runtime_error(
              "Poisson::residuo_da_exata: a linha referencia celula que nao esta' "
              "nas folhas locais nem na franja -- o residuo seria calculado com "
              "termos FALTANDO");
        soma += kv.second * exata(m_->centroide_do_indice(vx));
      }
      // NORMALIZA PELO VOLUME: em volumes finitos a linha foi multiplicada por
      // `|face|`, logo o residuo bruto cresce com o tamanho da celula.  Quem
      // sabe disso e' a `Discretizacao`.
      const double norm = disc_.normalizacao_do_residuo(m_->volume(b, i));
      pior = std::fmax(pior, std::fabs(soma - rhs) / norm);
    }
  return reduz_max(pior);
}

}  // namespace t8flow3d
