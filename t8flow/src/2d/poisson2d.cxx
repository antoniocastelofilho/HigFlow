// Implementacao do Poisson.  A DERIVACAO do fechamento de parede e a razao de o
// operador ser um balanco de fluxos estao no `poisson2d.hxx`.  Aqui ficam as
// decisoes de implementacao:
//
//   1. a linha de cada celula e' acumulada num mapa `id -> coeficiente` antes de
//      ir ao solver, porque o MESMO id pode receber contribuicao de mais de uma
//      faceta -- o vizinho interior usado no fechamento de uma parede e' tambem
//      vizinho por outra face.  Somar e' o correto; sobrescrever perderia um termo
//      em silencio;
//   2. o termo conhecido (o valor de Dirichlet) e' movido para o lado direito com
//      sinal trocado, uma vez, no ponto em que e' calculado.  Acumular e trocar o
//      sinal depois e' onde o sinal se perde;
//   3. a montagem RECUSA o que esta versao nao sabe: faceta com dois vizinhos
//      (transicao 2:1), celula de parede sem vizinho interior, e contorno que nao
//      seja Dirichlet.  Recusar e' o que impede um resultado plausivel de sair de
//      um caso nao suportado.
#include "poisson2d.hxx"

#include <cmath>
#include <map>
#include <utility>
#include <vector>
#include <stdexcept>
#include <string>

namespace t8flow2d {

Poisson::Poisson(const Malha &m, ContornoPorRegiao contorno)
    : m_(&m), contorno_(std::move(contorno)),
      primeiro_(m.primeiro_id()), n_(m.num_folhas())
{
  if (contorno_.vazio())
    throw std::invalid_argument("Poisson: nenhuma condicao de contorno instalada");
  sol_.assign(static_cast<size_t>(n_), 0.0);
}

namespace {

//! (bloco, folha) de um id global.  O(folhas) -- caminho de montagem, nao laco
//! quente; dito para nao ser suposto rapido.
bool localiza(const Malha &m, long id, long *b, long *i)
{
  for (long bb = 0; bb < m.num_blocos(); bb++)
    for (long ii = 0; ii < m.num_folhas_do_bloco(bb); ii++)
      if (m.id_global(bb, ii) == id) { *b = bb; *i = ii; return true; }
  return false;
}

}  // namespace

//! Indice `id global -> (bloco, folha)`, montado UMA vez por resolucao.
//
// Sem ele, cada consulta de vizinho custaria O(folhas) e a montagem viraria
// O(N^2): com 2239 celulas isso e' 5 milhoes de comparacoes por faceta.  O indice
// e' O(N) de memoria e O(1) por consulta.
std::vector<std::pair<long, long>> monta_indice(const Malha &m)
{
  std::vector<std::pair<long, long>> idx(static_cast<size_t>(m.num_folhas()), {-1, -1});
  const long p0 = m.primeiro_id();
  for (long b = 0; b < m.num_blocos(); b++)
    for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
      const long k = m.id_global(b, i) - p0;
      if (k >= 0 && k < static_cast<long>(idx.size())) idx[static_cast<size_t>(k)] = {b, i};
    }
  return idx;
}

//! Monta a linha de UMA celula: coeficientes por id, e o lado direito.
//
// Extraida de `resolve` para que o MESMO codigo sirva ao oraculo de consistencia
// (`residuo_da_exata`).  Duas montagens divergiriam, e divergir aqui significaria
// o teste medir um operador e o solver resolver outro -- que e' o pior resultado
// possivel: verde num operador que nao e' o que roda.
void Poisson::monta_linha(long b, long i, const void *indice,
                          std::map<long, double> *linha_out, double *rhs_out) const
{
  const std::vector<std::pair<long, long>> &idx =
      *static_cast<const std::vector<std::pair<long, long>> *>(indice);
  const long p0 = m_->primeiro_id();
  const long id = m_->id_global(b, i);
  const Ponto2 c = m_->centroide(b, i);
  const Ponto2 h = m_->lados_da_folha(b, i);
  const double area = h[0] * h[1];

  // `int_V lap(u) = sum_faces du/dn |face|`, logo a linha e' a soma dos fluxos
  // e o lado direito e' `f(c) * area`.
  std::map<long, double> linha;
  double rhs = f_(c) * area;

  for (int d = 0; d < 2; d++) {
    const double hd = h[d];
    const double larg = h[1 - d];     // comprimento da faceta normal a d
    for (int lado = 0; lado < 2; lado++) {
      const int fa = 2 * d + lado;
      const std::vector<long> vz = m_->ids_vizinhos(b, i, fa);

      if (!vz.empty()) {
        // FACETA INTERNA, UNIFORME OU EM TRANSICAO 2:1, pela MESMA formula.
        //
        // Para cada vizinho `v`, a face compartilhada tem o tamanho do MENOR
        // dos dois, e a distancia entre os centros e' a soma dos dois
        // meios-passos:
        //
        //     |face| = min(larg_self, larg_viz)
        //     dist   = (h_self + h_viz) / 2
        //     coef   = |face| / dist
        //
        // UMA formula e nao tres (uniforme, grossa->fina, fina->grossa): tres
        // divergem, e divergir aqui quebra a CONSERVACAO sem quebrar nada
        // visivel.  Conferido que ela e' conservativa -- o coeficiente que a
        // grossa usa para cada fina e' IGUAL ao que a fina usa para a grossa:
        //
        //     grossa->fina:  larg_C/2  /  (3 h_C/4)  =  (2/3) larg_C / h_C
        //     fina->grossa:  larg_F    /  (3 h_F/2)  =  (2/3) larg_F / h_F
        //                                            =  (2/3) larg_C / h_C
        //
        // Logo o que sai de uma entra na outra com o mesmo valor, e a matriz
        // fica simetrica.  Sem essa igualdade o esquema NAO conserva, e o
        // sintoma seria um erro que nao cai com o refino -- facil de confundir
        // com defeito de contorno.
        //
        // O tamanho do vizinho vem do NIVEL dele, nao de uma busca: com 2:1 a
        // razao e' 2^(nivel_self - nivel_viz), e os niveis `niveis_vizinhos`
        // ja' da'.  Buscar (bloco, folha) de cada vizinho seria O(folhas) por
        // faceta, dentro do laco de montagem.
        const std::vector<int> nvz = m_->niveis_vizinhos(b, i, fa);
        const int nivel_self = m_->nivel(b, i);
        const int t = 1 - d;              // direcao TRANSVERSAL a' faceta

        for (size_t k = 0; k < vz.size(); k++) {
          const int dn = nivel_self - (k < nvz.size() ? nvz[k] : nivel_self);
          const double fator = std::pow(2.0, static_cast<double>(dn));
          const double h_viz = hd * fator;
          const double face = larg * std::fmin(1.0, fator);
          const double dist = (hd + h_viz) / 2.0;
          const double coef = face / dist;

          linha[vz[k]] += coef;
          linha[id]    -= coef;

          // ---- A CORRECAO TRANSVERSAL, e ela NAO e' opcional ---------------
          //
          // Na transicao 2:1 os dois centros NAO estao alinhados com a normal:
          // o centro da celula fina esta' deslocado um quarto de celula do
          // centro da grossa.  A diferenca `(u_viz - u_self)/dist` mede a
          // derivada ao longo da LINHA DOS CENTROS, que tem componente
          // transversal -- e para a celula FINA esse erro NAO cancela, porque
          // ela ve' uma unica face contra a grossa.
          //
          // MEDIDO sem esta correcao: o grau 1, que o esquema reproduz
          // exatamente em malha uniforme, deu erro 3,9e-02 na malha adaptada, e
          // a solucao trigonometrica teve ordem NEGATIVA -- nao convergia.  Nao
          // era erro de O(h): era de O(1).
          //
          // A correcao desloca o valor do lado GROSSO para a linha transversal
          // do centro da sub-face, usando o gradiente transversal:
          //
          //     u_grosso*  =  u_grosso + delta * (du/dt)_grosso
          //
          // Com `du/dt` linear, isso e' EXATO para campo linear -- que e' o
          // minimo que se exige de um esquema de segunda ordem.
          const long kv = vz[k] - p0;
          if (kv < 0 || kv >= static_cast<long>(idx.size())) continue;
          const long vb = idx[static_cast<size_t>(kv)].first;
          const long vi = idx[static_cast<size_t>(kv)].second;
          if (vb < 0) continue;

          const Ponto2 cv = m_->centroide(vb, vi);
          const double delta = cv[t] - c[t];       // deslocamento transversal
          if (std::fabs(delta) < 1e-14) continue;  // alinhados: nada a corrigir

          // Quem esta' deslocado e' o lado GROSSO em relacao a' sub-face: o
          // gradiente transversal e' estimado na celula GROSSA (a de h maior),
          // e o deslocamento aplicado tem o sinal que leva o valor dela para a
          // linha do centro da fina.
          const bool self_grosso = (h_viz < hd);
          const long gb = self_grosso ? b : vb;
          const long gi = self_grosso ? i : vi;
          // O SINAL NAO DEPENDE DE QUEM E' O GROSSO, e supor que dependia foi o
          // erro da primeira versao.  Refazendo a algebra: o fluxo e'
          // `coef*(u_viz* - u_self*)`, com o valor do lado GROSSO levado para a
          // linha normal do centro da sub-face.  Nos DOIS casos a conta chega a
          //
          //     fluxo = coef*u_viz - coef*u_self - coef*delta*(du/dt)
          //
          // com `delta = c_viz[t] - c_self[t]`.  O que muda entre eles e' apenas
          // ONDE `du/dt` e' avaliado -- sempre na celula GROSSA, que e' a que
          // precisa ser deslocada.
          //
          // E o efeito do sinal trocado foi facil de reconhecer: o erro DOBROU em
          // vez de cair (7,8e-02 contra 3,9e-02 sem correcao nenhuma).  Erro que
          // dobra e' a assinatura de uma correcao invertida, e foi o que levou a
          // refazer a algebra em vez de ajustar numeros.
          const double sinal = delta;

          // GRADIENTE TRANSVERSAL na celula grossa, por diferenca entre os dois
          // lados transversais.  Cada lado contribui com a media dos vizinhos
          // dele (ou com o valor de Dirichlet, se for parede), na distancia
          // correspondente -- a mesma construcao do fechamento de parede.
          std::map<long, double> grad;
          double grad_rhs = 0.0;
          const Ponto2 cg = m_->centroide(gb, gi);
          const Ponto2 hg = m_->lados_da_folha(gb, gi);
          double dist_total = 0.0;
          bool ok = true;

          for (int lt = 0; lt < 2 && ok; lt++) {
            const int fat = 2 * t + lt;
            const double sgn = (lt == 0) ? -1.0 : 1.0;   // menos / mais
            const std::vector<long> vt = m_->ids_vizinhos(gb, gi, fat);
            if (!vt.empty()) {
              const std::vector<int> nvt = m_->niveis_vizinhos(gb, gi, fat);
              const int dnt = m_->nivel(gb, gi) - (nvt.empty() ? 0 : nvt[0]);
              const double ht = hg[t] * std::pow(2.0, static_cast<double>(dnt));
              const double dt = (hg[t] + ht) / 2.0;
              dist_total += dt;
              for (long w : vt) grad[w] += sgn / static_cast<double>(vt.size());
            } else {
              // PAREDE NO LADO TRANSVERSAL.  Tres casos, e tratar os tres e'
              // obrigatorio: desistir da correcao aqui traz de volta o erro de
              // O(1).  MEDIDO -- com tres paredes de Neumann, desistir fez o
              // residuo do operador para campo LINEAR voltar a 53,3, exatamente o
              // valor de quando nao havia correcao nenhuma.  O numero repetido foi
              // o que apontou a causa: nao era outro defeito, era o MESMO defeito
              // reaparecendo por um caminho que eu havia deixado aberto.
              Ponto2 pw = cg;
              pw[t] += sgn * hg[t] / 2.0;
              const Contorno &cw = contorno_.em(pw);

              if (cw.tipo == TipoContorno::Neumann) {
                // O MELHOR CASO, e ele e' EXATO: a normal daquela parede E' a
                // direcao transversal desta faceta, logo o dado de Neumann ali E'
                // `du/dt`, com o sinal da normal exterior.  Nao ha' diferenca a
                // montar -- o gradiente transversal e' CONHECIDO.
                grad.clear();
                // `g` e' `du/dn` com n EXTERIOR, logo `du/dt = sgn*g`: no lado
                // `mais` a normal e' +t e os dois coincidem; no lado `menos` a
                // normal e' -t e o sinal inverte.  A primeira versao escreveu
                // `g*sgn*sgn`, que cancela o sinal e deixa o erro em O(1).
                grad_rhs = sgn * cw.g(pw);
                dist_total = 1.0;                  // `g_t = grad_rhs / dist_total`
                break;
              }
              if (cw.tipo == TipoContorno::Dirichlet) {
                dist_total += hg[t] / 2.0;
                grad_rhs += sgn * cw.g(pw);
                continue;
              }
              // ROBIN nao da' `du/dt` direto.  Usa DIFERENCA DE UM LADO: este lado
              // passa a ser a propria celula, a distancia zero -- exato para campo
              // linear, que e' o piso que o esquema tem de garantir.
              grad[m_->id_global(gb, gi)] += sgn;
            }
          }

          if (!ok || dist_total <= 0.0) continue;
          const double escala = sinal * coef / dist_total;
          for (const auto &kvp : grad) linha[kvp.first] -= escala * kvp.second;
          rhs += escala * grad_rhs;
        }
        continue;
      }

      // PAREDE.  O ponto onde a condicao e' avaliada e' o centro da faceta.
      const double sinal = (lado == 0) ? -1.0 : 1.0;
      Ponto2 parede = c;
      parede[d] += sinal * hd / 2.0;

      const Contorno &cond = contorno_.em(parede);

      if (cond.tipo == TipoContorno::Neumann) {
        // NEUMANN E' O CASO TRIVIAL NO BALANCO DE FLUXOS, e isso vale dizer:
        // o que a condicao da' e' exatamente `du/dn`, que e' o proprio fluxo.
        // Nao ha' aproximacao nenhuma -- nem de dois pontos, nem de tres -- e
        // nenhum coeficiente novo entra na linha.  E' o lado em que o arranjo
        // deslocado paga: na forma pontual, Neumann exigiria fantasma.
        //
        // Com a normal EXTERIOR, `du/dn` e' o fluxo para fora direto.
        rhs -= cond.g(parede) * larg;
        continue;
      }

      if (cond.tipo != TipoContorno::Dirichlet)
        throw std::runtime_error(
            "Poisson: a parede em (" + std::to_string(parede[0]) + ", " +
            std::to_string(parede[1]) + ") e' Robin, e esta versao monta"
            " Dirichlet e Neumann.  Robin mistura valor e derivada na MESMA"
            " faceta: a linha ganha um coeficiente que depende de alfa, beta e"
            " h, e montar isso como se fosse Neumann daria resultado plausivel"
            " e errado");
      const double g = cond.g(parede);

      if (fechamento_ == FechamentoParede::DoisPontos) {
        // `du/dn` para fora = (g - u_c)/(h/2).  Classico, segunda ordem global,
        // NAO exato para grau 2 -- esta' aqui como contraste do fechamento.
        const double coef = larg / (hd / 2.0);
        linha[id] -= coef;
        rhs -= coef * g;                 // o conhecido vai para o lado direito
        continue;
      }

      // TRES PONTOS, EM POSICOES GENERICAS -- e generico nao e' luxo: na malha
      // adaptada o vizinho interior pode ser mais fino (DOIS vizinhos), mais
      // grosso, ou de mesmo nivel, e a distancia dele muda em cada caso.  A
      // primeira versao fixava `3h/2` e EXIGIA vizinho unico; na primeira malha
      // aleatoria ela lancou na celula 0, o que foi a guarda funcionando.
      //
      // Com a parede em `s=0`, o centro em `s1` e o terceiro ponto em `s2`,
      // Lagrange da'
      //
      //     du/ds(0) = g*(-s1-s2)/(s1*s2) + u_c*s2/(s1*(s2-s1))
      //                                   - u_3*s1/(s2*(s2-s1))
      //
      // e a soma dos tres coeficientes e' ZERO -- conferido, e e' a verificacao
      // que pega erro de sinal e de escala.  Com s1=h/2 e s2=3h/2 isso reduz a
      // `-8g/(3h) + 3u_c/h - u_v/(3h)`, que era a forma anterior.
      const int fa_oposta = 2 * d + (1 - lado);
      const std::vector<long> vo = m_->ids_vizinhos(b, i, fa_oposta);
      if (vo.empty())
        throw std::runtime_error(
            "Poisson: a celula " + std::to_string(id) + " nao tem vizinho"
            " interior na direcao " + std::to_string(d) + " para o fechamento de"
            " tres pontos (dominio de uma celula de espessura?).  Use"
            " FechamentoParede::DoisPontos, que precisa de dois");

      const std::vector<int> nvo = m_->niveis_vizinhos(b, i, fa_oposta);
      const int nivel_self = m_->nivel(b, i);
      const int dn = nivel_self - (nvo.empty() ? nivel_self : nvo[0]);
      const double h_viz = hd * std::pow(2.0, static_cast<double>(dn));

      // `s` E' MEDIDO DA PAREDE, nao do centro da celula -- e confundir as duas
      // origens foi exatamente o erro da primeira versao generica: ela escreveu
      // `s2 = hd/2 + h_viz/2`, que e' a distancia ENTRE CENTROS, e em malha
      // uniforme deu `hd` onde o certo e' `3hd/2`.  O caso de regressao da
      // malha uniforme (grau 2 exato) pegou na hora: 1,97e-01 em vez de 3e-13.
      //
      //     s1 = parede -> centro desta celula            = hd/2
      //     s2 = parede -> centro do vizinho interior     = s1 + (hd + h_viz)/2
      const double s1 = hd / 2.0;
      const double s2 = s1 + (hd + h_viz) / 2.0;

      const double cg = (-s1 - s2) / (s1 * s2);
      const double c1 = s2 / (s1 * (s2 - s1));
      const double c3 = -s1 / (s2 * (s2 - s1));

      // O TERCEIRO PONTO E' A MEDIA dos vizinhos interiores quando ha' mais de
      // um (lado de la' mais fino).  Eles estao simetricamente deslocados em
      // torno do eixo, entao a media aproxima o valor no centro daquela face e
      // os termos de deslocamento se cancelam -- o mesmo cancelamento que vale
      // para o fluxo nas sub-faces.  Dividir o peso igualmente e' o que
      // preserva a soma dos coeficientes em zero.
      //
      // `du/ds` e' para DENTRO, logo o fluxo para FORA e' o negativo: os sinais
      // abaixo ja' vem trocados.
      linha[id] += -c1 * larg;
      const double peso_cada = (-c3 * larg) / static_cast<double>(vo.size());
      for (long v : vo) linha[v] += peso_cada;
      rhs -= (-cg) * g * larg;
    }
  }


  *linha_out = linha;
  *rhs_out = rhs;
}

void Poisson::resolve(Motor motor, double tol, int max_it)
{
  if (!f_) throw std::invalid_argument("Poisson::resolve: termo de fonte nao instalado");
  const std::vector<std::pair<long, long>> idx = monta_indice(*m_);

  Solver s(motor, static_cast<int>(primeiro_), static_cast<std::size_t>(n_), 25);
  s.tolerancia(tol);
  s.max_iteracoes(max_it);

  for (long b = 0; b < m_->num_blocos(); b++)
    for (long i = 0; i < m_->num_folhas_do_bloco(b); i++) {
      std::map<long, double> linha;
      double rhs = 0.0;
      monta_linha(b, i, &idx, &linha, &rhs);

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

  s.monta();
  s.resolve();

  convergiu_ = s.convergiu();
  iteracoes_ = s.iteracoes();
  residuo_ = s.residuo();
  sol_ = s.solucao();
}

//! O ORACULO DIRETO DO OPERADOR: `max |A*u_exata - b|`, sem passar pelo solver.
//
// Mede a CONSISTENCIA do esquema, nao a da solucao -- e a diferenca importa: o
// erro da solucao mistura consistencia, condicionamento e tolerancia do solver, e
// quando ele e' grande nao se sabe qual dos tres.  Aqui, se o esquema reproduz a
// solucao exata, isto e' ZERO a menos de arredondamento, qualquer que seja o
// solver.
//
// E' o oraculo que isola defeito de MONTAGEM, que foi exatamente o que faltou
// quando o esquema de transicao 2:1 deu ordem negativa: medindo so' a solucao, nao
// havia como separar "o operador esta' errado" de "o sistema esta' mal
// condicionado".
double Poisson::residuo_da_exata(const std::function<double(const Ponto2 &)> &exata) const
{
  if (!f_) throw std::invalid_argument("Poisson::residuo_da_exata: fonte nao instalada");
  const std::vector<std::pair<long, long>> idx = monta_indice(*m_);
  const long p0 = m_->primeiro_id();

  double pior = 0.0;
  for (long b = 0; b < m_->num_blocos(); b++)
    for (long i = 0; i < m_->num_folhas_do_bloco(b); i++) {
      std::map<long, double> linha;
      double rhs = 0.0;
      monta_linha(b, i, &idx, &linha, &rhs);

      double soma = 0.0;
      for (const auto &kv : linha) {
        const long k = kv.first - p0;
        if (k < 0 || k >= static_cast<long>(idx.size())) continue;
        const long vb = idx[static_cast<size_t>(k)].first;
        const long vi = idx[static_cast<size_t>(k)].second;
        if (vb < 0) continue;
        soma += kv.second * exata(m_->centroide(vb, vi));
      }
      // Normaliza pela area, para o numero nao crescer com o refino so' porque as
      // linhas foram multiplicadas por |face|: o que se quer e' o residuo do
      // LAPLACIANO, nao do balanco integrado.
      const double area = m_->area(b, i);
      pior = std::fmax(pior, std::fabs(soma - rhs) / area);
    }
  return pior;
}

double Poisson::u(long id) const
{
  const long k = id - primeiro_;
  if (k < 0 || k >= static_cast<long>(sol_.size()))
    throw std::out_of_range("Poisson::u: id " + std::to_string(id) + " fora da faixa local");
  return sol_[static_cast<size_t>(k)];
}

double Poisson::erro_inf(const std::function<double(const Ponto2 &)> &exata) const
{
  double pior = 0.0;
  for (long b = 0; b < m_->num_blocos(); b++)
    for (long i = 0; i < m_->num_folhas_do_bloco(b); i++) {
      const Ponto2 c = m_->centroide(b, i);
      pior = std::fmax(pior, std::fabs(u(m_->id_global(b, i)) - exata(c)));
    }
  return pior;
}

double Poisson::erro_l2(const std::function<double(const Ponto2 &)> &exata) const
{
  double soma = 0.0, area_total = 0.0;
  for (long b = 0; b < m_->num_blocos(); b++)
    for (long i = 0; i < m_->num_folhas_do_bloco(b); i++) {
      const Ponto2 c = m_->centroide(b, i);
      const double e = u(m_->id_global(b, i)) - exata(c);
      const double a = m_->area(b, i);
      soma += e * e * a;
      area_total += a;
    }
  return (area_total > 0.0) ? std::sqrt(soma / area_total) : 0.0;
}

}  // namespace t8flow2d
