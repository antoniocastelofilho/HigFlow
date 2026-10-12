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

#include "estencil2d.hxx"

#include <cmath>
#include <map>
#include <memory>
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
          // O CENTRO DO VIZINHO VEM PELO INDICE, que cobre as folhas locais e as
          // de FRANJA.  A versao anterior fazia `vz[k] - primeiro_id()` e
          // `continue` quando saia da faixa local -- e em paralelo isso OMITIA a
          // correcao transversal em toda transicao 2:1 que cruzasse fronteira de
          // particao.
          //
          // E omitia em silencio, de um jeito especialmente difícil de notar: em
          // malha UNIFORME a correcao e' zero, e portanto np=2 e np=3 davam erro
          // de maquina do mesmo jeito.  So' a malha adaptada em paralelo
          // acusaria, e por um erro que parece de truncamento.
          const long vx = m_->indice_do_id(vz[k]);
          if (vx < 0)
            throw std::runtime_error(
                "Poisson::monta_linha: vizinho " + std::to_string(vz[k]) +
                " nao esta' nas folhas locais nem na franja -- a correcao "
                "transversal da transicao 2:1 seria OMITIDA em silencio");

          const Ponto2 cv = m_->centroide_do_indice(vx);
          const double delta = cv[t] - c[t];       // deslocamento transversal
          if (std::fabs(delta) < 1e-14) continue;  // alinhados: nada a corrigir

          // Quem esta' deslocado e' o lado GROSSO em relacao a' sub-face: o
          // gradiente transversal e' estimado na celula GROSSA (a de h maior),
          // e o deslocamento aplicado tem o sinal que leva o valor dela para a
          // linha do centro da fina.
          const bool self_grosso = (h_viz < hd);

          // A CELULA GROSSA E' IDENTIFICADA POR INDICE, e pode ser a de FRANJA.
          // A versao anterior guardava o par (arvore, folha) dela, que a folha de
          // franja nao tem -- e em paralelo, com o vizinho grosso do outro lado da
          // fronteira de particao, o par vinha de um indice fora da faixa e
          // apontava para outra celula.
          //
          // A franja precisa ser EXPANSIVEL aqui, porque o gradiente transversal
          // consulta os vizinhos da celula grossa.  Com `aneis + 1` camadas e um
          // anel de fluxo, a camada 1 e' expansivel e isto nunca falha; a guarda
          // existe para que, se algum dia falhar, falhe dizendo o que falta.
          const long gx = self_grosso ? m_->indice_do_id(m_->id_global(b, i)) : vx;
          if (gx < 0 || !m_->expansivel_do_indice(gx))
            throw std::runtime_error(
                "Poisson::monta_linha: a celula grossa da transicao 2:1 esta' na "
                "borda da franja e nao pode ser expandida -- o gradiente "
                "transversal precisa dos vizinhos dela.  Faltam camadas de franja");
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
          const Ponto2 cg = m_->centroide_do_indice(gx);
          const Ponto2 hg = m_->lados_do_indice(gx);
          double dist_total = 0.0;
          bool ok = true;

          for (int lt = 0; lt < 2 && ok; lt++) {
            const int fat = 2 * t + lt;
            const double sgn = (lt == 0) ? -1.0 : 1.0;   // menos / mais
            const std::vector<long> vt = m_->vizinhos_do_indice(gx, fat);
            if (!vt.empty()) {
              // O NIVEL DO VIZINHO TRANSVERSAL VIAJA JUNTO COM O ID.
              //
              // A versao anterior o buscava por `indice_do_id` e, quando o
              // vizinho nao estava na tabela, caia num FALLBACK silencioso -- o
              // nivel da propria celula grossa.  Isso dava `ht` errado, `dt`
              // errado, gradiente errado: medido, `grau1` em malha adaptada
              // passou de 3,97e-13 em np=1 para 1,57e-01 em np=2, num esquema que
              // e' EXATO para campo linear.  O fallback plausivel foi pior que
              // uma recusa.
              const std::vector<int> nvt = m_->niveis_vizinhos_do_indice(gx, fat);
              if (nvt.empty() || nvt[0] < 0)
                throw std::runtime_error(
                    "Poisson::monta_linha: nivel do vizinho transversal "
                    "desconhecido -- o gradiente transversal da transicao 2:1 "
                    "seria calculado com a distancia errada");
              const int dnt = m_->nivel_do_indice(gx) - nvt[0];
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
              grad[m_->id_global_do_indice(gx)] += sgn;
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

  // O MLS com ordem 3 e suporte largo da' linhas MUITO mais cheias: ~40 pontos de
  // suporte por ponto de quadratura, oito pontos de quadratura por celula.  A
  // pre-alocacao tem de acompanhar, senao o PETSc realoca a cada insercao.
  const int nnz = (esquema_ == Esquema::MlsGauss) ? 200 : 25;
  Solver s(motor, static_cast<int>(primeiro_), static_cast<std::size_t>(n_), nnz);
  s.tolerancia(tol);
  s.max_iteracoes(max_it);

  const bool mlsg = (esquema_ == Esquema::MlsGauss);
  const bool mlsp = (esquema_ == Esquema::MlsPontual);
  const bool compacto = (esquema_ == Esquema::Compacto4 ||
                         esquema_ == Esquema::Compacto4Lele);

  // O ESQUEMA DE FLUXO TAMBEM PRECISA DE FRANJA, e de DUAS camadas, nao uma.
  //
  // A linha dele toca o vizinho (camada 1) e, na transicao 2:1, consulta os
  // vizinhos da celula GROSSA para montar o gradiente transversal -- que, se a
  // grossa for a de franja, estao na camada 2.  Com uma camada so', essa consulta
  // devolveria vazio e a correcao transversal seria OMITIDA: em malha uniforme a
  // correcao e' zero e nada apareceria, e em malha adaptada o erro voltaria a ser
  // o de um esquema sem correcao -- medido em 3,9e-02 contra o exato.
  // AS CAMADAS DE FRANJA QUE CADA FAMILIA PEDE.  O laplaciano pontual consulta o
  // suporte no CENTRO da celula, e nao sobre a faceta -- logo a celula que
  // contem o ponto e' sempre LOCAL, e `aneis` camadas bastam.  O `MlsGauss`
  // precisa de uma a mais, porque o ponto de quadratura esta' sobre a faceta e a
  // celula que o contem pode ser de franja.
  if (!compacto)
    const_cast<Malha *>(m_)->exige_camadas(mlsg   ? mls_aneis_ + 1
                                           : mlsp ? mls_aneis_
                                                  : 2);

  long ng = 0;
  std::vector<long> grade;
  if (compacto) grade = grade_logica(&ng);
  const std::pair<const std::vector<long> *, long> gp{&grade, ng};

  std::unique_ptr<ConstrutorEstencil> ce;
  if (mlsg || mlsp) {
    ce.reset(new ConstrutorEstencil(*m_, mls_ordem_, mls_maxpts_, mls_aneis_));
    ce->contorno(contorno_);
    ce->mls().zera_contadores();
    ce->zera_contadores_de_franja();

    // A FRANJA TEM DE COBRIR OS ANEIS QUE O SUPORTE VAI PEDIR.
    //
    // Sem isso, em paralelo, o ponto de quadratura sobre uma faceta de fronteira
    // de particao cai numa celula de franja cujos vizinhos este processo nao
    // conhece, e o suporte sai com UM ponto -- contra os 10 que a ordem 3 exige.
    // Medido: `suporte de 1 pontos em 6 aneis, minimo 10 para ordem 3`.
    //
    // SAO `aneis + 1` CAMADAS, e nao `aneis`, e a diferenca de um foi medida: o
    // ponto de quadratura esta' SOBRE a faceta, logo a celula que o contem pode
    // ser de franja (camada 1) -- e a partir dela o suporte ainda precisa abrir
    // os `aneis` aneis inteiros.  Com `aneis` camadas, np=2 truncou 32 suportes e
    // np=3 falhou com `suporte de 9 pontos em 2 aneis, minimo 10`.
    //
    // Em np=1 nao ha' franja e a chamada nao cobra nada.
    const_cast<Malha *>(m_)->exige_camadas(mls_aneis_ + 1);
  }

  if (mlsg || mlsp) {
    for (long b = 0; b < m_->num_blocos(); b++)
      for (long i = 0; i < m_->num_folhas_do_bloco(b); i++) {
        std::map<long, double> linha;
        double rhs = 0.0;
        if (mlsp) monta_linha_mls_pontual(b, i, ce.get(), &linha, &rhs);
        else monta_linha_mls(b, i, ce.get(), &linha, &rhs);
        std::vector<int> cols;
        std::vector<double> vals;
        for (const auto &kv : linha) {
          cols.push_back(static_cast<int>(kv.first));
          vals.push_back(kv.second);
        }
        const long id = m_->id_global(b, i);
        s.linha(static_cast<int>(id), cols, vals);
        s.lado_direito(static_cast<int>(id), rhs);
      }
    ajustes_mls_ = ce->mls().ajustes();
    rebaixos_mls_ = ce->mls().rebaixamentos();
    truncados_franja_ = reduz_soma(static_cast<double>(ce->truncados_pela_franja()));
    s.monta();
    s.resolve();
    convergiu_ = s.convergiu();
    iteracoes_ = s.iteracoes();
    residuo_ = s.residuo();
    sol_ = s.solucao();
    return;
  }

  if (compacto) {
    for (long iy = 0; iy < ng; iy++)
      for (long ix = 0; ix < ng; ix++) {
        std::map<long, double> linha;
        double rhs = 0.0;
        monta_linha_compacta(ix, iy, &gp, &linha, &rhs);
        std::vector<int> cols;
        std::vector<double> vals;
        for (const auto &kv : linha) {
          cols.push_back(static_cast<int>(kv.first));
          vals.push_back(kv.second);
        }
        const long id = grade[static_cast<size_t>(iy) * ng + ix];
        s.linha(static_cast<int>(id), cols, vals);
        s.lado_direito(static_cast<int>(id), rhs);
      }
    s.monta();
    s.resolve();
    convergiu_ = s.convergiu();
    iteracoes_ = s.iteracoes();
    residuo_ = s.residuo();
    sol_ = s.solucao();
    return;
  }

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

  const bool mlsg = (esquema_ == Esquema::MlsGauss);
  const bool mlsp = (esquema_ == Esquema::MlsPontual);
  const bool compacto = (esquema_ == Esquema::Compacto4 ||
                         esquema_ == Esquema::Compacto4Lele);

  // O ESQUEMA DE FLUXO TAMBEM PRECISA DE FRANJA, e de DUAS camadas, nao uma.
  //
  // A linha dele toca o vizinho (camada 1) e, na transicao 2:1, consulta os
  // vizinhos da celula GROSSA para montar o gradiente transversal -- que, se a
  // grossa for a de franja, estao na camada 2.  Com uma camada so', essa consulta
  // devolveria vazio e a correcao transversal seria OMITIDA: em malha uniforme a
  // correcao e' zero e nada apareceria, e em malha adaptada o erro voltaria a ser
  // o de um esquema sem correcao -- medido em 3,9e-02 contra o exato.
  // AS CAMADAS DE FRANJA QUE CADA FAMILIA PEDE.  O laplaciano pontual consulta o
  // suporte no CENTRO da celula, e nao sobre a faceta -- logo a celula que
  // contem o ponto e' sempre LOCAL, e `aneis` camadas bastam.  O `MlsGauss`
  // precisa de uma a mais, porque o ponto de quadratura esta' sobre a faceta e a
  // celula que o contem pode ser de franja.
  if (!compacto)
    const_cast<Malha *>(m_)->exige_camadas(mlsg   ? mls_aneis_ + 1
                                           : mlsp ? mls_aneis_
                                                  : 2);

  long ng = 0;
  std::vector<long> grade;
  if (compacto) grade = grade_logica(&ng);
  const std::pair<const std::vector<long> *, long> gp{&grade, ng};

  std::unique_ptr<ConstrutorEstencil> ce;
  if (mlsg || mlsp) {
    ce.reset(new ConstrutorEstencil(*m_, mls_ordem_, mls_maxpts_, mls_aneis_));
    ce->contorno(contorno_);
    ce->mls().zera_contadores();
    ce->zera_contadores_de_franja();
    const_cast<Malha *>(m_)->exige_camadas(mls_aneis_ + 1);
  }

  double pior = 0.0;
  for (long b = 0; b < m_->num_blocos(); b++)
    for (long i = 0; i < m_->num_folhas_do_bloco(b); i++) {
      std::map<long, double> linha;
      double rhs = 0.0;
      if (mlsp) {
        monta_linha_mls_pontual(b, i, ce.get(), &linha, &rhs);
      } else if (mlsg) {
        monta_linha_mls(b, i, ce.get(), &linha, &rhs);
      } else if (compacto) {
        // A celula e' localizada na grade pela GEOMETRIA, para o residuo medir o
        // MESMO operador que o solver monta.
        const Ponto2 cc = m_->centroide(b, i);
        const Bloco &kk = m_->bloco(0);
        const double hh = (kk.hi[0] - kk.lo[0]) / ng;
        const long ix = static_cast<long>((cc[0] - kk.lo[0]) / hh);
        const long iy = static_cast<long>((cc[1] - kk.lo[1]) / hh);
        monta_linha_compacta(ix, iy, &gp, &linha, &rhs);
      } else {
        monta_linha(b, i, &idx, &linha, &rhs);
      }

      // A SOMA PERCORRE A LINHA INTEIRA, inclusive as colunas que caem na FRANJA.
      //
      // A versao anterior fazia `kv.first - primeiro_id()` e DESCARTAVA o que
      // saisse da faixa local.  Em np=1 nada sai, e o descarte nunca acontecia.
      // Em paralelo, cada celula vizinha a uma fronteira de particao perderia os
      // termos de la' -- e o residuo sairia MENOR, isto e', um operador
      // inconsistente se apresentando como consistente.  Era o oraculo mentindo
      // a favor do codigo, que e' o pior modo de falha possivel aqui.
      double soma = 0.0;
      bool fora = false;
      for (const auto &kv : linha) {
        const long vx = m_->indice_do_id(kv.first);
        if (vx < 0) { fora = true; continue; }
        soma += kv.second * exata(m_->centroide_do_indice(vx));
      }
      if (fora)
        throw std::runtime_error(
            "Poisson::residuo_da_exata: a linha referencia celula que nao esta' "
            "nas folhas locais nem na franja -- o residuo seria calculado com "
            "termos FALTANDO, e sairia menor que o verdadeiro");
      // Normaliza pela area no esquema de FLUXO, cujas linhas foram multiplicadas
      // por |face|; o compacto ja' esta' na forma do laplaciano pontual, e dividir
      // ali daria um numero que cresce com o refino sem significado.
      // A NORMALIZACAO E' DA FAMILIA, e quem a conhece e' a `Discretizacao`:
      // volumes finitos dividem pela area (a linha foi multiplicada por `|face|`),
      // diferencas finitas nao dividem (a linha ja' esta' na forma pontual).
      // Dividir pela area em DF daria um numero que CRESCE com o refino sem
      // significado -- e isso ja' aconteceu neste projeto.
      const double norm = discretizacao().normalizacao_do_residuo(m_->area(b, i));
      pior = std::fmax(pior, std::fabs(soma - rhs) / norm);
    }
  return reduz_max(pior);
}

// ---------------------------------------------------------------------------
// RECONSTRUCAO MLS COM QUADRATURA DE GAUSS
// ---------------------------------------------------------------------------

void Poisson::monta_linha_mls(long b, long i, const void *construtor,
                              std::map<long, double> *linha_out, double *rhs_out) const
{
  const ConstrutorEstencil &ce = *static_cast<const ConstrutorEstencil *>(construtor);

  const Ponto2 c = m_->centroide(b, i);
  const Ponto2 h = m_->lados_da_folha(b, i);
  const double area = h[0] * h[1];

  std::map<long, double> linha;

  // GAUSS DE DOIS PONTOS no intervalo [-1,1]: nos em +-1/sqrt(3), pesos 1 cada.
  // Normalizado para a MEDIA sobre a faceta, cada peso vira 1/2.  Exato para
  // integrando CUBICO -- contra o ponto medio, exato so' para linear.
  const double no = 1.0 / std::sqrt(3.0);
  const double wg = 0.5;

  // O LADO DIREITO TAMBEM E' QUADRATURA, e ignorar isso foi o que travou o esquema
  // em SEGUNDA ordem com o fluxo ja' exato: `f(c) * area` e' a regra do ponto
  // medio, exata so' para `f` linear, e portanto um erro `O(h^2)` relativo em
  // `integral f dA`.  No esquema de fluxo de segunda ordem esse erro se perdia no
  // erro do proprio fluxo; aqui ele passa a ser o termo DOMINANTE.
  //
  // Gauss 2x2 e' exato para `f` bicubico -- que cobre `f = lap(u)` de qualquer `u`
  // de grau <= 5 -- e deixa o erro do lado direito em `O(h^4)` relativo.
  double rhs = 0.0;
  if (quad_fonte_ == QuadFonte::PontoMedio) {
    rhs = f_(c);
  } else {
    for (int qy = 0; qy < 2; qy++)
      for (int qx = 0; qx < 2; qx++) {
        const Ponto2 xg{c[0] + (qx == 0 ? -no : no) * h[0] / 2.0,
                        c[1] + (qy == 0 ? -no : no) * h[1] / 2.0};
        rhs += wg * wg * f_(xg);
      }
  }
  rhs *= area;

  for (int d = 0; d < 2; d++) {
    const int t = 1 - d;
    for (int lado = 0; lado < 2; lado++) {
      const int fa = 2 * d + lado;
      const double sinal = (lado == 0) ? -1.0 : 1.0;   // normal exterior

      // AS SUB-FACETAS.  Em transicao 2:1 a faceta da celula grossa e' coberta por
      // DUAS sub-facetas -- cada uma com o tamanho da face da vizinha fina -- e o
      // fluxo total e' a soma.  Integrar a face inteira de uma vez daria o fluxo
      // medio sobre ela, que NAO e' a soma dos dois fluxos finos quando a derivada
      // varia ao longo dela.
      const std::vector<int> nvz = m_->niveis_vizinhos(b, i, fa);
      const int nivel_self = m_->nivel(b, i);
      int nsub = 1;
      double larg_sub = h[t];
      if (!nvz.empty() && nvz[0] > nivel_self) {
        nsub = 1 << (nvz[0] - nivel_self);        // 2 em 2:1
        larg_sub = h[t] / nsub;
      }

      for (int sub = 0; sub < nsub; sub++) {
        // centro desta sub-faceta
        Ponto2 cf = c;
        cf[d] += sinal * h[d] / 2.0;
        cf[t] += -h[t] / 2.0 + (sub + 0.5) * larg_sub;

        for (int q = 0; q < 2; q++) {
          Ponto2 xg = cf;
          xg[t] += (q == 0 ? -no : no) * larg_sub / 2.0;

          // O ESTENCIL DA DERIVADA NO PONTO DE QUADRATURA.  O suporte depende so'
          // do PONTO, logo a celula vizinha obtem o mesmo estencil para o mesmo
          // ponto -- e' daqui que vem a conservacao.
          const Estencil e = ce.derivada_em(xg, d);

          const double peso = sinal * wg * larg_sub;
          for (size_t k = 0; k < e.ids.size(); k++)
            linha[e.ids[k]] += peso * e.pesos[k];
          // o termo conhecido do fechamento vai para o lado direito
          rhs -= peso * e.rhs;
        }
      }
    }
  }

  *linha_out = linha;
  *rhs_out = rhs;
}

void Poisson::monta_linha_mls_pontual(long b, long i, const void *construtor,
                                      std::map<long, double> *linha_out,
                                      double *rhs_out) const
{
  const ConstrutorEstencil &ce = *static_cast<const ConstrutorEstencil *>(construtor);
  const Ponto2 c = m_->centroide(b, i);

  // A LINHA E' O LAPLACIANO NO PONTO, sem quadratura e sem area.  E' o que
  // diferencas finitas significa: a equacao e' imposta no centro, nao integrada
  // sobre a celula.
  //
  // E A FONTE E' `f(c)`, tambem sem integrar.  Misturar as duas coisas -- linha
  // pontual com fonte integrada, ou o contrario -- nao falharia: daria um
  // operador consistente com outra equacao, que converge para a coisa errada.
  // E' por isso que a escolha vive na `Discretizacao`, onde as duas decisoes sao
  // a MESMA decisao.
  const Estencil e = ce.laplaciano_em(c);

  std::map<long, double> linha;
  for (size_t k = 0; k < e.ids.size(); k++) linha[e.ids[k]] += e.pesos[k];

  *linha_out = linha;
  *rhs_out = f_(c) - e.rhs;
}

// ---------------------------------------------------------------------------
// O ESQUEMA COMPACTO DE QUARTA ORDEM
// ---------------------------------------------------------------------------

std::vector<long> Poisson::grade_logica(long *n) const
{
  // SO' MALHA UNIFORME DE UM BLOCO.  O estencil compacto usa vizinhos DIAGONAIS, e
  // em transicao 2:1 o vizinho diagonal pode nem existir como celula unica -- nao
  // e' limitacao de implementacao, e' o esquema que nao esta' definido ali sem uma
  // reconstrucao propria.  Recusar e' o que impede um resultado plausivel.
  if (m_->num_blocos() != 1)
    throw std::runtime_error("Poisson: o esquema compacto exige UM bloco; a malha tem " +
                             std::to_string(m_->num_blocos()));
  if (m_->nivel_minimo() != m_->nivel_maximo())
    throw std::runtime_error(
        "Poisson: o esquema compacto exige malha UNIFORME; esta tem niveis " +
        std::to_string(m_->nivel_minimo()) + " a " + std::to_string(m_->nivel_maximo()) +
        ".  O estencil usa vizinhos diagonais, e na transicao 2:1 eles nao estao"
        " definidos sem reconstrucao propria");

  const Bloco &k = m_->bloco(0);
  const int div = 1 << m_->nivel(0, 0);
  *n = div;

  // A grade e' montada pela GEOMETRIA (busca por ponto), nao pela ordem de
  // percurso: a ordem da curva de preenchimento nao e' lexicografica, e supor que
  // fosse daria uma grade embaralhada que ainda "funciona" -- com a solucao
  // permutada, que e' erro silencioso.
  std::vector<long> g(static_cast<size_t>(div) * static_cast<size_t>(div), -1);
  const double hx = (k.hi[0] - k.lo[0]) / div;
  const double hy = (k.hi[1] - k.lo[1]) / div;
  for (long iy = 0; iy < div; iy++)
    for (long ix = 0; ix < div; ix++) {
      const Ponto2 c = {k.lo[0] + (ix + 0.5) * hx, k.lo[1] + (iy + 0.5) * hy};
      const long id = m_->id_do_ponto(c);
      if (id < 0)
        throw std::runtime_error("Poisson: a grade logica nao fechou em (" +
                                 std::to_string(ix) + ", " + std::to_string(iy) + ")");
      g[static_cast<size_t>(iy) * div + ix] = id;
    }
  return g;
}

namespace {

//! Expressao linear nas incognitas: `sum coef*u[id] + cte`.
struct Expr {
  std::map<long, double> coef;
  double cte = 0.0;
};

//! O valor da celula logica `(ix, iy)`, como expressao -- com FANTASMA quando o
//! indice sai do dominio.
//
// O fantasma vem da interpolacao de QUATRO pontos que reproduz o valor de
// Dirichlet na parede: com a parede em `s=0` e os centros em `h/2, 3h/2, 5h/2`, o
// fantasma em `-h/2` satisfaz
//
//     g = (5/16) u_f + (15/16) u_0 - (5/16) u_1 + (1/16) u_2
//     u_f = (16/5) g - 3 u_0 + u_1 - (1/5) u_2
//
// (a soma dos coeficientes de Lagrange e' 1, conferido: 5/16+15/16-5/16+1/16 = 1.)
//
// QUATRO PONTOS E NAO DOIS, e a razao e' a ordem: o fantasma de dois pontos
// (`u_f = 2g - u_0`) e' exato para campo LINEAR e limitaria o esquema a segunda
// ordem na primeira celula -- jogando fora a quarta ordem do interior.  Com quatro
// pontos ele e' exato para grau 3, que e' o que o compacto precisa.
//
// O CANTO sai por COMPOSICAO: o fantasma diagonal e' o fantasma em `y` aplicado
// sobre a coluna de fantasmas em `x`.  A recursao faz isso sozinha, e e' por isso
// que a funcao chama a si mesma em vez de tratar o canto como caso especial --
// caso especial de canto e' onde o sinal se perde.
Expr valor(long ix, long iy, long n, const std::vector<long> &g,
           const Malha &m, const ContornoPorRegiao &bc, int prof = 0)
{
  Expr e;
  if (prof > 2) return e;                      // guarda contra recursao sem fim

  if (ix >= 0 && ix < n && iy >= 0 && iy < n) {
    e.coef[g[static_cast<size_t>(iy) * n + ix]] = 1.0;
    return e;
  }

  // Decide a direcao do fantasma: a que saiu do dominio (se as duas sairam, trata
  // `x` primeiro e deixa a recursao cuidar de `y`).
  const bool fora_x = (ix < 0 || ix >= n);
  const long passo = fora_x ? ((ix < 0) ? 1 : -1) : ((iy < 0) ? 1 : -1);

  // O ponto da PAREDE, para a condicao ser avaliada no lugar certo.
  const Bloco &k = m.bloco(0);
  const double hx = (k.hi[0] - k.lo[0]) / n;
  const double hy = (k.hi[1] - k.lo[1]) / n;
  const long jx = fora_x ? ((ix < 0) ? 0 : n - 1) : ix;
  const long jy = fora_x ? iy : ((iy < 0) ? 0 : n - 1);
  Ponto2 pw = {k.lo[0] + (jx + 0.5) * hx, k.lo[1] + (jy + 0.5) * hy};
  if (fora_x) pw[0] += (ix < 0 ? -hx / 2.0 : hx / 2.0);
  else        pw[1] += (iy < 0 ? -hy / 2.0 : hy / 2.0);

  const Contorno &cond = bc.em(pw);
  if (cond.tipo != TipoContorno::Dirichlet)
    throw std::runtime_error(
        "Poisson: o esquema compacto monta so' DIRICHLET; a parede em (" +
        std::to_string(pw[0]) + ", " + std::to_string(pw[1]) + ") e' outra condicao."
        "  Neumann e Robin no compacto precisam de um fechamento proprio, de ordem"
        " compativel -- montar com o de Dirichlet daria quarta ordem no interior e"
        " segunda na parede, sem aviso");

  e.cte = (16.0 / 5.0) * cond.g(pw);
  const double peso[3] = {-3.0, 1.0, -1.0 / 5.0};
  for (int r = 0; r < 3; r++) {
    const long qx = fora_x ? (jx + passo * r) : ix;
    const long qy = fora_x ? iy : (jy + passo * r);
    const Expr sub = valor(qx, qy, n, g, m, bc, prof + 1);
    for (const auto &kv : sub.coef) e.coef[kv.first] += peso[r] * kv.second;
    e.cte += peso[r] * sub.cte;
  }
  return e;
}

}  // namespace

void Poisson::monta_linha_compacta(long ix, long iy, const void *grade,
                                   std::map<long, double> *linha_out,
                                   double *rhs_out) const
{
  const std::pair<const std::vector<long> *, long> *gp =
      static_cast<const std::pair<const std::vector<long> *, long> *>(grade);
  const std::vector<long> &g = *gp->first;
  const long n = gp->second;

  const Bloco &k = m_->bloco(0);
  const double h = (k.hi[0] - k.lo[0]) / n;
  const Ponto2 c = {k.lo[0] + (ix + 0.5) * h, k.lo[1] + (iy + 0.5) * h};

  // O ESTENCIL DE NOVE PONTOS.  Na forma de Mehrstellen,
  //
  //     (1/(6h^2)) [ 4(faces) + 1(diagonais) - 20(centro) ]  =  f + (h^2/12) lap f
  //
  // A soma dos coeficientes e' 4*4 + 4*1 - 20 = 0, como tem de ser para o
  // laplaciano de uma constante dar zero -- a conferencia que pega escala e sinal.
  const double esc = 1.0 / (6.0 * h * h);
  const int dx[9] = {0, 1, -1, 0,  0,  1, -1,  1, -1};
  const int dy[9] = {0, 0,  0, 1, -1,  1,  1, -1, -1};
  const double w[9] = {-20.0, 4.0, 4.0, 4.0, 4.0, 1.0, 1.0, 1.0, 1.0};

  std::map<long, double> linha;
  double cte = 0.0;
  for (int t = 0; t < 9; t++) {
    const Expr e = valor(ix + dx[t], iy + dy[t], n, g, *m_, contorno_);
    for (const auto &kv : e.coef) linha[kv.first] += esc * w[t] * kv.second;
    cte += esc * w[t] * e.cte;
  }

  // O LADO DIREITO, nas duas variantes.
  double rhs = 0.0;
  if (esquema_ == Esquema::Compacto4) {
    // MEHRSTELLEN: exige `lap(f)` analitico.  Sem ele, rodar com zero daria
    // SEGUNDA ordem e pareceria quarta falhando -- por isso recusa.
    if (!lap_f_)
      throw std::invalid_argument(
          "Poisson: o esquema Compacto4 (Mehrstellen) exige `fonte_laplaciano`."
          "  Rodar sem ele daria segunda ordem e pareceria quarta falhando;"
          "  use Compacto4Lele, que usa apenas `f`");
    rhs = f_(c) + (h * h / 12.0) * lap_f_(c);
  } else {
    // LELE: `A_x A_y f / (6/5)^2`, isto e' `f` nos NOVE pontos com os pesos do
    // operador tridiagonal nas duas direcoes -- sem `lap(f)`.
    //
    // `f` e' avaliado inclusive FORA do dominio, e isso e' legitimo: ele e' uma
    // funcao conhecida, nao um campo discreto.  Interpolar `f` ali em vez de
    // avalia-lo introduziria erro onde nao ha' necessidade.
    const double a[3] = {1.0 / 10.0, 1.0, 1.0 / 10.0};
    double soma = 0.0;
    for (int jx = -1; jx <= 1; jx++)
      for (int jy = -1; jy <= 1; jy++) {
        const Ponto2 q = {c[0] + jx * h, c[1] + jy * h};
        soma += a[jx + 1] * a[jy + 1] * f_(q);
      }
    rhs = soma / (1.44);                        // (6/5)^2
  }

  // O termo conhecido dos fantasmas passa para o lado direito com sinal trocado.
  *linha_out = linha;
  *rhs_out = rhs - cte;
}

double Poisson::u(long id) const
{
  const long k = id - primeiro_;
  if (k < 0 || k >= static_cast<long>(sol_.size()))
    throw std::out_of_range("Poisson::u: id " + std::to_string(id) + " fora da faixa local");
  return sol_[static_cast<size_t>(k)];
}

// ---------------------------------------------------------------------------
// AS REDUCOES
// ---------------------------------------------------------------------------
//
// Ficam aqui, e nao na `Malha`, porque sao propriedades das NORMAS e do
// residuo: a malha nao tem opiniao sobre como se combina um erro.  O
// comunicador vem dela.

void Poisson::discretizacao(const Discretizacao &d)
{
  switch (d.variante()) {
    case Variante::FluxoBalanco:        esquema_ = Esquema::Fluxo2; break;
    case Variante::MlsGauss:            esquema_ = Esquema::MlsGauss; break;
    case Variante::CompactoMehrstellen: esquema_ = Esquema::Compacto4; break;
    case Variante::CompactoLele:        esquema_ = Esquema::Compacto4Lele; break;
    case Variante::MlsPontual:          esquema_ = Esquema::MlsPontual; break;
  }
  if (d.ordem_mls() > 0) mls(d.ordem_mls(), mls_maxpts_);
}

Discretizacao Poisson::discretizacao() const
{
  switch (esquema_) {
    case Esquema::Fluxo2:         return Discretizacao(Variante::FluxoBalanco);
    case Esquema::MlsGauss:       return Discretizacao(Variante::MlsGauss, mls_ordem_);
    case Esquema::Compacto4:      return Discretizacao(Variante::CompactoMehrstellen);
    case Esquema::Compacto4Lele:  return Discretizacao(Variante::CompactoLele);
    case Esquema::MlsPontual:     return Discretizacao(Variante::MlsPontual, mls_ordem_);
  }
  return Discretizacao::padrao();
}

std::map<long, double> Poisson::linha_montada(long b, long i, double *rhs) const
{
  if (!f_) throw std::invalid_argument("Poisson::linha_montada: fonte nao instalada");
  const std::vector<std::pair<long, long>> idx = monta_indice(*m_);
  std::map<long, double> linha;
  double r = 0.0;
  if (esquema_ == Esquema::MlsGauss || esquema_ == Esquema::MlsPontual) {
    ConstrutorEstencil ce(*m_, mls_ordem_, mls_maxpts_, mls_aneis_);
    ce.contorno(contorno_);
    if (esquema_ == Esquema::MlsPontual) monta_linha_mls_pontual(b, i, &ce, &linha, &r);
    else monta_linha_mls(b, i, &ce, &linha, &r);
  } else if (esquema_ == Esquema::Compacto4 || esquema_ == Esquema::Compacto4Lele) {
    long ng = 0;
    const std::vector<long> grade = grade_logica(&ng);
    const std::pair<const std::vector<long> *, long> gp{&grade, ng};
    const Ponto2 cc = m_->centroide(b, i);
    const Bloco &kk = m_->bloco(0);
    const double hh = (kk.hi[0] - kk.lo[0]) / ng;
    monta_linha_compacta(static_cast<long>((cc[0] - kk.lo[0]) / hh),
                         static_cast<long>((cc[1] - kk.lo[1]) / hh), &gp, &linha, &r);
  } else {
    monta_linha(b, i, &idx, &linha, &r);
  }
  if (rhs != nullptr) *rhs = r;
  return linha;
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

double Poisson::erro_inf(const std::function<double(const Ponto2 &)> &exata) const
{
  double pior = 0.0;
  for (long b = 0; b < m_->num_blocos(); b++)
    for (long i = 0; i < m_->num_folhas_do_bloco(b); i++) {
      const Ponto2 c = m_->centroide(b, i);
      pior = std::fmax(pior, std::fabs(u(m_->id_global(b, i)) - exata(c)));
    }

  // A REDUCAO E' PARTE DA NORMA, nao um detalhe de paralelizacao.  Sem ela cada
  // processo devolve o pior erro DA SUA PARTICAO, e o numero que sai e' menor
  // que o verdadeiro -- o que, num relatorio de verificacao, e' pior que um erro
  // grande: e' um erro grande que se apresenta como pequeno.
  //
  // Em np=1 o `Allreduce` e' identidade e nao cobra nada, e por isso nao ha'
  // ramo serial separado -- um ramo serial seria codigo que o paralelo nunca
  // exercita, que e' a classe de defeito desta frente.
  return reduz_max(pior);
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
  // AS DUAS SOMAS REDUZEM, e o quociente vem DEPOIS: reduzir o quociente daria a
  // media das normas por processo, que nao e' a norma.
  soma = reduz_soma(soma);
  area_total = reduz_soma(area_total);
  return (area_total > 0.0) ? std::sqrt(soma / area_total) : 0.0;
}

}  // namespace t8flow2d
