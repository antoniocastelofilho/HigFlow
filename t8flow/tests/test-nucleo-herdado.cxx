// O PRIMEIRO TIJOLO DO T8FLOW, e ele existe para afirmar as decisoes de projeto
// 1, 2 e 5 de uma vez, numa UNICA unidade de traducao:
//
//   D1  o nucleo independente de malha da HiGTree (solver e MLS) e' reusavel;
//   D2  a estrutura e' a floresta do t8code;
//   D5  tudo compila com o compilador e a norma do t8code (gnu++20).
//
// Por que juntos: compilar as tres coisas separadas nao prova que convivem.  A
// combinacao arriscada era justamente esta -- codigo C da HiGTree, compilado como
// C++20, na mesma unidade que o t8code.  Este arquivo e' o experimento, e os
// casos abaixo afirmam VALOR, nao forma.
//
// O CASO QUE IMPORTA e' o ultimo: o suporte sai da FLORESTA (centroides das
// folhas) e os pesos saem do MLS herdado.  Esse e' o contrato da secao 3 do
// doc/projeto-t8flow.md -- a malha entrega suporte, a discretizacao entrega
// pesos --, exercitado de ponta a ponta.
//
// ORACULO, e a escolha e' deliberada: reproducao polinomial COM CONTRASTE.  Um
// ajuste de ordem k reproduz exatamente polinomio de grau <= k.  Sozinho isso
// nao discrimina quase nada (ver secao 4.2 do projeto: nao pega franja ausente,
// nem parede errada com campo linear, nem suporte deslocado), entao todo caso de
// reproducao vem acompanhado de um campo de grau ACIMA da ordem, cujo erro tem
// de ser ordens de grandeza maior.  Se os dois derem o mesmo erro, nao ha'
// reproducao -- ha' coincidencia, e o teste reprova dizendo isso.

#include "testing.h"

#include <array>
#include <cstdlib>
#include <vector>

// SEM `extern "C"`, e isto e' consequencia direta da decisao 5: as unidades
// herdadas sao compiladas pelo MESMO compilador e na MESMA norma, entao os
// simbolos ja' vem com encadeamento C++ dos dois lados.  Envolver em extern "C"
// aqui daria `undefined reference` no link -- o sintoma classico de duas normas
// convivendo, que e' exatamente o que o t8Flow nao tem.
#include "coord.h"
#include "utils.h"
#include "wls.h"
#include "solver.h"

// A PRIMEIRA FRICCAO REAL DA DECISAO 5, e ela aparece no tipo mais basico que o
// nucleo herdado expoe: `typedef coordtype Point[DIM]` (coord.h:31) e' um ARRAY
// CRU.  Array nao e' copiavel nem construivel in-place, entao `std::vector<Point>`
// NAO COMPILA -- o erro sai dentro do <bits/stl_construct.h> e nao diz isto.
//
// A ponte e' `std::array<real, DIM>`, que a norma garante ter o mesmo layout de
// `real[DIM]`: um vetor dela e' contiguamente um `Point[]` valido, e a conversao
// e' de tipo, nao de dados.  E' o padrao a usar em toda a fronteira do t8Flow
// com o nucleo herdado.
using Pt = std::array<real, DIM>;

static Point *como_pontos(std::vector<Pt> &v)
{
    static_assert(sizeof(Pt) == sizeof(Point), "std::array<real,DIM> tem de ter o layout de Point");
    return reinterpret_cast<Point *>(v.data());
}

#include <petsc.h>

#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_schemes/t8_default/t8_default.hxx>

// ---------------------------------------------------------------------------
// Campos analiticos.  Coeficientes DISTINTOS por direcao de proposito: assim o
// valor identifica a posicao, e um suporte que troque eixos nao passa por acaso.
// ---------------------------------------------------------------------------
static real campo_grau2(const Point x)
{
    real v = 1.0 + 2.0 * x[0] - 3.0 * x[1] + 0.5 * x[0] * x[0] + 0.25 * x[0] * x[1];
#if DIM == 3
    v += 1.5 * x[2] - 0.75 * x[2] * x[2] + 0.125 * x[1] * x[2];
#else
    v += -0.6 * x[1] * x[1];
#endif
    return v;
}

static real campo_grau3(const Point x)
{
    real v = campo_grau2(x) + 4.0 * x[0] * x[0] * x[0] - 2.5 * x[1] * x[1] * x[1];
#if DIM == 3
    v += 3.0 * x[2] * x[2] * x[2];
#endif
    return v;
}

//! Interpola `f` em `x` com os pesos que o MLS produziu para `pts`.
static real interpola(wls_interpolator *wls, int numpts, Point pts[], const Point x,
                      real (*f)(const Point))
{
    std::vector<real> w(numpts > 0 ? numpts : 1, 0.0);
    Point xq;
    for (int d = 0; d < DIM; d++) xq[d] = x[d];
    wls_set_points_and_calc(wls, numpts, pts, xq, w.data());

    real v = 0.0;
    for (int i = 0; i < numpts; i++) v += w[i] * f(pts[i]);
    return v;
}

// ---------------------------------------------------------------------------
// D1: o MLS herdado reproduz polinomio, e o contraste prova que reproduz.
// ---------------------------------------------------------------------------
static void caso_mls_reproduz(void)
{
    const int ORDEM = 2;
    const int lado = 4;                       // 4^DIM pontos numa grade regular
    const real h = 0.25;

    std::vector<Pt> pts;
    for (int i = 0; i < lado; i++)
        for (int j = 0; j < lado; j++)
#if DIM == 3
            for (int k = 0; k < lado; k++)
#endif
            {
                Pt p{};
                p[0] = i * h; p[1] = j * h;
#if DIM == 3
                p[2] = k * h;
#endif
                pts.push_back(p);
            }

    const int numpts = (int) pts.size();
    wls_interpolator *wls = wls_create(DIM, ORDEM, numpts);

    // PONTO DE CONSULTA FORA DO CENTRO DE SIMETRIA, e isto o teste descobriu
    // sozinho: com x no centro exato da nuvem (0,375 nesta grade) o ajuste de
    // ordem 2 reproduziu o campo CUBICO com erro 1,1e-16 -- os termos impares se
    // cancelam contra pesos simetricos, e o contraste deu 0,01x em vez de 1000x.
    // O caso reprovou dizendo "o oraculo nao discrimina", que e' exatamente o seu
    // trabalho.  Coordenadas distintas e assimetricas; ponto interior, para nao
    // medir extrapolacao.
    Point x;
    x[0] = 0.30; x[1] = 0.44;
#if DIM == 3
    x[2] = 0.37;
#endif

    t_case("mls_reproduz_polinomio_do_grau_da_ordem");
    T_CHECK_MSG(numpts >= wls_num_min_points(DIM, ORDEM),
                "nuvem de %d pontos e' menor que o minimo %d para ordem %d",
                numpts, wls_num_min_points(DIM, ORDEM), ORDEM);
    const real obtido2 = interpola(wls, numpts, como_pontos(pts), x, campo_grau2);
    const real exato2  = campo_grau2(x);
    T_NEAR(obtido2, exato2, 1e-9, "campo de grau 2 com ajuste de ordem 2");

    // CONTRASTE: o mesmo ajuste num campo de grau 3 TEM de errar muito mais.  Sem
    // isto, "erro pequeno" nao distingue reproducao de coincidencia.
    t_case("mls_erra_mais_no_grau_acima_da_ordem");
    const real obtido3 = interpola(wls, numpts, como_pontos(pts), x, campo_grau3);
    const real exato3  = campo_grau3(x);
    const real err2 = std::fabs(obtido2 - exato2);
    const real err3 = std::fabs(obtido3 - exato3);
    T_CONTRASTE(err2 > 1e-14 ? err2 : 1e-14, err3, 1e3,
                "grau 2 contra grau 3 no mesmo ajuste de ordem 2");

    wls_destroy(wls);
}

// ---------------------------------------------------------------------------
// D1: o solver herdado resolve um sistema de solucao conhecida.
//
// Tridiagonal de Poisson 1D, -u'' = 2 com u(0)=u(1)=0 discretizada em n+2 nos:
// a solucao do SISTEMA LINEAR e' exatamente x_i = i*(n+1-i)*h^2, e isso e'
// identidade algebrica (nao e' a solucao do problema continuo truncada), entao a
// tolerancia pode ser apertada sem medo.
// ---------------------------------------------------------------------------
static void caso_solver_resolve(void)
{
    const int n = 32;
    const real h = 1.0 / (n + 1);

    solver *s = slv_create(SOLVER_PETSC, 0, (size_t) n);
    slv_set_maxnonzeros(s, 3);
    slv_set_maxiteration(s, 1000);
    slv_set_eps(s, 1e-14);

    for (int i = 0; i < n; i++) {
        int    js[3];
        real   vs[3];
        int    m = 0;
        if (i > 0)     { js[m] = i - 1; vs[m] = -1.0 / (h * h); m++; }
        /* diagonal */ { js[m] = i;     vs[m] =  2.0 / (h * h); m++; }
        if (i < n - 1) { js[m] = i + 1; vs[m] = -1.0 / (h * h); m++; }
        slv_set_Ai(s, i, m, js, vs);
        slv_set_bi(s, i, 2.0);
    }
    // A ORDEM NAO E' LIVRE: o `slv_assemble` do proprio solver faz
    // rhs -> output -> matrix (solver.c:204), e a pre-alocacao da matriz do PETSc
    // consome as linhas acumuladas nessa etapa.  Chamar as tres na ordem inversa
    // custou um segfault dentro do PETSc, sem mensagem que aponte para aqui.
    slv_assemble(s);
    slv_solve(s);

    t_case("solver_resolve_sistema_de_solucao_conhecida");
    T_CHECK_MSG(s->converged,
                "o solver nao convergiu: %d iteracoes, residuo %.3g",
                s->SolverIterations, (double) s->SolverResidual);

    real pior = 0.0;
    int  pior_i = -1;
    for (int i = 0; i < n; i++) {
        const real esperado = (i + 1) * (real) (n - i) * h * h;
        const real e = std::fabs(slv_get_xi(s, i) - esperado);
        if (e > pior) { pior = e; pior_i = i; }
    }
    T_CHECK_MSG(pior <= 1e-10,
                "pior desvio %.3g na linha %d (tolerancia 1e-10)", pior, pior_i);

    slv_destroy(s);
}

// ---------------------------------------------------------------------------
// D2 + D5 + o contrato da secao 3: o SUPORTE vem da floresta, os PESOS vem do
// MLS herdado.  Este e' o caso que nenhuma das duas metades prova sozinha.
// ---------------------------------------------------------------------------
static void caso_suporte_da_floresta(void)
{
    const int NIVEL = 3;                      // 8 por direcao
#if DIM == 2
    const t8_eclass_t ECLASSE = T8_ECLASS_QUAD;
#else
    const t8_eclass_t ECLASSE = T8_ECLASS_HEX;
#endif

    t8_cmesh_t cmesh;
    t8_cmesh_init(&cmesh);                    // obrigatorio ANTES do gerador
    t8_cmesh_new_hypercube(&cmesh, ECLASSE, sc_MPI_COMM_SELF, 0, 0, 0);
    const t8_scheme_c *esquema = t8_scheme_new_default();
    t8_forest_t f = t8_forest_new_uniform(cmesh, esquema, NIVEL, 0, sc_MPI_COMM_SELF);

    // A contagem de folhas e' a primeira afirmacao, e ela e' de VALOR: refino
    // uniforme de nivel L num hipercubo de uma arvore da' (2^L)^DIM folhas.
    long esperadas = 1;
    for (int d = 0; d < DIM; d++) esperadas *= (1L << NIVEL);
    const long folhas = (long) t8_forest_get_local_num_leaf_elements(f);

    t_case("floresta_uniforme_tem_as_folhas_do_nivel");
    T_CHECK_MSG(folhas == esperadas,
                "nivel %d em %dD: esperava %ld folhas, obtive %ld",
                NIVEL, DIM, esperadas, folhas);

    // Os centroides das folhas SAO o suporte.  Em malha uniforme eles formam a
    // grade deslocada de meia celula -- e e' nessa grade que o ajuste trabalha.
    std::vector<Pt> centros;
    const t8_locidx_t narv = t8_forest_get_num_local_trees(f);
    for (t8_locidx_t it = 0; it < narv; it++) {
        const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements(f, it);
        for (t8_locidx_t ie = 0; ie < ne; ie++) {
            const t8_element_t *e = t8_forest_get_leaf_element_in_tree(f, it, ie);
            double c[3] = {0.0, 0.0, 0.0};
            t8_forest_element_centroid(f, it, e, c);
            Pt p{};
            for (int d = 0; d < DIM; d++) p[d] = c[d];
            centros.push_back(p);
        }
    }

    t_case("suporte_colhido_da_floresta_tem_uma_entrada_por_folha");
    T_CHECK_MSG((long) centros.size() == folhas,
                "colhi %ld centroides para %ld folhas", (long) centros.size(), folhas);

    // O suporte LOCAL: as folhas mais proximas do ponto de consulta.  Nao e' a
    // floresta inteira de proposito -- busca sem parada devolve a malha toda, o
    // ajuste reproduz mesmo assim, e o teste passa enquanto a topologia deixa de
    // importar (armadilha registrada na secao 4.2 do projeto).
    const int ORDEM = 2;
    const int MAXPTS = 4 * wls_num_min_points(DIM, ORDEM);
    // Mesma razao do caso anterior: 0,5 e' o centro do hipercubo unitario, e a
    // nuvem de centroides e' simetrica em torno dele.
    Point x;
    x[0] = 0.41; x[1] = 0.56;
#if DIM == 3
    x[2] = 0.47;
#endif

    std::vector<int> ordem(centros.size());
    for (size_t i = 0; i < centros.size(); i++) ordem[i] = (int) i;
    // ordenacao parcial por distancia ao ponto de consulta
    for (int k = 0; k < MAXPTS && k < (int) ordem.size(); k++) {
        int melhor = k;
        for (int i = k + 1; i < (int) ordem.size(); i++)
            if (co_distance(centros[ordem[i]].data(), x) < co_distance(centros[ordem[melhor]].data(), x))
                melhor = i;
        const int t = ordem[k]; ordem[k] = ordem[melhor]; ordem[melhor] = t;
    }

    const int numpts = (MAXPTS < (int) centros.size()) ? MAXPTS : (int) centros.size();
    std::vector<Pt> suporte(numpts);
    for (int i = 0; i < numpts; i++)
        for (int d = 0; d < DIM; d++) suporte[i][d] = centros[ordem[i]][d];

    wls_interpolator *wls = wls_create(DIM, ORDEM, numpts);

    t_case("mls_reproduz_polinomio_sobre_suporte_da_floresta");
    const real obtido2 = interpola(wls, numpts, como_pontos(suporte), x, campo_grau2);
    const real exato2  = campo_grau2(x);
    T_NEAR(obtido2, exato2, 1e-9,
           "campo de grau 2 nos centroides das folhas do t8code");

    t_case("suporte_da_floresta_erra_mais_no_grau_acima");
    const real obtido3 = interpola(wls, numpts, como_pontos(suporte), x, campo_grau3);
    const real err2 = std::fabs(obtido2 - exato2);
    const real err3 = std::fabs(obtido3 - campo_grau3(x));
    T_CONTRASTE(err2 > 1e-14 ? err2 : 1e-14, err3, 1e3,
                "grau 2 contra grau 3 sobre o suporte da floresta");

    wls_destroy(wls);
    t8_forest_unref(&f);
}

// A ORDEM DE ARRANQUE NAO E' LIVRE, e descobri isto com um SIGSEGV dentro do
// PETSc, sem uma linha que apontasse para ca'.  O mecanismo, lido do codigo:
//
//   `_slv_create_petsc` chama `_try_initialize_petsc()` (utils.c:34), que
//   desreferencia um `static char **_argv` -- e QUEM preenche esse estatico e'
//   o `higtree_initialize`.  Sem ele, `*_argv` e' desreferencia de ponteiro nulo.
//   E a flag `is_petsc_initialized` e' `static`: inicializar o PETSc por fora NAO
//   a desarma, entao o solver tenta de novo e cai do mesmo jeito.
//
// Logo o nucleo herdado TEM estado global de arranque, e a consequencia para o
// t8Flow e' de projeto: o arranque e' UM, com esta ordem --
//
//   1. `higtree_initialize`  faz o MPI_Init e guarda argc/argv;
//   2. `sc_init` + `t8_init` por processo, com o MPI JA' inicializado (e' o que o
//      adaptador existente ja' faz: "O MPI ja' vem inicializado");
//   3. o PETSc entra preguicoso, pelo proprio solver, agora com argv valido.
//
// NAO chamar `sc_MPI_Init` aqui: o MPI_Init do passo 1 ja' aconteceu, e dois
// MPI_Init sao erro.  Nem `MPI_Finalize`/`PetscFinalize`: o `higtree_initialize`
// registra `_higtree_finalize` no `atexit`, que fecha os dois.
int main(int argc, char **argv)
{
    higtree_initialize(&argc, &argv);
    sc_init(sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ERROR);
    t8_init(SC_LP_ERROR);

    caso_mls_reproduz();
    caso_solver_resolve();
    caso_suporte_da_floresta();

    const int rc = t_end();

    sc_finalize();
    return rc;
}
