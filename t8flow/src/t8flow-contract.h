// O CONTRATO DO T8FLOW.
//
// Isto nao e' prosa: o `ci/run_t8flow_tests.py` LE as clausulas deste arquivo e
// exige, para cada uma, um caso de teste que tenha saido PASS na execucao.
// Clausula sem `COBERTA-POR` sai como "SEM TESTE"; clausula cujo caso declarado
// nao apareceu na saida sai como "CASO AUSENTE"; clausula cujo caso reprovou sai
// como "REPROVADA".  As tres contam como FALHA da suite -- nao como omissao.
//
// A razao de ser assim esta' na HiGTree, de onde o padrao vem: a suite ATF de
// 2020 tinha 71 casos que verificavam FORMA (`numelems >= 3`, `rhs != 0`) e
// deixariam passar os nove defeitos de fechamento de 2026-09.  Contrato que o
// driver nao cobra e' documentacao, e documentacao envelhece em silencio.
//
// REGRA DE CRESCIMENTO, e ela e' deliberada: **uma clausula entra aqui quando
// tem caso que a afirme**, nunca antes.  As clausulas previstas e ainda sem
// codigo vivem no `doc/projeto-t8flow.md`, nao aqui -- se entrassem, a suite
// nasceria vermelha e o vermelho deixaria de significar "quebrou".
//
// NOMENCLATURA: N = nucleo herdado (decisao 1), B = construcao e norma
// (decisao 5), F = floresta do t8code (decisao 2), D = dominio unico (decisao 4),
// P = particao e franja, S = suite.
#ifndef T8FLOW_CONTRACT_H
#define T8FLOW_CONTRACT_H

// ---------------------------------------------------------------------------
// N -- O NUCLEO HERDADO DA HIGTREE (decisao 1)
// ---------------------------------------------------------------------------

//! N1  A aproximacao por minimos quadrados moveis reproduz EXATAMENTE polinomio
//!     de grau <= a ordem do ajuste, e NAO reproduz grau acima dela.  As duas
//!     metades sao uma clausula so': sem a segunda, "erro pequeno" nao distingue
//!     reproducao de coincidencia.
//! COBERTA-POR: test-nucleo-herdado:mls_reproduz_polinomio_do_grau_da_ordem
//! COBERTA-POR: test-nucleo-herdado:mls_erra_mais_no_grau_acima_da_ordem

//! N2  O solver linear herdado resolve um sistema de solucao CONHECIDA, e
//!     reporta convergencia.  O criterio e' o valor da solucao, nao o residuo
//!     que o proprio solver calcula.
//! COBERTA-POR: test-nucleo-herdado:solver_resolve_sistema_de_solucao_conhecida

// ---------------------------------------------------------------------------
// B -- CONSTRUCAO E NORMA (decisao 5)
// ---------------------------------------------------------------------------

//! B1  O nucleo herdado e o t8code convivem na MESMA unidade de traducao,
//!     compilada com o compilador e a norma do t8code (gnu++20), sem
//!     `extern "C"` e sem `-fpermissive`.  Afirmada por um caso que usa os dois
//!     juntos: se nao convivessem, ele nao linkaria.
//! COBERTA-POR: test-nucleo-herdado:mls_reproduz_polinomio_sobre_suporte_da_floresta

// ---------------------------------------------------------------------------
// F -- A FLORESTA DO T8CODE COMO ESTRUTURA (decisao 2)
// ---------------------------------------------------------------------------

//! F1  Refino uniforme de nivel L sobre um hipercubo de uma arvore produz
//!     (2^L)^DIM folhas.  Contagem exata, nos dois DIM.
//! COBERTA-POR: test-nucleo-herdado:floresta_uniforme_tem_as_folhas_do_nivel

//! F2  O SUPORTE de um ponto e' colhido da floresta -- um centroide por folha --
//!     e alimenta o ajuste herdado, que reproduz polinomio sobre ele.  E' o
//!     contrato malha/discretizacao: a malha entrega suporte, a discretizacao
//!     entrega pesos.
//! COBERTA-POR: test-nucleo-herdado:suporte_colhido_da_floresta_tem_uma_entrada_por_folha
//! COBERTA-POR: test-nucleo-herdado:mls_reproduz_polinomio_sobre_suporte_da_floresta
//! COBERTA-POR: test-nucleo-herdado:suporte_da_floresta_erra_mais_no_grau_acima

// ---------------------------------------------------------------------------
// O QUE ESTE CONTRATO AINDA NAO AFIRMA -- e dizer isto e' parte do contrato.
//
// Medido em 2026-10-11, por sabotagem controlada:
//
//   * deslocando TODOS os centroides do suporte meia celula, os tres casos de F2
//     continuam VERDES.  A reproducao polinomial verifica a DISCRETIZACAO, nao a
//     malha.  Quem verificaria a malha e' "o suporte cai na grade e e' local", que
//     ainda nao existe -- sera' a clausula F3.
//   * baixando a ordem do ajuste de 2 para 1, N1 reprova nas duas metades.  O
//     oraculo discrimina, e isso foi conferido e nao suposto.
//
// Tambem ausentes, e previstas no doc/projeto-t8flow.md: numeracao global
// contigua por rank alcancada por busca por ponto (nao pelo percurso que a
// gerou, que e' tautologia); conteudo da franja depois do sync na profundidade de
// PRODUCAO (5) e em np=3, que e' o menor np em que pareamento transposto
// aparece; a banda fina da decisao 4; e posse distribuida.
// ---------------------------------------------------------------------------

#endif
