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
// REGRA DE CRESCIMENTO, deliberada: **uma clausula entra aqui quando tem caso
// que a afirme**, nunca antes.  As previstas e ainda sem codigo vivem no
// `doc/projeto-t8flow.md` -- se entrassem, a suite nasceria vermelha e o vermelho
// deixaria de significar "quebrou".
//
// O contrato vive em `src/comum/` porque vale para as duas dimensoes; os casos
// que o cobrem hoje sao todos do lado 2D, que e' por onde a decisao 7 manda
// comecar.  Quando o 3D existir, cada clausula tera' caso nos dois lados.
//
// NOMENCLATURA: N = nucleo herdado e suas fachadas (decisoes 1 e 7), B =
// construcao e norma (decisao 5), F = floresta do t8code (decisao 2), D =
// dominio unico (decisao 4), P = particao e franja.
#ifndef T8FLOW_CONTRACT_H
#define T8FLOW_CONTRACT_H

// ---------------------------------------------------------------------------
// N -- O NUCLEO HERDADO E AS FACHADAS 2D (decisoes 1 e 7)
// ---------------------------------------------------------------------------

//! N1  O discretizador reproduz EXATAMENTE polinomio de grau <= a ordem do
//!     ajuste, e NAO reproduz grau acima dela.  As duas metades sao uma clausula
//!     so': sem a segunda, "erro pequeno" nao distingue reproducao de
//!     coincidencia.
//! COBERTA-POR: test-mls2d:mls_reproduz_polinomio_do_grau_da_ordem
//! COBERTA-POR: test-mls2d:mls_erra_mais_no_grau_acima_da_ordem

//! N2  Os pesos do discretizador somam 1.  E' a reproducao da constante, e e'
//!     independente do campo de teste -- pega erro de ESCALA nos pesos, que a
//!     reproducao de um polinomio especifico pode absorver.
//! COBERTA-POR: test-mls2d:pesos_somam_um

//! N3  O discretizador RECUSA uso invalido em vez de devolver numero ruim:
//!     suporte abaixo do minimo da ordem, e numero de valores diferente do
//!     numero de pontos.  Numeros ruins em silencio sao o modo de falha que
//!     este projeto mais paga.
//! COBERTA-POR: test-mls2d:suporte_insuficiente_e_recusado
//! COBERTA-POR: test-mls2d:valores_de_tamanho_errado_sao_recusados

//! N4  O solver resolve um sistema de solucao CONHECIDA e reporta convergencia.
//!     O criterio e' o valor da solucao, nao o residuo que o proprio solver
//!     calcula -- um motor que pare na iteracao maxima tambem devolve.
//! COBERTA-POR: test-solver2d:solver_resolve_sistema_de_solucao_conhecida
//! COBERTA-POR: test-solver2d:solucao_inteira_bate_com_a_consulta_individual

//! N5  O nome do motor pedido bate com a tabela do herdado.  As duas listas do
//!     `solver.h` estao DESLOCADAS DE UMA POSICAO e nada confere isso em tempo de
//!     compilacao: errar faz o programa ANUNCIAR um motor e rodar outro.
//! COBERTA-POR: test-solver2d:o_nome_do_motor_pedido_bate_com_a_tabela

//! N6  O solver RECUSA, antes de montar, os usos que abortariam dentro do
//!     herdado: linha definida duas vezes, id fora da faixa local, e linha local
//!     nao definida.  A ultima e' a que importa -- o caminho de pre-alocacao do
//!     PETSc afirma `col_ids[i] && vals[i]` para toda linha, e linha faltando nao
//!     e' sistema singular, e' abort tres camadas abaixo.
//! COBERTA-POR: test-solver2d:linha_repetida_e_recusada
//! COBERTA-POR: test-solver2d:id_fora_da_faixa_local_e_recusado
//! COBERTA-POR: test-solver2d:linha_faltando_e_recusada_antes_de_montar

// ---------------------------------------------------------------------------
// B -- CONSTRUCAO E NORMA (decisao 5)
// ---------------------------------------------------------------------------

//! B1  O nucleo herdado e o t8code convivem na MESMA unidade de traducao,
//!     compilada com o compilador e a norma do t8code (gnu++20), sem
//!     `extern "C"` e sem `-fpermissive`.  Afirmada por um caso que usa os dois
//!     juntos: se nao convivessem, ele nao linkaria.
//! COBERTA-POR: test-floresta2d:mls_reproduz_polinomio_sobre_suporte_da_floresta

// ---------------------------------------------------------------------------
// F -- A FLORESTA DO T8CODE COMO ESTRUTURA (decisao 2)
// ---------------------------------------------------------------------------

//! F1  Refino uniforme de nivel L sobre um hipercubo de uma arvore produz
//!     (2^L)^2 folhas em 2D.  Contagem exata.
//! COBERTA-POR: test-floresta2d:floresta_uniforme_tem_as_folhas_do_nivel

//! F2  O SUPORTE de um ponto e' colhido da floresta -- um centroide por folha --
//!     e alimenta o discretizador, que reproduz polinomio sobre ele.  E' o
//!     contrato malha/discretizacao: a malha entrega suporte, a discretizacao
//!     entrega pesos.
//! COBERTA-POR: test-floresta2d:suporte_colhido_da_floresta_tem_uma_entrada_por_folha
//! COBERTA-POR: test-floresta2d:mls_reproduz_polinomio_sobre_suporte_da_floresta
//! COBERTA-POR: test-floresta2d:suporte_da_floresta_erra_mais_no_grau_acima

// ---------------------------------------------------------------------------
// O QUE ESTE CONTRATO AINDA NAO AFIRMA -- e dizer isto e' parte do contrato.
//
// Medido por sabotagem controlada em 2026-10-11:
//
//   * deslocando TODOS os centroides do suporte meia celula, os casos de F2
//     continuam VERDES.  A reproducao polinomial verifica a DISCRETIZACAO, nao a
//     malha.  Quem verificaria a malha e' "o suporte cai na grade e e' local" --
//     sera' a clausula F3, e ainda nao existe;
//   * baixando a ordem do ajuste de 2 para 1, N1 reprova nas duas metades.  O
//     oraculo discrimina, e isso foi conferido e nao suposto.
//
// E o herdado tem uma lacuna que nenhum teste pode fechar daqui: NAO existe
// `slv_get_engine`, entao qual motor foi INSTALADO nao e' consultavel.  Com
// `Motor::Qualquer` a escolha vem de variavel de ambiente ou heuristica, e a
// unica pista e' a linha que o herdado imprime.  A N5 cobra a tabela, nao a
// instalacao.
//
// Tambem ausentes, previstas no doc/projeto-t8flow.md: malha MULTIBLOCOS (uma
// arvore por bloco, com nivel por bloco); numeracao global contigua por rank
// alcancada por busca por ponto (nao pelo percurso que a gerou, que e'
// tautologia); conteudo da franja depois do sync na profundidade de PRODUCAO (5)
// e em np=3, que e' o menor np em que pareamento transposto aparece; a banda fina
// da decisao 4; posse distribuida; e o lado 3D inteiro.
// ---------------------------------------------------------------------------

#endif
