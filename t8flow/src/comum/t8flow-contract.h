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
// M -- A MALHA MULTIBLOCOS (decisoes 2 e 7)
// ---------------------------------------------------------------------------

//! M1  Uma arvore por bloco, e as folhas de cada bloco seguem o NIVEL DAQUELE
//!     bloco: 4^nivel em 2D.  So' discrimina com niveis DIFERENTES entre blocos --
//!     com todos iguais, uma malha uniforme global passaria.
//! COBERTA-POR: test-malha2d:uma_arvore_por_bloco
//! COBERTA-POR: test-malha2d:folhas_por_bloco_seguem_o_nivel_do_bloco
//! COBERTA-POR: test-malha2d:total_de_folhas_e_a_soma_dos_blocos

//! M2  A geometria de cada bloco e' a caixa pedida, e a ORIENTACAO do mapeamento
//!     e' a do t8code.  Sao tres afirmacoes porque cada uma pega um erro que as
//!     outras nao veem: centroide dentro da caixa pega geometria que SAI da caixa;
//!     o lado da celula amarra nivel e caixa juntos; e a ordem de percurso (a
//!     folha 1 deslocada em X da folha 0) e' a UNICA que pega espelhamento --
//!     medido, nao suposto.
//! COBERTA-POR: test-malha2d:centroide_cai_na_caixa_do_seu_bloco
//! COBERTA-POR: test-malha2d:celulas_do_bloco_largo_sao_maiores_na_proporcao_certa
//! COBERTA-POR: test-malha2d:a_ordem_de_percurso_varre_x_primeiro

//! M3  As interfaces coincidentes entre blocos sao UNIDAS, e o que a uniao
//!     significa e' afirmado onde importa: face que da' para outro bloco TEM
//!     vizinho, face que da' para fora NAO tem.  Sem a uniao a malha tambem
//!     commita, tambem tem as folhas certas, e toda face de borda de bloco passa
//!     a ser "sem vizinho" -- que e' o criterio de CONTORNO, no meio do dominio.
//!     A CONTAGEM de juncoes nao basta: medido, ela continua certa com o
//!     `t8_cmesh_set_join` removido.
//! COBERTA-POR: test-malha2d:as_interfaces_coincidentes_foram_unidas
//! COBERTA-POR: test-malha2d:face_entre_blocos_tem_vizinho_e_face_externa_nao

//! M4  A malha RECUSA o que nao pode representar, em vez de aproximar: lista
//!     vazia, bloco degenerado, blocos que se sobrepoem em area (o HiGFlow assume
//!     dominios disjuntos, e com sobreposicao a resposta dependeria da ORDEM dos
//!     blocos), e faces que se sobrepoem PARCIALMENTE -- que o join do t8code nao
//!     representa e cujo silencio daria contorno no meio do dominio.
//! COBERTA-POR: test-malha2d:lista_vazia_e_recusada
//! COBERTA-POR: test-malha2d:bloco_degenerado_e_recusado
//! COBERTA-POR: test-malha2d:blocos_sobrepostos_sao_recusados
//! COBERTA-POR: test-malha2d:faces_parcialmente_coincidentes_sao_recusadas

//! M5  O contrato malha/discretizacao vale sobre a malha MULTIBLOCOS: o suporte
//!     sai dos centroides das folhas e os pesos do discretizador, com o alvo no
//!     bloco mais fino e o suporte atravessando resolucoes diferentes.
//! COBERTA-POR: test-malha2d:mls_reproduz_polinomio_sobre_suporte_da_malha
//! COBERTA-POR: test-malha2d:suporte_da_malha_erra_mais_no_grau_acima

// ---------------------------------------------------------------------------
// R -- REFINO LOCAL DENTRO DE UM BLOCO
// ---------------------------------------------------------------------------

//! R1  O refino local sobe ao nivel pedido na regiao pedida, e nao rebaixa nada.
//!     Afirmado com salto de DOIS niveis de uma vez -- um nivel so' nunca
//!     produziria violacao de 2:1, e o R3 passaria sem o equilibrio.
//! COBERTA-POR: test-refino2d:a_malha_comeca_uniforme_em_cada_bloco
//! COBERTA-POR: test-refino2d:o_refino_subiu_ao_nivel_pedido
//! COBERTA-POR: test-refino2d:o_refino_nao_rebaixou_nada
//! COBERTA-POR: test-refino2d:o_refino_aumentou_as_folhas_so_onde_pediu

//! R2  O refino CONSERVA A AREA de cada bloco, folha a folha.  Contagem de
//!     folhas nao substitui isto: celula perdida ou com geometria errada entra na
//!     contagem do mesmo jeito.
//! COBERTA-POR: test-refino2d:o_refino_conserva_a_area_de_cada_bloco

//! R3  A razao entre niveis de celulas vizinhas e' no maximo 2:1, em TODAS as
//!     folhas e TODAS as faces.  E' o que o `t8_forest_set_balance` existe para
//!     garantir; sem ele a malha COMMITA, as contagens fecham, e o salto
//!     atravessa o estencil do discretizador sem aviso.  Medido: com o balance
//!     removido, so' esta clausula reprova (salto de 2 niveis).
//! COBERTA-POR: test-refino2d:a_razao_entre_niveis_vizinhos_e_no_maximo_2_para_1

//! R4  O refino nao VAZA alem de uma celula da caixa pedida, e nao transforma
//!     interface entre blocos em contorno -- a face que da' para o outro bloco
//!     continua com vizinho depois do refino.  O segundo e' o modo de falha mais
//!     caro desta frente.
//! COBERTA-POR: test-refino2d:o_refino_maximo_nao_vaza_alem_de_uma_celula_da_caixa
//! COBERTA-POR: test-refino2d:o_bloco_que_a_caixa_nao_alcanca_ganhou_so_transicao
//! COBERTA-POR: test-refino2d:a_interface_entre_blocos_continua_com_vizinho_depois_do_refino

//! R5  O discretizador reproduz polinomio sobre suporte que ATRAVESSA a transicao
//!     de nivel, e erra mais no grau acima.  E o refino RECUSA o que nao faz
//!     sentido: caixa degenerada, caixa fora de todos os blocos, e nivel alvo que
//!     nao aprofundaria nada -- este ultimo porque aceito calado deixaria a malha
//!     PARECER refinada.
//! COBERTA-POR: test-refino2d:mls_reproduz_polinomio_sobre_suporte_com_transicao_de_nivel
//! COBERTA-POR: test-refino2d:suporte_com_transicao_erra_mais_no_grau_acima
//! COBERTA-POR: test-refino2d:caixa_degenerada_e_recusada
//! COBERTA-POR: test-refino2d:caixa_fora_de_todos_os_blocos_e_recusada
//! COBERTA-POR: test-refino2d:nivel_alvo_que_nao_aprofunda_e_recusado

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
// E uma licao de oraculo desta rodada, que vale guardar: a primeira versao do
// test-malha2d AFIRMAVA que "centroide dentro da caixa" pegava a troca da ordem
// dos vertices do QUAD.  Nao pega -- a troca espelha o mapeamento DENTRO da mesma
// caixa, e todo centroide continua dentro dela.  So' a ordem de percurso acusou.
// O comentario estava errado ate' a sabotagem ser rodada.
//
// Duas licoes de oraculo da rodada do refino, as duas de testes MEUS que estavam
// errados antes do codigo estar:
//
//   * "toda folha de nivel maximo tem o centro dentro da caixa" e' FALSO, e
//     reprovou com 7 folhas fora.  O criterio e' aplicado a' celula ANTES de
//     dividir e os filhos herdam a decisao do pai, entao a regiao refinada e' a
//     caixa DILATADA por ate' uma celula do nivel final.  A afirmacao que vale e'
//     a do R4: o refino nao vaza ALEM de uma celula;
//   * o `t8_forest_commit` CONSOME a floresta passada em `set_adapt`/`set_balance`.
//     O `unref` extra que eu havia escrito deixou a malha VAZIA -- sem erro, sem
//     aviso, contagens em zero.  Quem pegou foi a area (0 contra 1).
//
// Tambem ausentes, previstas no doc/projeto-t8flow.md: numeracao global contigua
// por rank alcancada por busca por ponto (nao pelo percurso que a gerou, que e'
// tautologia); conteudo da franja depois do sync na profundidade de PRODUCAO (5)
// e em np=3, que e' o menor np em que pareamento transposto aparece; a guarda
// contra arvore de UMA celula (que nao tem interior, e a franja declarada e' de
// cinco); a banda fina da decisao 4; posse distribuida; e o lado 3D inteiro.
// ---------------------------------------------------------------------------

#endif
