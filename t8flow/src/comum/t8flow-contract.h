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

//! N7  PESOS DE DERIVADA: `sum w_k u_k` aproxima `du/dx_dir` no alvo, EXATO para
//!     polinomio de grau <= a ordem, e os pesos SOMAM ZERO (a derivada de uma
//!     constante) -- a conferencia barata que pega escala e sinal sem depender de
//!     campo nenhum.  Os de valor somam um; estes somam zero.
//!
//!     E A IMPLEMENTACAO PROPRIA CONCORDA COM O NUCLEO HERDADO no caso de valor --
//!     o herdado resolve por QR, a propria por equacoes normais, e a diferenca
//!     medida e' 1,2e-16.  Sem essa concordancia nao haveria oraculo para a
//!     derivada, que so' existe na propria: `wls_calc` devolve pesos de VALOR e
//!     deixa a matriz fatorada interna.
//!
//!     E' a peca que a reconstrucao de fluxo de alta ordem em malha adaptada
//!     precisa, e que estava declarada como limite desde o inicio.
//! COBERTA-POR: test-mls2d:os_pesos_proprios_de_valor_concordam_com_o_nucleo_herdado
//! COBERTA-POR: test-mls2d:a_derivada_e_exata_para_polinomio_do_grau_da_ordem
//! COBERTA-POR: test-mls2d:a_derivada_erra_no_grau_acima_da_ordem
//! COBERTA-POR: test-mls2d:os_pesos_de_derivada_somam_zero
//! COBERTA-POR: test-mls2d:derivada_malformada_e_recusada

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

//! F3  O SUPORTE CAI NA GRADE E E' LOCAL.  "Na grade": cada ponto do suporte e' o
//!     centroide de uma celula que EXISTE -- conferido pedindo a' malha o id
//!     daquele ponto e comparando com o id que o suporte declara -- ou o centro de
//!     uma faceta de CONTORNO, e esse esta' sobre a borda do dominio.  "Local": o
//!     diametro e' de poucas celulas e o numero de pontos e' menor que a malha.
//!
//!     E' a clausula que a REPRODUCAO POLINOMIAL NAO DA'.  Medido: deslocando os
//!     centroides do suporte meia celula, F2 continua verde -- a reproducao
//!     verifica a discretizacao, nao a malha.  Com a F3, o mesmo deslocamento
//!     reprova com "29 de 50 pontos do suporte fora da grade".
//! COBERTA-POR: test-estencil2d:F3_todo_ponto_do_suporte_esta_na_grade
//! COBERTA-POR: test-estencil2d:F3_o_suporte_e_local_e_nao_a_malha_inteira

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

//! R6  O REFINO ALEATORIO e' REPRODUTIVEL pela semente, e a semente FAZ
//!     DIFERENCA.  As duas metades: reprodutivel sozinho seria propriedade de
//!     qualquer codigo sem aleatoriedade; e a reproducao e' afirmada pelos
//!     CENTROIDES, nao pela contagem -- duas malhas diferentes com o mesmo numero
//!     de folhas sao perfeitamente possiveis.  Medido: semente 12345 da' 121
//!     folhas, semente 999 da' 115.
//! COBERTA-POR: test-refino2d:o_refino_aleatorio_e_reprodutivel_pela_semente
//! COBERTA-POR: test-refino2d:sementes_diferentes_dao_malhas_diferentes
//! COBERTA-POR: test-refino2d:a_fracao_controla_quanto_se_refina

//! R7  A malha aleatoria respeita 2:1, CONSERVA A AREA, e tem TRANSICAO DE
//!     VERDADE.  A ultima parte e' o que impede um teste de esquema na face de
//!     transicao de rodar numa malha sem transicao nenhuma e passar sem testar
//!     nada.  Medido com duas passadas: 310 folhas, niveis 3 a 5, 236 folhas em
//!     transicao (76%).
//! COBERTA-POR: test-refino2d:o_refino_aleatorio_respeita_2_para_1_e_tem_transicao_de_verdade

//! R8  Refino aleatorio malformado e' recusado: fracao fora de (0,1], e
//!     `nivel_max` que nao aprofunda -- refino que nao refina nada deixaria a
//!     malha PARECER adaptada.
//! COBERTA-POR: test-refino2d:refino_aleatorio_malformado_e_recusado

// ---------------------------------------------------------------------------
// G -- A NUMERACAO GLOBAL (o que o solver pede da malha)
// ---------------------------------------------------------------------------

//! G1  A numeracao local e' CONTIGUA e comeca em `primeiro_id()`; o total global
//!     fecha com a soma dos locais.
//! COBERTA-POR: test-numeracao2d:a_numeracao_local_e_contigua_e_comeca_no_primeiro_id

//! G2  Os ids cobrem [primeiro, primeiro + n_local) exatamente uma vez cada --
//!     sem furo e sem repeticao.  E' o defeito mais comum de numeracao, e nenhuma
//!     contagem de folhas o ve.
//! COBERTA-POR: test-numeracao2d:os_ids_cobrem_o_intervalo_local_sem_furo_nem_repeticao

//! G3  A ordem local E' a curva de preenchimento: num bloco uniforme de nivel L,
//!     a folha de indice k esta' na posicao obtida DESINTERCALANDO os bits de k.
//!     A conta sai da definicao da curva, feita no teste sem olhar a malha, e e'
//!     comparada com a geometria da folha -- via independente, nao contagem.
//! COBERTA-POR: test-numeracao2d:a_ordem_local_e_a_curva_de_morton

//! G4  A busca POR PONTO devolve o id da folha que contem o ponto, devolve -1
//!     fora do dominio, e DESEMPATA explicitamente: ponto sobre a face comum
//!     pertence a' celula de coordenada MAIOR (intervalo semiaberto), com o limite
//!     superior do dominio como excecao.  O desempate nao pode depender da ordem
//!     de visita -- na HiGTree a clausula analoga passava por SORTE, porque a
//!     curva visita a celula de menor coordenada primeiro.
//! COBERTA-POR: test-numeracao2d:a_busca_por_ponto_devolve_o_id_da_folha_que_contem_o_ponto
//! COBERTA-POR: test-numeracao2d:ponto_fora_de_todos_os_blocos_nao_tem_id
//! COBERTA-POR: test-numeracao2d:o_desempate_na_face_comum_e_a_celula_de_coordenada_maior
//! COBERTA-POR: test-numeracao2d:o_limite_superior_do_dominio_pertence_a_ultima_celula

//! G5  Quem ESCREVE por id de percurso le' o mesmo valor por id de GEOMETRIA --
//!     um sistema identidade e' montado com os ids do percurso e lido com os ids
//!     que `id_do_ponto` devolve.  E' o defeito do sdmult/sdp do HiGFlow, que
//!     sobreviveu anos porque em malha uniforme os dois indices coincidem celula a
//!     celula.  Vale tambem DEPOIS do refino local, que e' quando a numeracao muda
//!     inteira e os niveis deixam de ser uniformes.
//! COBERTA-POR: test-numeracao2d:o_valor_escrito_por_id_de_percurso_e_lido_por_id_de_geometria
//! COBERTA-POR: test-numeracao2d:a_consistencia_sobrevive_ao_refino_local

// ---------------------------------------------------------------------------
// A -- O ARRANJO DESLOCADO: onde cada propriedade mora
// ---------------------------------------------------------------------------

//! A1  Propriedade ESCALAR vive no CENTRO da celula: um valor por folha,
//!     indexado pelo id GLOBAL -- o mesmo do solver, nao um indice proprio.
//!     Escrito por id, lido por geometria, e id fora da faixa local e' recusado.
//! COBERTA-POR: test-campos2d:campo_escalar_tem_um_valor_por_celula
//! COBERTA-POR: test-campos2d:escalar_escrito_por_id_e_lido_por_geometria
//! COBERTA-POR: test-campos2d:escalar_recusa_id_fora_da_faixa_local

//! A2  Propriedade TENSORIAL vive no CENTRO, com as componentes INDEPENDENTES --
//!     escrever numa nao altera as outras, e a componente vista como campo
//!     escalar e' a MESMA memoria.  Layout sobreposto passa em qualquer teste de
//!     valor unico; so' a independencia o pega.
//! COBERTA-POR: test-campos2d:componentes_do_tensor_sao_independentes
//! COBERTA-POR: test-campos2d:tensor_recusa_componente_inexistente

//! A3  As FACETAS tem posse UNICA: a faceta da direcao d e' a face "menos" de
//!     exatamente uma folha, e as faces "mais" de contorno ganham faceta propria.
//!     A contagem segue a conta classica do MAC, e a posse unica sobrevive ao
//!     refino 2:1 -- onde a face "mais" de uma celula grossa NAO vira faceta,
//!     porque as duas finas possuem as suas "menos" e o lado fino fica com dois
//!     valores.  Contagem certa nao garante posse unica: duas trocas que se
//!     compensam dariam o mesmo total.
//! COBERTA-POR: test-campos2d:a_contagem_de_facetas_segue_a_conta_do_mac
//! COBERTA-POR: test-campos2d:cada_faceta_tem_exatamente_um_dono
//! COBERTA-POR: test-campos2d:no_refino_2_para_1_o_lado_fino_tem_duas_facetas

//! A4  A GEOMETRIA da faceta: o centro esta' a meia celula do centro da celula,
//!     na direcao da normal e so' nela; a MEDIDA e' o lado TRANSVERSAL.  O caso
//!     tem um bloco de celulas NAO-QUADRADAS, sem o que usar o lado errado como
//!     medida seria invisivel.
//! COBERTA-POR: test-campos2d:geometria_e_medida_da_faceta

//! A5  Propriedade VETORIAL vive nas FACETAS, uma componente por direcao, cada
//!     uma nas facetas daquela direcao -- e as direcoes sao independentes.  E o
//!     oraculo que amarra posse, geometria, medida e sinal de uma vez: um campo
//!     CONSTANTE tem divergencia ZERO em toda celula interna.
//! COBERTA-POR: test-campos2d:campo_vetorial_tem_um_valor_por_faceta_de_cada_direcao
//! COBERTA-POR: test-campos2d:as_duas_direcoes_do_vetorial_sao_independentes
//! COBERTA-POR: test-campos2d:divergencia_de_campo_constante_e_zero

// ---------------------------------------------------------------------------
// L -- O LOTE DE ESTENCEIS, com CPU e GPU
// ---------------------------------------------------------------------------

//! L1  O lote de CPU devolve um estencil por consulta, com os ids PRESERVADOS na
//!     ordem do suporte, e reproduz o grau da ordem pedida -- errando no grau
//!     acima.  A CPU e' implementacao de PRIMEIRA CLASSE: funciona sem OpenCL
//!     instalado, sem OpenCL compilado e sem placa.
//!
//!     E O PADRAO DE BACK-END E' DE QUEM COMPILA: `make GPU=sim` faz a GPU ser o
//!     padrao, `GPU=nao` tira o codigo de placa do binario, e sem a variavel a
//!     placa fica disponivel com padrao CPU.  Tres estados, porque "detectei" e
//!     "o usuario pediu" sao diferentes.  O caso afirma a COERENCIA entre a flag
//!     e o comportamento -- o unico lado afirmavel de dentro de um build --, e a
//!     saida declara qual build foi testado.
//! COBERTA-POR: test-lote2d:o_back_end_padrao_segue_a_flag_de_compilacao
//! COBERTA-POR: test-lote2d:o_lote_de_cpu_devolve_um_estencil_por_consulta_com_ids_preservados
//! COBERTA-POR: test-lote2d:cpu_reproduz_o_grau_da_ordem
//! COBERTA-POR: test-lote2d:cpu_erra_mais_no_grau_acima_da_ordem

//! L2  Os pesos somam 1 -- reproducao da constante, independente do campo de
//!     teste.  Pega erro de ESCALA que a reproducao de um polinomio especifico
//!     absorve.
//! COBERTA-POR: test-lote2d:os_pesos_somam_um

//! L3  O lote RECUSA uso invalido ANTES de escolher o back-end, identificando a
//!     consulta: suporte abaixo do minimo da ordem, e numero de ids diferente do
//!     numero de pontos.  Validar depois de escolher faria "funciona na CPU e
//!     quebra na placa" ser possivel.
//! COBERTA-POR: test-lote2d:suporte_insuficiente_e_recusado_identificando_a_consulta
//! COBERTA-POR: test-lote2d:ids_em_numero_diferente_do_suporte_e_recusado

//! L4  Sem GPU o lote CONTINUA FUNCIONANDO e DIZ por que nao ha' placa; e pedir
//!     GPU explicitamente sem GPU LANCA, em vez de cair na CPU em silencio -- quem
//!     pediu placa quer saber que nao ha'.
//! COBERTA-POR: test-lote2d:sem_gpu_o_lote_continua_funcionando_e_diz_por_que

//! L5  Com GPU, os PESOS batem com os da CPU (nao o valor interpolado: valor pode
//!     bater por cancelamento com pesos errados), e a GPU passa no MESMO oraculo
//!     de reproducao -- bater com a CPU nao basta se as duas estiverem erradas
//!     juntas.  Lote heterogeneo cai na CPU com o resultado CERTO, e ordem que o
//!     kernel nao cobre se declara indisponivel em vez de calcular errado.
//!
//!     **CONDICIONAL, e a saida declara:** sem placa estes casos SALTAM e passam.
//!     O teste imprime `gpu-exercitada sim|nao` e o driver reporta, para
//!     "contrato 100%" nunca significar "a GPU nao foi testada".
//! COBERTA-POR: test-lote2d:gpu_bate_com_cpu_nos_pesos
//! COBERTA-POR: test-lote2d:lote_heterogeneo_cai_na_cpu_com_o_resultado_certo
//! COBERTA-POR: test-lote2d:ordem_que_o_kernel_nao_cobre_nao_vai_para_a_gpu_em_silencio

// ---------------------------------------------------------------------------
// E -- O ESTENCIL COM A ORDEM DESEJADA (como na HiGTree)
// ---------------------------------------------------------------------------

//! E1  O estencil reproduz EXATAMENTE o grau da ordem pedida no interior, e erra
//!     no grau acima.  A ordem e' a pedida, nao a que o suporte permitiu: suporte
//!     insuficiente faz a chamada LANCAR em vez de devolver um ajuste de ordem
//!     menor, que daria numeros plausiveis e errados.
//! COBERTA-POR: test-estencil2d:o_estencil_reproduz_o_grau_da_ordem_no_interior
//! COBERTA-POR: test-estencil2d:o_estencil_erra_mais_no_grau_acima_da_ordem
//! COBERTA-POR: test-estencil2d:sem_condicao_de_contorno_o_suporte_perto_da_parede_e_recusado
//! COBERTA-POR: test-estencil2d:ponto_fora_do_dominio_e_recusado

//! E2  O FECHAMENTO DE CONTORNO funciona: junto a' parede, parte do suporte vem
//!     das facetas de contorno, o peso daquele ponto multiplica o valor de
//!     Dirichlet e vai para o `rhs` -- e o ponto NAO aparece em `ids`.  A
//!     reproducao vale ate' a parede, que e' o unico oraculo que exercita o `rhs`.
//!
//!     O campo de contorno VARIA nas duas direcoes de proposito: com valor
//!     constante -- parede de nao-deslizamento, que e' quase todo exemplo --
//!     fechar contra a parede ERRADA da' o mesmo numero, e foi assim que nove
//!     defeitos de fechamento sobreviveram cinco anos na HiGTree.  Medido: fazendo
//!     o ponto de contorno entrar como CELULA, este caso reprova com "nenhum
//!     estencil junto a' parede teve rhs".
//! COBERTA-POR: test-estencil2d:o_estencil_reproduz_junto_a_parede_pelo_fechamento

//! E3  O ANEL E' FECHADO: o crescimento do suporte nunca para no meio de um anel
//!     de vizinhanca, e o conjunto de ids nao depende da ordem de visita.  Parar
//!     no meio faz o suporte -- e os pesos -- depender do percurso, que e' como a
//!     clausula de desempate da HiGTree passava por SORTE.  Medido: decidindo a
//!     parada dentro do laco de faces, o suporte cai para 3 pontos onde o anel
//!     fechado da' 5.
//! COBERTA-POR: test-estencil2d:o_anel_e_fechado_e_o_suporte_nao_depende_da_ordem_de_percurso

//! E4  O LOTE de estenceis da' o mesmo que um-por-um, com o `rhs` preservado nos
//!     dois caminhos.  A selecao do suporte fica na CPU (e' navegacao de malha,
//!     nao aritmetica densa); o que vai para a placa sao os pesos.
//! COBERTA-POR: test-estencil2d:o_lote_de_estenceis_da_o_mesmo_que_um_por_um

// ---------------------------------------------------------------------------
// C -- AS CONDICOES DE FRONTEIRA: Dirichlet, Neumann e Robin
// ---------------------------------------------------------------------------

//! C1  Os coeficientes afins `(A, B)` de `u_fantasma = A*u_celula + B` batem com a
//!     algebra, nos TRES tipos, e a conta e' UMA so' parametrizada por (alfa,
//!     beta) -- tres contas separadas divergiriam, e divergir aqui significa Robin
//!     certo e Neumann errado, ou pior, os dois plausiveis.  Neumann e' (1, h*g);
//!     Dirichlet, no limite beta=0, e' (-1, 2g) -- a reflexao antissimetrica.
//! COBERTA-POR: test-contorno2d:os_coeficientes_afins_batem_com_a_algebra

//! C2  Condicao DEGENERADA e' recusada, nao aproximada: alfa=beta=0 nao e'
//!     condicao, `alfa*h/2 + beta = 0` nao determina o fantasma, e h=0 nao e'
//!     distancia.  Devolver numero ali seria inventar dado.
//! COBERTA-POR: test-contorno2d:condicao_degenerada_e_recusada

//! C3  As TRES condicoes reproduzem um campo analitico ATE' A PAREDE, alimentadas
//!     com o dado EXATO que aquele campo satisfaz -- `u` para Dirichlet, `du/dn`
//!     para Neumann, `alfa*u + beta*du/dn` para Robin.  E Robin com beta=0 fecha
//!     com o caso conhecido de Dirichlet, o que amarra a formula.
//! COBERTA-POR: test-contorno2d:dirichlet_reproduz_o_campo_ate_a_parede
//! COBERTA-POR: test-contorno2d:neumann_reproduz_o_campo_ate_a_parede_pelo_fantasma
//! COBERTA-POR: test-contorno2d:robin_reproduz_o_campo_ate_a_parede
//! COBERTA-POR: test-contorno2d:robin_com_beta_zero_concorda_com_dirichlet_na_reproducao

//! C4  O MECANISMO e' observavel e DADO ERRADO PIORA MUITO.  Dirichlet amostra a
//!     PAREDE e nao poe fantasma; Neumann poe FANTASMA do lado de FORA e nao
//!     amostra a parede -- se os dois fizessem o mesmo, a conta de (A,B) nunca
//!     seria exercitada.  E alimentar a condicao com o dado de OUTRO campo (ou
//!     com outros alfa/beta) tem de piorar o erro em ordens de grandeza E em
//!     VALOR ABSOLUTO: medido 5,8e-16 contra 1,9e-01.
//! COBERTA-POR: test-contorno2d:o_suporte_de_neumann_usa_fantasma_e_nao_amostra_a_parede
//! COBERTA-POR: test-contorno2d:neumann_com_a_derivada_de_outro_campo_piora
//! COBERTA-POR: test-contorno2d:robin_com_alfa_e_beta_errados_piora
//! COBERTA-POR: test-contorno2d:sem_condicao_instalada_o_suporte_curto_e_recusado

//! C5  NO PROPRIO CENTRO DA CELULA O ESTENCIL E' A IDENTIDADE, e isto esta'
//!     AFIRMADO em vez de evitado em silencio: o peso do MLS e' `1/dist^2`, logo o
//!     ponto coincidente recebe `1/EPSMACH` e domina -- peso 1,0 na propria celula,
//!     ~0 nas outras, `rhs` ~0.  E' propriedade do ajuste movel, nao defeito; o
//!     caso existe para ninguem voltar a medir fechamento ali e achar que mediu.
//! COBERTA-POR: test-contorno2d:no_proprio_centro_da_celula_o_estencil_e_a_identidade

//! C6  CONDICAO POR REGIAO: paredes diferentes com condicoes diferentes, cada
//!     regiao com nome e predicado de POSICAO.  Afirmada no ensaio que o Poisson
//!     vai precisar -- quatro paredes, TRES tipos (Dirichlet no sul, Neumann no
//!     norte, Robin com coeficientes distintos no oeste e no leste), todas
//!     alimentadas com o dado exato do mesmo campo, e a reproducao valendo nas
//!     quatro.  O caso tambem afirma QUAL regiao responde em cada parede: sem
//!     isso, "a condicao certa foi aplicada" nao e' verificavel.
//!
//!     O predicado e' GEOMETRICO e nao por id de faceta: id depende da numeracao,
//!     que muda com refino e particao -- contorno declarado por id teria de ser
//!     refeito a cada adaptacao.  (O HiGFlow declara por `sim_boundary` com id, e
//!     paga esse preco.)
//! COBERTA-POR: test-contorno2d:quatro_paredes_com_tres_tipos_de_condicao_reproduzem_o_campo
//! COBERTA-POR: test-contorno2d:trocar_a_condicao_de_uma_so_parede_piora

//! C7  PAREDE SEM REGIAO e' RECUSADA, e PAREDE EM DUAS REGIOES tambem.
//!
//!     A primeira porque contorno sem condicao faz a conta rodar com a parede
//!     ausente -- resultado plausivel e errado, e dos defeitos mais caros em
//!     simulacao.  Quem quiser cobrir o resto instala um PADRAO, explicitamente.
//!
//!     A segunda porque a alternativa -- "a primeira que casa vence" -- faria a
//!     condicao depender da ORDEM DE REGISTRO.  E' o oposto do que o
//!     `cell_find_in_center` do HiGFlow faz, de proposito: lá o `break` no
//!     primeiro bloco e' justamente o que a clausula C8 proibe para localizacao.
//! COBERTA-POR: test-contorno2d:parede_sem_regiao_e_recusada
//! COBERTA-POR: test-contorno2d:parede_em_duas_regioes_e_recusada
//! COBERTA-POR: test-contorno2d:o_padrao_cobre_o_que_nenhuma_regiao_cobre

//! C8  REGIAO MALFORMADA e' recusada no registro: nome repetido (o diagnostico
//!     ficaria ambiguo justamente quando importa), predicado ausente, dado `g`
//!     ausente.
//! COBERTA-POR: test-contorno2d:regiao_malformada_e_recusada

// ---------------------------------------------------------------------------
// P -- A EQUACAO DE POISSON, por solucoes manufaturadas
// ---------------------------------------------------------------------------

//! P1  Com solucao de GRAU 2 e fechamento de TRES pontos, o erro e' de MAQUINA nas
//!     DUAS normas -- do maximo e L2 discreta ponderada por area.  Vale com todas
//!     as paredes de Dirichlet, com duas de Neumann, e com TRES de Neumann (uma so'
//!     de Dirichlet, o minimo para o problema nao ser singular).
//!
//!     E A EXATIDAO NAO DEPENDE DO REFINO: se o erro caisse como h^2 seria
//!     convergencia, nao exatidao, e o caso estaria passando por a malha ser fina.
//!     O caso exige que a razao entre os niveis 3 e 5 fique entre 0,25 e 4.
//! COBERTA-POR: test-poisson2d:grau2_com_dirichlet_da_erro_de_maquina_nas_duas_normas
//! COBERTA-POR: test-poisson2d:a_exatidao_do_grau2_nao_depende_do_refino
//! COBERTA-POR: test-poisson2d:grau2_com_neumann_em_duas_paredes_tambem_da_erro_de_maquina
//! COBERTA-POR: test-poisson2d:grau2_com_uma_parede_dirichlet_e_tres_neumann_da_erro_de_maquina

//! P2  O FECHAMENTO DE DOIS PONTOS perde a exatidao, e converge em ~2.  As duas
//!     metades importam: se ele tambem desse erro de maquina, a P1 nao
//!     distinguiria esquema exato de sistema montado como IDENTIDADE -- que
//!     tambem daria 1e-13.  Medido: 1,16e-03 contra 3,91e-13, dez ordens de
//!     grandeza, com ordem 1,82 -> 1,93.
//! COBERTA-POR: test-poisson2d:o_fechamento_de_dois_pontos_perde_a_exatidao
//! COBERTA-POR: test-poisson2d:o_fechamento_de_dois_pontos_converge_em_segunda_ordem

//! P3  ORDEM DE CONVERGENCIA ~2 onde HA' erro para medir: solucao trigonometrica,
//!     com Dirichlet em tudo e com tres Neumann.  E a ordem so' e' calculada se o
//!     erro estiver acima do ruido do solver -- `log2` de ruido sobre ruido e' um
//!     numero aleatorio que pareceria uma taxa.
//!
//!     E O GRAU 3 COM TRES NEUMANN cai para ordem ~2, e isso NAO e' defeito: o
//!     dado de Neumann e' o fluxo exato NO CENTRO da faceta, e o balanco precisa
//!     do fluxo MEDIO.  A quadratura do ponto medio e' exata para integrando
//!     LINEAR -- grau 2 tem `du/dn` linear na faceta, grau 3 tem quadratico.  E' a
//!     razao de a exatidao do grau 2 nao ser sorte: a quadratura e' exata nos DOIS
//!     lugares, no fluxo e na fonte.
//! COBERTA-POR: test-poisson2d:a_solucao_trigonometrica_converge_em_segunda_ordem
//! COBERTA-POR: test-poisson2d:com_tres_neumann_o_grau3_converge_em_segunda_ordem
//! COBERTA-POR: test-poisson2d:com_tres_neumann_a_trigonometrica_tambem_converge_em_segunda_ordem

//! P4  A MONTAGEM RECUSA o que nao sabe, em vez de produzir resultado plausivel:
//!     Robin (mistura valor e derivada na MESMA faceta), fonte nao instalada,
//!     contorno vazio, e nome de solucao inexistente -- este com a LISTA dos
//!     validos na mensagem, porque errar o nome e receber a solucao padrao em
//!     silencio seria medir uma coisa e reportar outra.
//!
//!     A transicao 2:1 SAIU desta lista: ela era recusada e hoje e' montada, e a
//!     familia T cobre o que isso exige.  O driver pegou a mudanca sozinho -- a
//!     clausula ficou apontando para um caso removido e saiu como "CASO AUSENTE",
//!     que e' exatamente o trabalho dele.
//! COBERTA-POR: test-poisson2d:fonte_nao_instalada_e_recusada
//! COBERTA-POR: test-poisson2d:contorno_vazio_e_recusado
//! COBERTA-POR: test-poisson2d:nome_de_solucao_inexistente_e_recusado_com_a_lista

// ---------------------------------------------------------------------------
// Q -- O ESQUEMA COMPACTO DE QUARTA ORDEM (malha uniforme)
// ---------------------------------------------------------------------------

//! Q1  O compacto e' EXATO ATE' GRAU 3, e o LIMITE E' O FECHAMENTO -- nao o
//!     interior.  O estencil de nove pontos tem erro proporcional a derivadas de
//!     SEXTA ordem, logo seria exato ate' grau 5; o fechamento de parede usa
//!     quatro pontos (parede e tres celulas) e e' exato ate' grau 3.  Medido, e o
//!     corte e' nitido:
//!
//!         grau1 1,44e-14   grau2 1,12e-13   grau3 3,58e-13   grau4 3,99e-05
//!
//!     A clausula afirma o LIMITE e nao so' a exatidao: sem o grau 4 reprovando,
//!     nao se saberia ONDE o esquema deixa de ser exato.
//! COBERTA-POR: test-poisson2d:o_compacto_e_exato_ate_grau_3_e_o_limite_e_o_fechamento

//! Q2  O compacto converge em QUARTA ordem (medido ~3,7 a 4,0 nas duas normas), e
//!     e' ordens de grandeza mais preciso que o de segunda no MESMO nivel -- 70x
//!     na trigonometrica, nivel 5.  As duas metades: "ordem 4" sozinha poderia
//!     conviver com um erro absoluto pior.
//! COBERTA-POR: test-poisson2d:o_compacto_converge_em_QUARTA_ordem
//! COBERTA-POR: test-poisson2d:o_compacto_e_muito_mais_preciso_que_o_de_segunda_ordem

//! Q3  MEHRSTELLEN E LELE SAO O MESMO ESQUEMA, verificado em vez de suposto.  O
//!     primeiro usa `f + (h^2/12) lap(f)` e exige `lap(f)` ANALITICO; o segundo e'
//!     implicito na segunda derivada e usa `f` nos nove pontos.  A algebra:
//!     `A_x f / (6/5)` expande em `f + (h^2/12) f_xx + O(h^4)` porque
//!     `(1/10)/(6/5) = 1/12`.
//!
//!     Medido: em polinomio de grau <= 4 dao valores IDENTICOS (ali o termo de
//!     O(h^4) do lado direito envolve derivadas quartas de `f`, nulas); na
//!     trigonometrica diferem 1,8%, que e' a ordem do proprio truncamento.
//!
//!     Ter as duas variantes e' o que permite esta verificacao -- e tem valor
//!     pratico: `lap(f)` existe na manufatura e nao existe em problema real.
//! COBERTA-POR: test-poisson2d:mehrstellen_e_lele_concordam_na_ordem_do_truncamento

//! Q4  O compacto RECUSA o que nao sabe: malha ADAPTADA (o estencil usa vizinhos
//!     DIAGONAIS, que na transicao 2:1 nao estao definidos sem reconstrucao
//!     propria), mais de um bloco, NEUMANN (o fechamento e' de Dirichlet; usa-lo
//!     para Neumann daria quarta ordem no interior e segunda na parede, sem
//!     aviso), e Mehrstellen SEM `lap(f)` -- este ultimo porque rodar com zero
//!     daria SEGUNDA ordem e pareceria quarta falhando, que e' o pior resultado
//!     possivel: o numero sai plausivel.
//! COBERTA-POR: test-poisson2d:o_compacto_recusa_o_que_nao_sabe

// ---------------------------------------------------------------------------
// T -- A TRANSICAO 2:1 NO OPERADOR (quadtree adaptada)
// ---------------------------------------------------------------------------

//! T1  Na malha ADAPTADA o operador e' EXATO PARA CAMPO LINEAR -- residuo
//!     `max |A*u_exata - b| / area` abaixo de 1e-9.  Campo linear e' o PISO: se o
//!     operador nao reproduz linear, ele nao e' consistente e nenhuma ordem de
//!     convergencia medida depois significa nada.
//!
//!     E O ORACULO E' O RESIDUO DO OPERADOR, nao o erro da solucao.  A razao: o
//!     erro da solucao mistura consistencia, condicionamento e tolerancia do
//!     solver, e quando ele e' grande nao se sabe qual dos tres.  Foi exatamente o
//!     que faltou quando a primeira versao do esquema deu ordem NEGATIVA -- havia
//!     um defeito de O(1) na montagem, e medir so' a solucao nao o isolava.
//! COBERTA-POR: test-poisson2d:na_malha_adaptada_o_operador_e_exato_para_campo_linear
//! COBERTA-POR: test-poisson2d:na_malha_adaptada_com_neumann_o_linear_continua_exato

//! T2  Na malha adaptada o GRAU 2 DEIXA DE SER EXATO e converge em ~2, e dizer
//!     isso e' parte da clausula: a exatidao era propriedade da malha UNIFORME.
//!     Na transicao, a distancia entre os centros (3h/4) nao e' a distancia do
//!     centro a' face (h/2), e o termo de curvatura nao cancela -- o que a analise
//!     da quadratura do ponto medio ja' previa.  Nao e' defeito; afirmar que deixa
//!     de ser exato e' o que impede alguem de "consertar" o que esta' certo.
//! COBERTA-POR: test-poisson2d:na_malha_adaptada_o_grau2_deixa_de_ser_exato_e_converge_em_segunda_ordem
//! COBERTA-POR: test-poisson2d:na_malha_adaptada_a_trigonometrica_converge_em_segunda_ordem

//! T3  O RESIDUO DO OPERADOR DISTINGUE malha uniforme de adaptada no grau 2 --
//!     maquina contra O(1).  E' o par que da' sentido a T1 e T2: se os dois dessem
//!     o mesmo numero, o esquema de transicao nao estaria sendo exercitado em
//!     nenhum deles.
//! COBERTA-POR: test-poisson2d:o_residuo_do_operador_distingue_malha_uniforme_de_adaptada_no_grau2

// ---------------------------------------------------------------------------
// O QUE ESTE CONTRATO AINDA NAO AFIRMA -- e dizer isto e' parte do contrato.
//
// Medido por sabotagem controlada em 2026-10-11:
//
//   * deslocando TODOS os centroides do suporte meia celula, os casos de F2
//     continuam VERDES.  A reproducao polinomial verifica a DISCRETIZACAO, nao a
//     malha.  Quem verifica a malha e' a F3, que hoje EXISTE -- e com o mesmo
//     deslocamento ela reprova, "29 de 50 pontos fora da grade";
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
// LIMITE DA G5, e dizer isto importa porque e' tentador ler mais nela do que ha':
// ela pega DIVERGENCIA entre as duas vias, nao erro CONSISTENTE na numeracao.
// Medido: deslocando todos os ids de 1, a G5 continua verde (escrita e leitura
// erram igual) e quem reprova e' a G2 (id fora do intervalo).  Ignorando o
// deslocamento da arvore, a G2 acusa 4 ids repetidos e o processo MORRE na guarda
// de linha repetida do solver -- duas defesas em cadeia, e sem linha de resumo,
// que e' o que a guarda 3 do driver pega.
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
// TRES DEFEITOS QUE A FAMILIA L PEGOU, e os tres sao de ORACULO -- nenhum era do
// calculo em si:
//
//   1. o teste chamava de "cpu" um vetor calculado com o back-end PADRAO, e com
//      `make GPU=sim` o padrao e' a placa.  A comparacao GPU-contra-CPU comparava
//      GPU contra GPU e dava diferenca EXATAMENTE ZERO.  Zero nao era boa noticia:
//      era o oraculo deixando de existir.  Hoje o caso forca a CPU e RECUSA
//      diferenca zero, dizendo o motivo;
//   2. trocar `GPU=sim` por `GPU=nao` nao reconstruia nada -- a flag nao entrava
//      na dependencia do make, e o binario do build anterior ficava no lugar
//      anunciando a configuracao errada.  Hoje a configuracao e' um SELO em
//      `obj/.config`, conferido no parse: mudou, limpa;
//   3. o kernel fazia minimos quadrados SEM PESO, e o `wls` herdado e' WLSMOVING: o kernel
//      peso `D_i = 1/sqrt(|x-p_i|^2 + EPSMACH)` nas linhas, resolvido por QR.  Com 16 pontos
//      e 3 monomios o sistema tem 13 graus de liberdade, entao os dois ajustes
//      REPRODUZEM OS MESMOS POLINOMIOS com pesos completamente diferentes -- a
//      reproducao nao podia acusar, e nao acusou.  Quem acusou foi a comparacao
//      GPU-contra-CPU: 0,311.  Com a ponderacao no kernel, 8,8e-11.
//
// A tolerancia de 1e-10 nao e' arbitraria: e' o patamar de dois ALGORITMOS
// diferentes (QR no herdado, equacoes normais no kernel) resolvendo o mesmo ajuste.
// E em dispositivo sem fp64 ela e' 1e-4, porque o kernel roda em float -- medido
// 5,8e-07 na Radeon via rusticl.
//
// O DEFEITO DE ORACULO QUE A FAMILIA C EXPOS, e ele atingia DUAS familias: os
// casos de E e de C avaliavam o estencil NO CENTRO das celulas, onde o ajuste
// movel e' degenerado -- peso 1,0 na propria celula e 6,4e-11 nas outras, `rhs`
// de -3,2e-11.  Os casos passavam com 1e-16 SEM EXERCITAR O FECHAMENTO, porque o
// estencil ali e' a identidade.
//
// Quem denunciou foi o CONTRASTE, e nao o verde: alimentar Neumann com a derivada
// de outro campo piorava o erro para so' 1e-10 -- pequeno demais para um dado
// grosseiramente errado.  Corrigido, o mesmo contraste da' 1,9e-01.
//
// Duas guardas ficaram: os alvos sao DESLOCADOS do centro (que e' tambem o uso
// real -- interpolar em faceta, canto, ponto de projecao), e os casos de contraste
// exigem erro grande em VALOR ABSOLUTO, nao apenas uma razao.
//
// LIMITE DA FAMILIA R QUE O TESTE SERIAL NAO ALCANCA, e dizer isto importa: o
// sorteio do refino aleatorio vem de um HASH DA POSICAO da celula, e nao de um
// gerador sequencial, justamente para a decisao NAO depender da ordem de visita
// nem da particao -- a mesma semente tem de dar a mesma malha em np=1 e em np=4.
//
// Mas em SERIAL isso nao e' distinguivel: um gerador sequencial daria a mesma
// malha, porque a ordem de visita e' sempre a mesma.  A clausula R6 afirma
// reprodutibilidade, nao independencia de ordem.  O caso que discriminaria e'
// paralelo, e paralelo ainda nao existe -- quando existir, comparar a malha entre
// np diferentes e' o primeiro teste a escrever.
//
// OS TRES DEFEITOS DA FAMILIA T, e os tres eram de SINAL ou de caminho esquecido
// -- nenhum de formulacao.  Guardados porque o que os achou foi o mesmo oraculo:
//
//   1. a correcao transversal com o sinal dependendo de quem era o grosso.  O erro
//      DOBROU em vez de cair (7,8e-02 contra 3,9e-02 sem correcao), e erro que
//      dobra e' a assinatura de uma correcao invertida -- foi o que levou a refazer
//      a algebra em vez de ajustar numeros.  A conta certa nao depende de quem e' o
//      grosso: `fluxo = coef*u_viz - coef*u_self - coef*delta*(du/dt)` nos dois
//      casos, e o que muda e' ONDE `du/dt` e' avaliado;
//   2. `s2` medido do CENTRO em vez da PAREDE no fechamento generico.  Em malha
//      uniforme deu `h` onde o certo e' `3h/2`, e o caso de regressao do grau 2
//      pegou na hora: 1,97e-01 em vez de 3e-13;
//   3. com parede de Neumann no lado TRANSVERSAL, a correcao era ABANDONADA -- e o
//      residuo voltava a 53,3, exatamente o valor de quando nao havia correcao
//      nenhuma.  O numero repetido foi o que apontou a causa: nao era outro
//      defeito, era o MESMO reaparecendo por um caminho deixado aberto.  Ali o
//      dado de Neumann E' `du/dt` (a normal da parede e' a direcao transversal), e
//      usa-lo e' exato -- com `du/dt = sgn*g`, nao `g`.
//
// LIMITE DA FAMILIA Q, e ele e' o fechamento: para o compacto ser exato ate' grau
// 5 -- o que o interior permite -- o fechamento de parede precisaria de CINCO
// pontos (parede e quatro celulas).  Com quatro, o corte fica em grau 3.  A ordem
// de convergencia NAO sofre com isso (o erro do fechamento e' local e tambem
// O(h^4)); o que se perde e' a exatidao em grau 4 e 5.  Subir o fechamento e' o
// proximo passo natural desta familia, e e' trabalho isolado: a formula de
// Lagrange e' a mesma com um no' a mais.
//
// LIMITES DAS FAMILIAS P e T, declarados: Robin nao montado no operador (so'
// Dirichlet e Neumann), serial, e o PISO do erro e' o residuo do solver -- pedir menos de 1e-13 seria exigir do solver linear algo que
// ele nao promete com tolerancia 1e-14.
//
// LIMITE DA FAMILIA C: a regiao e' consultada por VARREDURA de todas as regioes
// registradas, O(regioes) por faceta de contorno.  Sao poucas, e varrer todas e'
// o que permite RECUSAR sobreposicao em vez de esconde-la -- parar na primeira
// seria mais rapido e dependeria da ordem.
//
// LIMITE DA FAMILIA E, declarado: so' DIRICHLET no construtor;  Neumann e Robin entram pela familia C,
// por celula FANTASMA, e isso resolveu o limite que estava declarado aqui.  E a
// selecao de suporte varre a malha para traduzir id
// em (bloco, folha): O(folhas) por consulta, aceitavel em caminho de montagem e
// nao em laco quente.  Dito para nao ser suposto rapido.
//
// LIMITE DA FAMILIA A: o criterio de "face mais e' contorno" e' `num_vizinhos == 0`,
// que vale em np=1.  Sob particionamento, face sem vizinho LOCAL inclui face de
// FRANJA, e a distincao exige a camada de ghost -- que o t8Flow nao tem.  A posse
// das facetas tera' de ser reexaminada quando a franja entrar.
//
// Tambem ausentes, previstas no doc/projeto-t8flow.md: numeracao global contigua
// conteudo da franja depois do sync na profundidade de PRODUCAO (5)
// e em np=3, que e' o menor np em que pareamento transposto aparece; a guarda
// contra arvore de UMA celula (que nao tem interior, e a franja declarada e' de
// cinco); a banda fina da decisao 4; posse distribuida; e o lado 3D inteiro.
// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------
// D -- DERIVADA DE ORDEM ALTA E QUADRATURA DE GAUSS
//
// Este grupo cobre o esquema `MlsGauss`: derivada normal por ajuste MLS na faceta
// e quadratura de Gauss nos DOIS lugares onde havia regra do ponto medio -- a
// faceta e o volume.  E' o unico esquema do sistema que passa de segunda ordem em
// malha ADAPTADA, e cada afirmacao aqui foi medida, nao deduzida.
// ---------------------------------------------------------------------------

//! D1  A DERIVADA MLS E' EXATA ATE' O GRAU DA ORDEM, inclusive na ordem 4 -- que
//! e' o teto da base propria.  Com `pesos_derivada` de ordem 4 um polinomio de
//! grau 4 tem a derivada reproduzida em erro de maquina, e um de grau 5 nao.
//!
//! O CONTRASTE E' A METADE QUE IMPORTA: sem o grau acima da ordem falhando, o
//! teste passaria com uma base truncada e ninguem saberia.
//
//! COBERTA-POR: test-mls2d:a_derivada_de_ordem_quatro_e_exata_no_grau_quatro

//! D2  A BASE E' ADIMENSIONALIZADA pelo raio do suporte, e sem isso a ordem 4
//! degrada EM SILENCIO.  A matriz normal da ordem 4 em base nao escalada tem
//! condicionamento da ordem de `h^-8`: com `h = 1e-3` isso passa de `1e24` e o
//! ajuste perde todos os digitos -- devolvendo pesos, nao erro.
//!
//! O oraculo e' um suporte MINUSCULO (`h = 1e-3`) onde a derivada tem de sair
//! exata na mesma; em base nao escalada ele falha.
//
//! COBERTA-POR: test-mls2d:a_derivada_de_ordem_quatro_sobrevive_a_suporte_minusculo

//! D3  O REBAIXAMENTO DE ORDEM E' REGISTRADO, nunca silencioso.  Quando o suporte
//! nao sustenta a ordem pedida, o ajuste cai para a ordem que ele sustenta e
//! `Mls::rebaixamentos()` conta -- e os pesos continuam VALIDOS.
//!
//! A versao anterior devolvia um vetor de ZEROS e seguia.  O efeito medido em
//! malha adaptada com ordem 4 foi uma linha da matriz silenciosamente errada: erro
//! de 3,3e-01 num nivel cercado por 1e-10 nos vizinhos.  Oito ajustes em 5172.
//
//! COBERTA-POR: test-mls2d:o_rebaixamento_de_ordem_e_contado_e_os_pesos_continuam_validos

//! D4  O ESQUEMA MlsGauss E' EXATO ATE' O GRAU DA ORDEM MLS em malha UNIFORME --
//! grau 3 com ordem 3, grau 4 com ordem 4 -- e o esquema de fluxo de segunda
//! ordem nao e', no mesmo grau e na mesma malha.
//
//! COBERTA-POR: test-poisson2d:o_mlsgauss_e_exato_no_grau_tres_onde_o_fluxo2_nao_e

//! D5  NA MALHA ADAPTADA O MlsGauss MANTEM A EXATIDAO que o fluxo2 perde.  A
//! clausula T2 registra que o grau 2 DEIXA de ser exato em malha adaptada com o
//! esquema de fluxo; com `MlsGauss` ele volta a ser exato, e o grau 3 tambem.
//!
//! E' a razao de existir do esquema, e o par de medidas na MESMA malha -- mesma
//! fracao, mesma semente -- e' o que separa a afirmacao da esperanca.
//
//! COBERTA-POR: test-poisson2d:na_malha_adaptada_o_mlsgauss_recupera_a_exatidao_do_grau2

//! D6  A QUADRATURA DA FONTE E' PARTE DO ESQUEMA, e nao um detalhe: com o fluxo da
//! faceta em ordem alta e `integral f dA` ainda na regra do ponto medio, o esquema
//! TRAVA em segunda ordem.  Medido no grau 4 com ordem 4 em malha uniforme:
//!
//!     ponto medio   9,58e-05  ordem 1,99
//!     Gauss 2x2     2,52e-14  exato
//!
//! O interruptor `quadratura_fonte` existe para esta medida -- verde que nunca foi
//! vermelho nao prova nada.
//
//! COBERTA-POR: test-poisson2d:a_fonte_no_ponto_medio_trava_o_esquema_em_segunda_ordem

//! D7  O FLUXO DA FACETA DE TRANSICAO 2:1 E' O MESMO VISTO DOS DOIS LADOS, e isso
//! e' consequencia do desenho e nao de um acerto: o suporte do estencil depende
//! apenas do PONTO de quadratura, nao da celula que pergunta.
//!
//! O oraculo calcula os pontos de quadratura das sub-facetas da celula GROSSA e os
//! da faceta da celula FINA por caminhos separados e exige que coincidam, e depois
//! compara os estenceis de derivada neles.  Um erro de indexacao na sub-faceta --
//! o lugar classico -- apareceria aqui e em nenhum dos testes de convergencia,
//! porque um esquema nao conservativo ainda converge.
//
//! COBERTA-POR: test-poisson2d:o_fluxo_da_faceta_de_transicao_e_o_mesmo_dos_dois_lados


// ---------------------------------------------------------------------------
// X -- A FRANJA E O PARALELO
//
// Este grupo cobre o que o caminho SERIAL nao exercita, e por isso e' o grupo
// onde a disciplina de contraste mais importa: em np=1 quase toda afirmacao
// daqui e' trivialmente verdadeira, e um caso escrito sem cuidado passa sem
// medir nada.  Todos os casos deste grupo rodam em np=2 e np=3.
// ---------------------------------------------------------------------------

//! X1  O ID DE UM VIZINHO DE FRANJA NAO E' `primeiro_id() + indice`.  O t8code
//! indexa as folhas de franja logo depois das locais, e a soma do deslocamento
//! deste processo a um indice de franja da' o id de uma folha local que EXISTE --
//! so' que outra.
//!
//! Em np=1 nao ha' franja e a soma acerta sempre; e' por isso que o caminho
//! serial nunca acusaria isto, e por isso o caso roda em np>=2.
//!
//! O oraculo e': a uniao dos ids que cada processo ve como vizinhos coincide com
//! a vizinhanca que o serial ve na mesma malha.
//
//! COBERTA-POR: test-paralelo2d:o_id_do_vizinho_de_franja_nao_e_deslocamento_mais_indice

//! X2  A FRANJA TEM O NUMERO DE CAMADAS PEDIDO, e `exige_camadas` e' idempotente.
//! A franja nativa do t8code tem UMA camada; as demais vem de rodadas de troca,
//! uma camada por rodada.
//!
//! Medido: com uma camada so', o suporte do MLS de ordem 3 saiu com UM ponto
//! contra os 10 exigidos, porque o ponto de quadratura cai numa celula de franja
//! cujos vizinhos o processo nao conhecia.
//
//! COBERTA-POR: test-paralelo2d:a_franja_tem_as_camadas_pedidas_e_exigir_e_idempotente

//! X3  CAMADAS = ANEIS + 1, e a diferenca de um foi MEDIDA, nao deduzida: o ponto
//! de quadratura esta' SOBRE a faceta, logo a celula que o contem pode ser de
//! franja, e dela o suporte ainda precisa abrir os aneis inteiros.
//!
//! Com `aneis` camadas: np=2 truncou 32 suportes e np=3 falhou com
//! `suporte de 9 pontos em 2 aneis, minimo 10 para ordem 3`.
//! Com `aneis + 1`: zero truncamentos em np ate' 4.
//
//! COBERTA-POR: test-paralelo2d:com_camadas_de_menos_o_suporte_trunca_e_com_as_certas_nao

//! X4  BORDA DA FRANJA NAO E' BORDA DO DOMINIO.  As duas dao lista de vizinhos
//! VAZIA, e confundi-las faria o esquema aplicar condicao de contorno no MEIO do
//! dominio -- ao longo de toda a fronteira de particao, e em silencio.
//!
//! `expansivel_do_indice` e' o que distingue, e o caso afirma que nenhuma folha
//! da fronteira de particao e' tratada como contorno: o numero de paredes vistas
//! pelo estencil em np>1 e' o MESMO que em np=1.
//
//! COBERTA-POR: test-paralelo2d:a_borda_da_franja_nao_e_confundida_com_borda_do_dominio

//! X5  AS NORMAS E O RESIDUO REDUZEM.  Sem reducao cada processo devolve o pior
//! erro da SUA particao, e o numero que sai e' menor que o verdadeiro -- num
//! relatorio de verificacao isso e' pior que um erro grande: e' um erro grande
//! que se apresenta como pequeno.
//!
//! O residuo tinha uma segunda falha, do mesmo tipo: ele descartava as colunas
//! da linha que caiam fora da faixa local, de modo que cada celula vizinha a uma
//! fronteira de particao perdia os termos de la' -- e o residuo saia MENOR, isto
//! e', um operador inconsistente se apresentando como consistente.  Era o oraculo
//! mentindo a favor do codigo.
//
//! COBERTA-POR: test-paralelo2d:as_normas_e_o_residuo_sao_os_mesmos_em_qualquer_np

//! X6  A SOLUCAO NAO DEPENDE DO NUMERO DE PROCESSOS.  A afirmacao mais forte do
//! grupo, e a que resume as outras: o mesmo problema em np=1, 2, 3 e 4 da' o
//! mesmo erro ate' o residuo do solver linear -- que legitimamente varia com a
//! particao, porque o precondicionamento varia.
//!
//! Vale para malha UNIFORME e para ADAPTADA, e a adaptada e' o caso que importa:
//! e' nela que a transicao 2:1 cruza fronteira de particao.
//
//! COBERTA-POR: test-paralelo2d:a_solucao_e_a_mesma_em_qualquer_np_na_malha_uniforme
//! COBERTA-POR: test-paralelo2d:a_solucao_e_a_mesma_em_qualquer_np_na_malha_adaptada

//! X7  A REPARTICAO PRESERVA A MALHA e refaz a franja.  O refino local
//! desbalanceia -- medido 2,8x --, e reparticionar devolve o balanco.  O que o
//! caso afirma e' que ela nao muda a MALHA: mesma contagem global de folhas,
//! mesma area total, mesmos niveis minimo e maximo.
//!
//! E que a franja foi refeita: esquecer isso deixaria a franja velha apontando
//! para folhas que mudaram de dono -- e isso nao falha, apenas responde errado.
//
//! COBERTA-POR: test-paralelo2d:a_reparticao_preserva_a_malha_e_refaz_a_franja

//! X8  A CARGA DA TROCA NAO SATUROU.  O canal de troca do t8code e' um array de
//! itens de tamanho FIXO, e a vizinhanca publicada por folha e' cortada quando
//! nao cabe.  Cortar e' inevitavel; cortar em silencio nao e'.
//!
//! Foi a saturacao invisivel que fez quatro camadas de franja ainda deixarem
//! niveis ausentes, e mandou a investigacao atras da camada errada por meia hora.
//! `saturacoes_da_franja()` conta, e esta clausula exige ZERO nas malhas da
//! suite.
//
//! COBERTA-POR: test-paralelo2d:a_carga_da_troca_nao_saturou

//! X9  OS IDS DE FRANJA SEGUEM A ORDEM DO T8CODE, e nao a ordem de id.
//!
//! A documentacao do t8code e' explicita: `leaf_face_neighbors` numera
//! `0 .. num_local-1` para as locais e `num_local .. num_local+num_ghosts-1`
//! para as da camada 1.  Ordenar a franja por id -- que e' determinista, e
//! parecia a escolha limpa -- faz o indice do terceiro fantasma designar OUTRA
//! folha, e `ids_vizinhos` entrega um id que EXISTE, so' que de outra celula.
//!
//! Medido: `grau1` em malha adaptada, que o esquema de fluxo resolve
//! EXATAMENTE, deu 3,97e-13 em np=1 e 1,70e-01 em np=2.  O oraculo e' a
//! correspondencia id -> geometria, que nao pode depender da particao.
//
//! COBERTA-POR: test-paralelo2d:o_mapa_de_id_para_geometria_nao_depende_da_particao


// ---------------------------------------------------------------------------
// W -- O MODULO DE DISCRETIZACAO: diferencas finitas ou volumes finitos
//
// Este grupo cobre a ESCOLHA entre as duas familias, que nao e' de implementacao
// e sim de significado da incognita.  Os solvers de Navier-Stokes que virao --
// um por familia -- serao construidos sobre esta base, e e' aqui que o contrato
// entre eles e' escrito.
// ---------------------------------------------------------------------------

//! W1  A FAMILIA SAI DA VARIANTE, nunca e' informada junto.  Informar as duas
//! abriria a possibilidade de informa-las INCONSISTENTES, e um par
//! (familia, variante) incoerente e' o defeito que o modulo existe para impedir.
//!
//! O caso afirma a tabela inteira: cada variante tem a familia certa, e o nome
//! de cada uma traz o prefixo da familia.
//
//! COBERTA-POR: test-discretizacao2d:a_familia_sai_da_variante_e_o_nome_traz_o_prefixo

//! W2  AS QUATRO CONSEQUENCIAS DA FAMILIA ANDAM JUNTAS: o que a incognita
//! significa, como a fonte entra, como o residuo normaliza, e o que o operador
//! impoe.  Separa-las e' um erro que NAO FALHA -- um operador de volumes finitos
//! com fonte pontual continua convergindo, em ordem reduzida e com a constante
//! errada.  Nenhum teste de convergencia o pega.
//
//! COBERTA-POR: test-discretizacao2d:as_consequencias_da_familia_andam_juntas

//! W3  O LAPLACIANO PONTUAL DO MLS e' exato ate' o grau da ordem, e exige ordem
//! pelo menos 2 -- um ajuste linear nao TEM segunda derivada, e devolver zero
//! seria responder a uma pergunta mal posta.
//
//! COBERTA-POR: test-mls2d:o_laplaciano_pontual_e_exato_ate_o_grau_da_ordem

//! W4  O PESO DO LAPLACIANO E' REGULARIZADO, e sem isso ele degrada em silencio.
//!
//! O peso herdado e' `1/(d^2 + EPSMACH)`.  Quando o alvo COINCIDE com um ponto do
//! suporte -- o que acontece SEMPRE no laplaciano pontual, porque o alvo e' o
//! centro da celula e o centro esta' no suporte -- esse ponto recebe peso da
//! ordem de 1e32 e o ajuste passa a INTERPOLAR.  Para o funcional VALOR isso e'
//! o comportamento certo (clausula C5); para a segunda derivada e' fatal.
//!
//! MEDIDO: com o peso herdado, as ordens 3 e 4 do laplaciano REBAIXAVAM sempre e
//! davam exatamente o resultado da ordem 2 -- o que apareceu, no Poisson, como
//! ordem 2 de convergencia num esquema que deveria dar mais.
//!
//! O contraste e' o caso: com o alvo DENTRO do suporte, ordem 4 tem de reproduzir
//! o laplaciano de um grau 4; sem o regularizador, nao reproduz.
//
//! COBERTA-POR: test-mls2d:o_peso_do_laplaciano_e_regularizado_quando_o_alvo_esta_no_suporte

//! W5  AS DUAS FAMILIAS DAO RESULTADOS DIFERENTES COM TUDO O MAIS IGUAL -- mesmo
//! ajuste, mesma ordem, mesmo suporte, mesmo fechamento de contorno.  A unica
//! diferenca e' integrar contra avaliar, e ela aparece:
//!
//!     malha uniforme, MLS ordem 3, solucao trigonometrica
//!       vf-mlsgauss     ordem 2,83
//!       df-mlspontual   ordem 2,03
//!
//! Afirmar que diferem e' o que impede alguem de tratar a escolha como detalhe de
//! implementacao.
//
//! COBERTA-POR: test-poisson2d:as_duas_familias_diferem_com_tudo_o_mais_igual

//! W6  A INSTABILIDADE DO LAPLACIANO PONTUAL EM MALHA ADAPTADA e' DECLARADA, e
//! nao descoberta em producao.
//!
//! O balanco de fluxos garante a estrutura da matriz por construcao -- cada
//! faceta contribui `+coef` ao vizinho e `-coef` a' diagonal, com `coef > 0` --,
//! e o laplaciano pontual nao garante nada.  Medido em malha adaptada com ordem
//! 3: 20 mil iteracoes sem convergir, erro 1,3e+06.
//!
//! `estavel_em_malha_adaptada()` responde, e a clausula exige que ele responda
//! NAO para a variante pontual: o que se proibe aqui e' a resposta otimista.
//
//! COBERTA-POR: test-discretizacao2d:a_instabilidade_do_pontual_em_malha_adaptada_e_declarada


// ---------------------------------------------------------------------------
// Y -- CUSTO: o que tem de crescer linearmente, e o balanco da particao
// ---------------------------------------------------------------------------

//! Y1  A MONTAGEM DA MATRIZ E' LINEAR no numero de celulas, e afirmar isso e'
//! necessario porque ela JA' FOI QUADRATICA -- num sistema cujo objetivo
//! declarado e' escalonamento.
//!
//! `id_do_ponto` varria todas as folhas (O(n)) e a montagem do esquema MLS a
//! chama uma vez por ponto de quadratura, oito por celula: O(n^2).
//!
//! E O SINTOMA APARECEU COMO MENTIRA SOBRE O PARALELO, nao como lentidao.  A
//! montagem levava 8,74 s em np=1 contra 0,285 s em np=8 -- ganho de 30,7x em
//! oito processos.  Montagem paraleliza SEM comunicacao, logo o teto e' 8x; 30,7x
//! so' e' possivel num algoritmo quadratico, onde dividir `n` por 8 da' 64x.  Esse
//! numero foi atribuido a cache e a precondicionador por tres rodadas de
//! investigacao antes de alguem separar os cronometros de montagem e resolucao.
//!
//! Com o indice espacial a consulta e' O(niveis): 0,679 s em np=1, e o ganho de
//! np=1 para np=8 caiu para 5,2x -- sublinear, como tem de ser.
//!
//! O ORACULO E' A RAZAO ENTRE DOIS NIVEIS: quadruplicar as celulas tem de
//! quadruplicar o tempo, nao multiplica-lo por dezesseis.
//
//! COBERTA-POR: test-poisson2d:a_montagem_cresce_linearmente_com_as_celulas

//! Y2  O DESBALANCEAMENTO DA PARTICAO E' MEDIDO, e nao presumido em nenhuma das
//! duas direcoes.
//!
//! O t8code JA' tenta balancear: a documentacao de `t8_forest_set_balance` diz
//! que as florestas intermediarias sao reparticionadas e que o resultado fica
//! "load-balanced per default", e este codigo o chama com `no_repartition = 0`.
//!
//! So' que as medidas nao fecham.  Todas do MESMO `set_balance(..., 0)`:
//!
//!     refino em CAIXA, nivel 6, np 2/4/8     desbalanceamento 1,00 (perfeito)
//!     refino ALEATORIO, nivel 6, np 8        1,11
//!     `bench/escalabilidade`, outra config   2,8
//!
//! A diferenca nao esta' explicada.  Por isso a reparticao deste sistema e'
//! EXPLICITA e nao automatica: acrescentar um segundo mecanismo automatico em
//! cima de um primeiro que nao se entende daria dois caminhos reparticionando --
//! redundante num caso, insuficiente noutro.
//!
//! O que a clausula exige e' que `desbalanceamento()` RESPONDA, e que
//! `reparticiona()` o leve a 1,0 exatamente.  Sem a medida, quem adapta decide
//! as cegas nos dois regimes.
//
//! COBERTA-POR: test-paralelo2d:o_desbalanceamento_e_medido_e_a_reparticao_o_zera

//! Y3  ARVORE DE UMA CELULA E' RECUSADA NA CONSTRUCAO.  Uma arvore de uma celula
//! NAO TEM INTERIOR -- as quatro facetas dela sao de fronteira --, e nenhum
//! estencil centrado nela encontra vizinho dentro do proprio bloco.
//!
//! O levantamento do t8code registra isto como defeito MEDIDO no sistema
//! antecessor, e com o preco: ele "so' apareceu como escoamento divergindo numa
//! simulacao de meia hora".  Recusar na construcao custa uma excecao imediata.
//
//! COBERTA-POR: test-malha2d:arvore_de_uma_celula_e_recusada

//! Y4  O LOTE DE ESTENCEIS EM GPU SUPORTA VARIOS PROCESSOS na mesma placa, sem
//! corromper resultado.  Tres processos abrindo o mesmo dispositivo dao a MESMA
//! diferenca contra a CPU (8,812e-11), que e' o que se quer: a contencao atrasa,
//! nao altera.
//!
//! O QUE ESTA CLAUSULA NAO AFIRMA, e dizer isso faz parte dela: o lote NAO foi
//! exercitado sobre malha PARTICIONADA.  O arreio usa `sc_MPI_COMM_SELF`, logo
//! cada processo resolve o problema inteiro -- o que verifica a placa com varios
//! clientes, e nao o lote com suporte vindo da franja.
//!
//! A distincao importa porque hoje o Poisson nao usa o lote; quando usar, o
//! suporte virá' da franja e esta clausula precisara' de um caso novo.
//
//! COBERTA-POR: test-lote2d:gpu_bate_com_cpu_nos_pesos


// ---------------------------------------------------------------------------
// H -- O LADO TRIDIMENSIONAL (hexaedros)
//
// As clausulas deste grupo cobrem o que e' PROPRIO de tres dimensoes.  O que
// vale igual nos dois lados nao se repete aqui: os quatro modos de falha do
// paralelo (X1..X9), o rebaixamento registrado (D3) e a regularizacao do peso
// do laplaciano (W4) ja' estao afirmados, e o codigo 3D nasceu com eles
// consertados em vez de redescobri-los.
//
// O ESCOPO DO LADO 3D e' "so' o que o solver vai usar", decidido pelo usuario:
// malha, MLS, estencil, contorno, campos e Poisson de volumes finitos.  Ficam
// de fora os esquemas compactos, o laplaciano pontual e o lote em GPU.
// ---------------------------------------------------------------------------

//! H1  A ORDEM DOS OITO VERTICES DO HEXAEDRO, que o `cmesh` espera em ordem de
//! bits (`x` mais rapido, depois `y`, depois `z`).  Troca-la produz um hexaedro
//! TORCIDO que commita sem erro nenhum.
//!
//! EM DUAS DIMENSOES A INTUICAO DE "SENTIDO ANTI-HORARIO" AJUDA a acertar; com
//! oito vertices ela nao se aplica, e o unico recurso e' verificar a geometria
//! resultante.  O oraculo e' duplo, porque cada metade deixa passar o que a
//! outra pega: o centroide de cada folha cai dentro da caixa do bloco dela, E a
//! soma dos volumes da' a soma dos blocos.  Um hexaedro torcido pode ter o
//! centroide no lugar certo e volume errado.
//!
//! Dois blocos, e o segundo DESLOCADO: com um bloco so' na origem, uma troca
//! que espelhe o cubo em torno do centro passaria.
//
//! COBERTA-POR: test-malha3d:o_centroide_de_cada_folha_cai_dentro_da_caixa_do_bloco

//! H2  NA TRANSICAO 2:1 UMA FACE TEM QUATRO VIZINHOS, nao dois: a face de um
//! hexaedro e' um quadrado, e refinar o vizinho a parte em quatro.  Afirmar o
//! numero e' o que distingue HEX de QUAD tratado como HEX.
//
//! COBERTA-POR: test-malha3d:o_hex_tem_seis_faces_e_a_transicao_parte_cada_uma_em_quatro

//! H3  A FRANJA SE DIMENSIONA PELOS VIZINHOS MEDIDOS, e nao pelo teto teorico.
//!
//! O teto e' 24 -- seis faces vezes quatro --, e o maximo REAL e' 9 em malha
//! adaptada e 6 em uniforme, porque uma celula tem no maximo uma face com
//! vizinhos refinados em cada direcao, e nao todas as seis.
//!
//! MEDIDO ANTES DO CODIGO EXISTIR (`bench/franja3d`), e o resultado mudou a
//! arquitetura: dimensionar pelo teto custaria 2,7 vezes mais memoria -- em
//! tres dimensoes, a diferenca entre 55 KB e 24 KB por folha publicada.
//
//! COBERTA-POR: test-malha3d:o_hex_tem_seis_faces_e_a_transicao_parte_cada_uma_em_quatro

//! H4  A BASE DO MLS EM TRES DIMENSOES TEM 4, 10, 20 E 35 MONOMIOS nas ordens 1
//! a 4 -- contra 3, 6, 10 e 15 em duas.  A contagem e' a conferencia mais barata
//! de que a base e' a certa, e o suporte minimo sai dela.
//!
//! A CONSEQUENCIA E' DE PROJETO, nao de contagem: a ordem 3 pede 20 pontos de
//! suporte onde em 2D pedia 10, o que significa mais aneis, que significa mais
//! camadas de franja -- e em tres dimensoes a franja custa o cubo do raio.
//
//! COBERTA-POR: test-mls3d:a_base_de_tres_dimensoes_tem_4_10_20_e_35_monomios

//! H5  O LAPLACIANO SOMA AS TRES SEGUNDAS DERIVADAS, e o erro de somar duas --
//! o que uma traducao descuidada do codigo 2D produziria -- daria um operador
//! que AINDA soma zero nos pesos e AINDA reproduz a constante.
//!
//! Nenhuma das duas conferencias baratas o pegaria; o que o pega e' o VALOR do
//! laplaciano de um campo quadratico com coeficientes distintos nas tres
//! direcoes.
//
//! COBERTA-POR: test-mls3d:o_laplaciano_soma_as_TRES_segundas_derivadas

//! H6  EM SUPORTE SIMETRICO A ORDEM 2 JA' ALCANCA O GRAU 3 NO LAPLACIANO, e esta
//! clausula existe porque a MEDIDA corrigiu uma expectativa.
//!
//! Esperava-se que a ordem 3 reproduzisse o laplaciano de um cubico e a ordem 2
//! nao.  As duas reproduzem, em 4e-13.  A razao e' real: num suporte simetrico
//! em torno do alvo os monomios IMPARES nao contribuem para o laplaciano no
//! centro -- a contribuicao de `x^3` se cancela contra pesos simetricos.
//!
//! E' SUPERCONVERGENCIA DO ESTENCIL SIMETRICO, e importa porque em malha
//! uniforme o centro de uma celula ve' exatamente esse suporte.  O CONTRASTE
//! existe fora da simetria: com o alvo deslocado, ordem 3 da' 7,3e-13 e ordem 2
//! da' 6,4e-01.
//
//! COBERTA-POR: test-mls3d:o_peso_do_laplaciano_e_regularizado_quando_o_alvo_esta_no_suporte

//! H7  O SUPORTE EM CRUZ NAO SUSTENTA A ORDEM 3 EM TRES DIMENSOES, e o suporte
//! por CAIXA sustenta.  Medido na mesma malha, mesma ordem, mesmo teto:
//!
//!     aneis de face     504 de 512 ajustes rebaixam   (98,4%)
//!     caixa, raio 2       0 de 512                    (0,0%)
//!
//! A razao e' a DISTRIBUICAO e nao a contagem: o crescimento por faces entrega
//! media de 24,4 pontos para 20 monomios -- numero de sobra --, mas quase tudo
//! alinhado aos eixos, e os monomios mistos `xy`, `xz`, `yz` e `xyz` ficam mal
//! determinados.  Em duas dimensoes isso nao aparece, porque ha' menos termos
//! mistos: a cruz de dois aneis da' 13 pontos para 10 monomios e sustenta.
//
//! COBERTA-POR: test-estencil3d:a_caixa_sustenta_a_ordem_tres_onde_a_cruz_nao_sustenta

//! H8  O RAIO MINIMO DA CAIXA SAI DA ALGEBRA, e nao da contagem de pontos.
//!
//! Uma caixa de raio `r` tem `2r+1` valores DISTINTOS de coordenada por direcao,
//! e um polinomio de grau `g` precisa de `g+1` para ser determinado naquela
//! direcao.  Com menos, os monomios ficam linearmente DEPENDENTES nos dados:
//! num cubo 3x3x3, onde `x` so' assume `-h`, `0` e `h`, vale `x^3 = h^2 x`
//! EXATAMENTE, e a matriz normal e' singular.
//!
//! MEDIDO, e o sintoma apontou a causa: parando pela contagem, a caixa de raio 1
//! dava 27 pontos para 20 monomios e 216 dos 512 ajustes rebaixavam.  Os 216
//! eram exatamente `6^3` -- as celulas do INTERIOR, que paravam no raio 1; as de
//! borda cresciam para o raio 2 por falta de pontos e nao rebaixavam.
//
//! COBERTA-POR: test-estencil3d:a_caixa_sustenta_a_ordem_tres_onde_a_cruz_nao_sustenta

//! H9  O PONTO FORA DA FRANJA E' CONTADO; o fora do DOMINIO nao.  As duas
//! respostas negativas de `id_do_ponto` sao diferentes, e confundi-las esconde a
//! que importa -- um suporte truncado pela franja da' pesos diferentes dos do
//! serial, logo um resultado que depende da particao.
//!
//! ESTA CLAUSULA EXISTE PORQUE A PRIMEIRA VERSAO DESCARTAVA OS DOIS EM SILENCIO,
//! e a medida de quantas camadas bastam dava ZERO truncamento com qualquer
//! numero -- bom demais, e portanto sintoma.  Com a distincao, a medida ficou:
//!
//!     camadas       2      3      4      5
//!     truncados   150     48      2      0      (np=3, nivel 3)
//!     KB/folha    5,2   26,6   83,0  200,7
//
//! COBERTA-POR: test-estencil3d:o_ponto_fora_da_franja_e_contado_e_o_fora_do_dominio_nao

//! H10  A FRANJA PUBLICA A CAIXA GEOMETRICA, e nao a vizinhanca topologica.
//!
//! O SUPORTE PRECISA DE UMA E A FRANJA PROPAGAVA A OUTRA, e sao coisas
//! diferentes: a celula diagonal `(2,2,2)` esta' a raio geometrico 2 e a SEIS
//! faces de distancia.  Propagando por busca em largura, cobri-la exigia cinco
//! camadas e 200,7 KB por folha publicada.
//!
//! UMA RODADA BASTA quando se publica a caixa, e a razao fica clara olhando de
//! quem publica para quem recebe: o dono de um fantasma meu tem como LOCAIS
//! justamente as celulas que me faltam.  A caixa de raio `R` dele cobre a
//! minha, porque a celula mais distante da minha, `(R,R,R)`, esta' a distancia
//! `max(R-1, R, R) = R` do fantasma `(1,0,0)`.
//!
//! MEDIDO, np=3, nivel 3:
//!
//!     camadas          2      3      5
//!     topologica     150     48      0    truncados
//!     geometrica      69      0      0
//!     KB/folha       4,3   19,6  114,0
//!
//! Zero truncamento passou de 5 camadas e 200,7 KB para 3 camadas e 19,6 KB --
//! DEZ VEZES menos memoria.
//
//! COBERTA-POR: test-estencil3d:a_franja_geometrica_cobre_o_suporte_com_raio_mais_dois

//! H11  SATURACAO E' PONTO DESCARTADO, e nao ultimo slot usado.
//!
//! A primeira versao contava `slot >= MAXV` depois do laco, o que marcava
//! saturacao sempre que a caixa preenchia exatamente os slots -- 216 falsos
//! positivos com raio 1, onde a caixa tem 27 pontos e o buffer tem 27.
//!
//! UM CONTADOR QUE ACUSA QUANDO NADA SE PERDEU DESACREDITA O ORACULO INTEIRO:
//! quem o le' passa a ignora-lo, e com ele a saturacao verdadeira.
//
//! COBERTA-POR: test-estencil3d:a_franja_geometrica_cobre_o_suporte_com_raio_mais_dois

//! H18  A MASCARA DE FACES SO' AFIRMA O QUE COUBE INTEIRO, e e' escrita DEPOIS
//! da carga util -- nunca antes.
//!
//! A folha publicada na franja leva, alem da geometria, a lista de vizinhos por
//! face e uma MASCARA de quais faces o remetente conhece.  A mascara existe
//! porque lista vazia e' ambigua: contorno do dominio, ou "nao me contaram"?
//!
//! A FORMA ERRADA escreve a mascara inteira primeiro e trunca a carga por
//! capacidade.  A folha entao viaja afirmando conhecer uma face cuja lista
//! ficou para tras -- e lista vazia com bit ligado e' exatamente a assinatura
//! de uma PAREDE DO DOMINIO.  O fechamento aplica Dirichlet no MEIO do dominio.
//!
//! MEDIDO em tres dimensoes, antes do conserto: `CAMPOS` era `8 + max_viz*2`,
//! mas a mascara ocupava o campo 8, de modo que cabiam `max_viz - 1` vizinhos.
//! Numa malha uniforme, onde a folha tem seis vizinhos, o SEXTO nao cabia -- e
//! o sexto e' a face 5, a ultima do laco.  Resultado: a linha interior somava
//! -1,907e-02 onde um laplaciano tem de somar zero, e o erro contra a solucao
//! exata mudava 11% com o numero de processos, sem nada acusar.
//!
//! UMA FACE PARCIAL E' PIOR QUE UMA FACE DESCONHECIDA: desconhecida, o
//! fechamento se RECUSA a aplicar condicao e o esquema para; parcial, ele
//! aplica a condicao errada e o resultado sai plausivel.
//
//! COBERTA-POR: test-poisson3d:o_operador_nao_depende_da_particao

//! H19  AS RODADAS DE TROCA DA FRANJA SAO ITERADAS ATE' COBRIR, com parada
//! MEDIDA e COLETIVA -- e nao fixadas em uma.
//!
//! A franja nativa do t8code e' `T8_GHOST_FACES`, e a documentacao de
//! `t8_forest_set_ghost` diz que so' esse tipo e' suportado.  Logo UMA RODADA
//! SO' ALCANCA VIZINHO DE FACE: o processo diagonalmente oposto numa quina nao
//! e' vizinho de face de ninguem que fale comigo, e a celula dele nunca chega.
//!
//! AUMENTAR AS CAMADAS NAO RESOLVE, e foi o que confundiu o diagnostico: a
//! informacao nao precisava ir mais LONGE, precisava ir por outro CAMINHO.  Com
//! 2, 3 e 4 camadas o truncamento ficava preso em 192 suportes, todos nas
//! celulas que tocam o centro do dominio -- onde as oito particoes se encontram.
//!
//! CADA RODADA AVANCA UM SALTO DE POSSE.  Medido, nivel 3: np=2 pede 1 rodada,
//! np=4 pede 2, np=8 pede 3.  A parada e' a reducao global de "quantos pontos
//! da caixa ainda faltam": um numero fixo de rodadas estaria certo para esta
//! particao e errado para a proxima.
//!
//! E O ESQUEMA PEDE `raio + 2` CAMADAS, nao `raio + 1`: a caixa do suporte e'
//! centrada na celula que contem o PONTO DE QUADRATURA, que fica sobre uma
//! faceta e pode cair na celula vizinha.  Com `+1` o esquema sobrevivia por
//! acidente -- a publicacao em torno dos fantasmas cobria o excedente sem que
//! nada afirmasse isso.
//
//! COBERTA-POR: test-estencil3d:a_franja_geometrica_cobre_o_suporte_com_raio_mais_dois

//! H12  A CONTAGEM DE FACETAS EM TRES DIMENSOES E' `N^2 (N+1)` POR DIRECAO --
//! `N^2` colunas de `N+1` planos --, contra `N(N+1)` em duas.
//!
//! E' a conferencia mais barata da regra de propriedade: a faceta pertence a'
//! celula para a qual e' a face "menos", e exatamente a uma.  O numero errado
//! que a traducao descuidada produziria e' `6 N^3` -- "cada celula tem suas seis
//! facetas" --, e ele conta cada faceta interna DUAS vezes.
//!
//! Contar pega o erro de propriedade antes que ele apareca como divergente que
//! nao fecha, que e' onde apareceria sem esta clausula.
//
//! COBERTA-POR: test-campos3d:a_contagem_de_facetas_segue_a_regra_de_propriedade
//! COBERTA-POR: test-campos3d:cada_faceta_tem_exatamente_um_dono

//! H13  A MEDIDA DE UMA FACETA E' UMA AREA, e nao um comprimento.
//!
//! E' A TROCA QUE MAIS FACILMENTE PASSA DESPERCEBIDA numa traducao de duas
//! dimensoes: la' a medida e' `h[1-d]`, um lado so'; aqui e' o produto dos DOIS
//! lados transversais.  Um fluxo com a dimensao errada AINDA CONVERGE -- para a
//! equacao errada --, e nenhum teste de convergencia o pega.
//!
//! O ORACULO E' A SOMA: as facetas de uma direcao cobrem o dominio `N+1` vezes
//! (os `N+1` planos perpendiculares a ela), cada plano com area 1, logo a soma
//! das medidas e' exatamente `N+1`.  Com a medida errada daria `(N+1) N` -- um
//! fator `N`, que CRESCE com o refino.
//
//! COBERTA-POR: test-campos3d:a_medida_da_faceta_e_uma_AREA_e_nao_um_comprimento

//! H14  O TENSOR TEM NOVE COMPONENTES e o campo vetorial TRES, e as componentes
//! sao independentes.
//!
//! O indice de duas dimensoes era `2*lin + col`; manter aquele em 3x3 faria
//! componentes COLIDIREM -- `(1,1)` e `(0,3)` cairiam no mesmo lugar, e `(2,*)`
//! sairia do vetor.  O caso escreve um valor distinto nas nove e confere todas.
//!
//! No vetorial, o risco e' o laco que nao acompanha a dimensao: ele deixaria a
//! componente `z` sem inicializar, em silencio.
//
//! COBERTA-POR: test-campos3d:o_tensor_tem_nove_componentes_e_o_indice_recusa_o_que_nao_existe
//! COBERTA-POR: test-campos3d:o_campo_vetorial_tem_uma_componente_por_direcao_de_faceta

//! H15  O POISSON 3D E' EXATO ATE' O GRAU DA ORDEM DO AJUSTE, em malha uniforme
//! e em ADAPTADA, e o esquema de fluxo nao e'.  Medido:
//!
//!                       vf-fluxo2    vf-mlsgauss
//!     grau 1 uniforme    1,40e-13          --
//!     grau 2 uniforme    8,95e-03     7,73e-13
//!     grau 2 ADAPTADA    1,55e-01     7,55e-12
//!     grau 3, ordem 3         --      1,16e-13
//!
//! A linha da malha adaptada e' a razao de existir do esquema: dez ordens de
//! grandeza, e e' onde a face de um hexaedro e' partida em QUATRO.
//!
//! O ORACULO E' A EXATIDAO POLINOMIAL e nao a ordem de convergencia: a ordem
//! pode estar certa com a constante errada; a exatidao nao.
//
//! COBERTA-POR: test-poisson3d:o_mlsgauss_e_exato_ate_o_grau_da_ordem_do_ajuste
//! COBERTA-POR: test-poisson3d:o_mlsgauss_de_ordem_tres_e_exato_no_grau_tres
//! COBERTA-POR: test-poisson3d:na_malha_adaptada_o_mlsgauss_mantem_a_exatidao

//! H16  A QUADRATURA EM TRES DIMENSOES E' 2x2 NA FACETA e 2x2x2 NO VOLUME --
//! quatro pontos num retangulo e oito num hexaedro, contra dois e quatro em
//! duas dimensoes.
//!
//! E A FONTE POR GAUSS NO VOLUME NAO E' DETALHE: em duas dimensoes, com o fluxo
//! da faceta ja' exato e a fonte na regra do ponto medio, o esquema TRAVAVA em
//! segunda ordem.  A regra do ponto medio e' erro `O(h^2)` RELATIVO em
//! `integral f dV`, e no esquema de segunda ordem ele se escondia atras do erro
//! do proprio fluxo.
//!
//! Na transicao 2:1 a face e' coberta por QUATRO sub-facetas, e o fluxo total e'
//! a soma.  Integrar a face inteira de uma vez daria o fluxo MEDIO sobre ela,
//! que nao e' a soma dos quatro quando a derivada varia.
//
//! COBERTA-POR: test-poisson3d:na_malha_adaptada_o_mlsgauss_mantem_a_exatidao

//! H17  O LADO 3D RECUSA DISCRETIZACAO QUE NAO SEJA DE VOLUMES FINITOS, e a
//! recusa e' explicita e com a lista dos validos -- nao um numero produzido em
//! silencio.
//!
//! O escopo foi decidido pelo usuario: "so' o que o solver vai usar".  Os dois
//! esquemas compactos recusam malha adaptada (sao contraste, nao ferramenta) e o
//! laplaciano pontual e' instavel nela.  Eles entram quando houver medida que os
//! sustente em tres dimensoes.
//
//! COBERTA-POR: test-poisson3d:o_poisson_3d_recusa_discretizacao_que_nao_e_de_volumes_finitos

// ---------------------------------------------------------------------------
// S -- O AVANCO SEMI-IMPLICITO (o termo viscoso na matriz)
// ---------------------------------------------------------------------------

//! S1  O OPERADOR DE DIFUSAO TEM DUAS IMPLEMENTACOES DE PROPOSITO -- avaliacao
//! explicita e pesos de matriz --, e elas concordam em precisao de maquina em
//! malha UNIFORME e ADAPTADA, com `mu` variavel e com os TRES tipos de
//! contorno.
//!
//! Colapsa-las numa so' trocaria a ordem das somas de ponto flutuante, e com
//! ela os ultimos bits de todo caso monofasico ja' gravado -- o criterio de
//! "o caminho geral reproduz o particular bit a bit" e' o que pegou duas
//! refatoracoes minhas mentindo, e perde-lo e' troca ruim.  A gemea testada e'
//! o preco; e' o mesmo arranjo de `amr_tabela` contra `amr_tabela_exaustiva`.
//!
//! O CASO REFINA ATE' A PAREDE de proposito: o defeito analogo anterior -- o
//! gradiente transversal com duas implementacoes -- concordava no interior e
//! discordava NA PAREDE, e so' apareceu quando o refino chegou la'.
//
//! COBERTA-POR: test-ns2d:o_operador_de_difusao_concorda_pelos_dois_caminhos

//! S2  O FECHAMENTO DE CONTORNO E' AFIM NA PROPRIA INCOGNITA (`A u_f + B`), e o
//! implicito o trata assim -- nao como valor conhecido no lado direito.
//!
//! Jogar o fantasma inteiro no lado direito trata um NEUMANN como Dirichlet com
//! valor velho.  Nao falha: converge, da' numero, e a condicao passa a valer
//! outra coisa.  Nos lados de deslizamento livre da bolha de Hysing a condicao
//! tangencial E' Neumann, entao isto nao e' hipotese remota.
//
//! COBERTA-POR: test-ns2d:o_operador_de_difusao_concorda_pelos_dois_caminhos

//! S3  OS DOIS AVANCOS CONVERGEM PARA O MESMO CAMPO quando `dt` cai, e DIFEREM
//! quando `dt` e' grande.
//!
//! AS DUAS METADES SAO NECESSARIAS.  So' a convergencia nao distingue um
//! semi-implicito correto de um que nao faz nada -- este ultimo daria diferenca
//! zero em qualquer `dt`.  So' a diferenca nao distingue um esquema correto de
//! um errado.
//
//! COBERTA-POR: test-ns2d:os_dois_avancos_temporais_convergem_para_o_mesmo_campo

//! S4  O SEMI-IMPLICITO SOBREVIVE ACIMA DO LIMITE VISCOSO EXPLICITO, e o
//! explicito NAO.
//!
//! O ganho declarado, medido, com o contraste junto: se os dois sobrevivessem,
//! o limite nao seria o que se diz e o ganho seria afirmacao vazia.  O que
//! SOBRA e' o limite capilar -- propriedade da tensao superficial, nao da
//! difusao, e implicitar o viscoso nao o remove.
//
//! COBERTA-POR: test-ns2d:o_semi_implicito_e_estavel_acima_do_limite_viscoso_explicito

// ---------------------------------------------------------------------------
// V -- O NUCLEO GEOMETRICO DO VOF (PLIC)
// ---------------------------------------------------------------------------

//! V1  A AREA CORTADA DA CELULA PELA RETA SAI IGUAL POR FORMA FECHADA E POR
//! RECORTE, em precisao de maquina, para qualquer normal e qualquer razao de
//! aspecto.
//!
//! A forma fechada reduz ao quadrado unitario canonico e inverte uma quadratica
//! por pedacos; o recorte faz Sutherland-Hodgman contra o semiplano e integra
//! pelo sapateiro.  NAO PARTILHAM ALGEBRA, e e' isso que faz a concordancia
//! afirmar alguma coisa.
//!
//! A VARREDURA TEM DE INCLUIR O REGIME ESTREITO, e nao so' o caso bonito.
//! Medido por mutacao deliberada (trocar `t <= a` por `t <= b` no primeiro
//! ramo): o defeito aparece em `n = (1, 1e-14)` com celula `(1; 0,125)` e
//! NAO aparece em normal a 45 graus com celula quadrada.  Varrer so' o caso
//! bonito passaria com metade dos ramos errados.
//
//! COBERTA-POR: test-vof2d:a_area_concorda_pelas_duas_vias

//! V2  A DISTANCIA E' A INVERSA EXATA DA AREA, e isso e' o argumento de
//! correspondencia com o HiGFlow.
//!
//! A area e' ESTRITAMENTE DECRESCENTE em `d` no intervalo ativo, logo a inversa
//! e' UNICA.  Quem reproduz a area do HiGFlow e inverte exatamente TEM de dar o
//! `d` dele -- sem precisar linkar o codigo dele para comparar numero a numero.
//! A ida e volta e' o que sustenta esse argumento.
//
//! COBERTA-POR: test-vof2d:a_distancia_e_a_inversa_exata_da_area
//! COBERTA-POR: test-vof2d:a_distancia_concorda_com_a_bisseccao

//! V3  A CONVENCAO DE LADO E' AFIRMADA CONTRA NUMERO CONHECIDO A MAO, e nao so'
//! contra a outra via.
//!
//! Duas vias podem concordar na MESMA convencao errada: trocar o lado daria
//! concordancia perfeita e area complementar, e os casos V1 e V2 passariam.
//! Com `n = +y` e `d = 0,3` numa celula unitaria a area e' 0,2 -- a da regiao
//! `y >= 0,3`, para onde a normal aponta.  Com o lado trocado daria 0,8.
//
//! COBERTA-POR: test-vof2d:os_casos_que_tem_area_conhecida_a_mao

//! V4  NORMAL NULA E' RECUSADA, e AREA FORA DE `[0, h0 h1]` TAMBEM.
//!
//! Normal nula nao tem lado para apontar -- e e' exatamente o que o estimador
//! devolve numa celula PURA, onde `grad f` e' zero.  Area zero, ao contrario, e'
//! legitima e passa.
//!
//! A RECUSA DA AREA FORA DO INTERVALO E' DELIBERADA.  A fracao volumetrica SAI
//! de `[0,1]` em adveccao direcional partida -- e' conhecido e e' parte do
//! metodo --, mas o ajuste pertence a quem advecta e tem de ser CONTADO.
//! Saturar aqui apagaria a unica evidencia de que o transporte passou do
//! limite.
//
//! COBERTA-POR: test-vof2d:a_normal_nula_e_a_area_fora_do_intervalo_sao_RECUSADAS

//! V5  A FRACAO INICIAL E' EXATA PARA INTERFACE PLANA, e isso e' verificado
//! CRUZANDO com o PLIC -- duas implementacoes que nao partilham nada.
//!
//! Com a funcao de forma AFIM a interface dentro da celula E' uma reta, logo a
//! marcha de quadrados nao aproxima: tem de dar o mesmo que `plic_area`, que
//! chega la' por reducao canonica sem amostrar canto nenhum.  Duas vias exatas
//! concordando afirma mais que qualquer teste de ordem.
//
//! COBERTA-POR: test-vof2d:a_fracao_inicial_e_EXATA_para_interface_plana

//! V6  A SUBDIVISAO 16x16 NAS CELULAS CORTADAS E' PARTE DO METODO, e o ganho e'
//! MEDIDO contra o caso sem ela.
//!
//! O envoltorio do HiGFlow (`example2d_VOF/ns-example-2d.c`) acrescenta essa
//! subdivisao por cima de `get_fracvolN`, e sem ela a interface fica
//! linearizada na celula INTEIRA: para uma curva de curvatura `1/R` isso erra
//! `O(h^3/R)` por celula cortada.  Medido no circulo de Hysing: 8,1e-06 em
//! `h = 0,025` e 5,1e-07 em `h = 0,00625`, caindo como `h^2`.
//!
//! O LIMIAR VEM DE NUMERO MEDIDO, nao escolhido: duas ordens abaixo da deriva
//! de massa que o proprio metodo tem ate' `t = 3` (2,44e-03 no FT do HiGFlow).
//! A primeira versao deste caso exigia 1e-9 -- numero que eu nao derivei de
//! nada -- e falhou.
//!
//! IMPORTA PORQUE O VOLUME INICIAL E' O DENOMINADOR DE TUDO: area de
//! referencia, circularidade e centro de massa partem dele, e um erro ali
//! desloca o benchmark inteiro sem aparecer como erro.
//
//! COBERTA-POR: test-vof2d:a_subdivisao_e_o_que_faz_o_volume_inicial_ficar_certo

//! V7  O CASO AMBIGUO DA DIAGONAL E' RECUSADO, e nao resolvido por convencao.
//!
//! Com dois cantos dentro na diagonal, a celula tanto pode conter duas quinas
//! separadas quanto uma garganta que as liga -- os quatro cantos nao
//! distinguem, e as duas tem areas diferentes.  O HiGFlow aborta com "malha
//! pouco refinada"; aqui lanca, porque quem chama pode estar numa SUB-celula e
//! ai' a saida e' subdividir mais.  Dar um numero seria ESCOLHER, e a escolha
//! ficaria invisivel.
//
//! COBERTA-POR: test-vof2d:o_caso_ambiguo_da_diagonal_e_RECUSADO

// ---------------------------------------------------------------------------
// O -- O PAR DIVERGENCIA/GRADIENTE
//
// A letra O, e nao C: `C` ja' nomeia as CONDICOES DE FRONTEIRA.  Esta colisao
// existiu e a suite nao reclamava -- contava as duas e dizia o total certo.
// Quem a pegou foi a guarda contra nome repetido, posta depois de a mesma
// coisa acontecer com `N`.
// ---------------------------------------------------------------------------

//! O1  A DIVERGENCIA E O GRADIENTE SAO ADJUNTOS EM MALHA UNIFORME (precisao de
//! maquina) e falham em SEGUNDA ORDEM na adaptada.
//!
//! AFIRMA SOBRE O PAR, e nao sobre cada operador.  Os dois tem clausulas
//! proprias -- exatidao em campo linear, fechamento da projecao -- e podem
//! passar com o par desencontrado: a adjuncao e' o que diz se a troca discreta
//! de energia entre velocidade e pressao se conserva, e nenhum teste de um
//! operador so' a ve'.
//!
//! `<Du, p>_V = -<u, Gp>_W` com `W_f = medida_f x dist_f` -- o volume de
//! controle dual.  O peso nao e' escolha: sai de igualar os coeficientes de
//! `u_f p_i` nas duas formas bilineares, e com outro peso a igualdade nao sai
//! nem no caso uniforme.
//!
//! NA ADAPTADA NAO VALE, E ERA PREVISIVEL: o gradiente leva correcao
//! transversal na transicao 2:1 -- o centro da celula grossa esta' deslocado da
//! linha da faceta -- e a divergencia nao tem termo que a case.  Sem a correcao
//! o Poisson tinha ordem NEGATIVA, entao ela nao e' negociavel.
//!
//! O QUE IMPORTA E' A ORDEM, NAO O VALOR.  Medido: 8,47e-06 / 1,78e-06 /
//! 4,18e-07 nos niveis 4/5/6, razoes 4,75 e 4,26.  Segunda ordem, logo
//! perturbacao que some e nao defeito estrutural.  Um valor pequeno num nivel
//! so' nao distinguiria "perturbacao que some" de "defeito pequeno que fica",
//! e e' por isso que a assercao e' a razao entre niveis.
//!
//! NO METODO SEGREGADO ISSO E' INVISIVEL: o Poisson da projecao e' montado
//! literalmente como `D G`, e e' consistente por construcao.  A clausula
//! existe porque a propriedade caracteriza a DISCRETIZACAO, e porque num
//! sistema acoplado -- se vier a haver um -- seria dela que sairia a escolha do
//! Krylov: a matriz seria nao simetrica nessa mesma ordem.
//
//! COBERTA-POR: test-ns2d:a_divergencia_e_o_gradiente_sao_adjuntos

// ---------------------------------------------------------------------------
// U -- O NAVIER-STOKES EM TRES DIMENSOES
//
// A letra U, e nao N: `N` ja' nomeia as clausulas do NUCLEO HERDADO, la' em
// cima.  Duas clausulas com o mesmo nome tornam ambigua qualquer citacao --
// e o arreio contava as duas e dizia 123/123, sem reclamar.  Hoje reclama.
// ---------------------------------------------------------------------------

//! U1  A DIVERGENCIA 3D E' EXATA PARA CAMPO LINEAR, e o zero vem com o par que
//! NAO e' zero.
//!
//! A divergencia de volumes finitos e' a regra do ponto medio sobre a
//! superficie, exata para campo linear: um campo solenoidal da' zero de
//! maquina em qualquer resolucao, e um de divergencia conhecida da' o valor.
//! Ordem certa com constante errada existe; exatidao nao.
//!
//! O PAR E' NECESSARIO.  "A divergencia deu 1e-16" nao distingue um operador
//! correto de um que devolve zero sempre -- e um operador que descarta faces
//! devolve justamente isso em metade dos casos.  Dois campos de divergencia
//! nula (um deles CISALHANTE, que varia ao longo da face que o carrega) e um
//! de divergencia 3.
//
//! COBERTA-POR: test-ns3d:a_divergencia_e_exata_em_malha_uniforme

//! U2  O BALANCO DE FLUXOS CONSULTA AS SEIS FACES, DONAS OU NAO.
//!
//! `DominioFacetas::id_faceta` devolve -1 quando a folha nao e' dona da
//! faceta, e isso e' informacao e nao falha: e' assim que o chamador sabe que
//! ela pertence ao vizinho.  Mas o balanco precisa das seis, e as que nao sao
//! dela vem do vizinho pela face oposta.
//!
//! MEDIDO AO ESQUECER DISTO, na primeira execucao do caso: campo de
//! divergencia 3 devolvia 18, e campo solenoidal devolvia 10 em vez de zero.
//! Nao e' erro pequeno -- e' o operador somando um subconjunto arbitrario de
//! faces --, mas num oraculo de ORDEM seria facil confundir com malha grossa.
//! E' por isso que o primeiro incremento do solver 3D e' a divergencia
//! sozinha, com oraculo de exatidao: o defeito estava na varredura de facetas,
//! e encontra-lo depois significaria procura-lo atraves de quatro operadores.
//
//! COBERTA-POR: test-ns3d:a_divergencia_e_exata_em_malha_uniforme

//! U3  A MEDIDA DA FACETA E' UMA AREA, e o caso CUBICO nao mostra isso.
//!
//! Numa celula cubica `area/volume` e `comprimento/area` dao o MESMO numero:
//! trocar a medida da faceta por engano passaria despercebido, e o fluxo com
//! dimensao errada ainda converge -- para a equacao errada.  O caso usa um
//! dominio 1 x 2 x 4 com contagens diferentes nas tres direcoes e um campo de
//! divergencia `2 + 5 - 1 = 6`.
//
//! COBERTA-POR: test-ns3d:a_medida_da_faceta_e_uma_AREA_e_o_caso_cubico_nao_mostra

//! U4  A TRANSICAO 2:1 E O PARALELO SAO RECUSADOS NA CONSTRUCAO, e a recusa
//! DIZ por que tres dimensoes nao e' duas com um indice a mais.
//!
//! Recusar na construcao e nao no primeiro passo: quem monta uma malha
//! adaptada e so' descobre no passo mil ja' gastou a montagem, e pode ter
//! visto numeros no caminho.
//!
//! E A MENSAGEM E' PARTE DA CLAUSULA.  O caso exige que ela diga o que muda:
//! a face da celula grossa e' coberta por QUATRO facetas finas, nao duas, e
//! cada faceta tem DUAS direcoes transversais, nao uma.  Uma recusa generica
//! passaria num teste de "lancou?" e nao ensinaria nada a quem a encontrasse.
//!
//! ESTA CLAUSULA JA' EXIGIU "OITO CELULAS NA QUINA", e era falso: o ponto onde
//! `mu` e' avaliado fica sobre uma ARESTA -- dois planos de face cruzando --,
//! tocada por QUATRO celulas, como a quina em duas dimensoes.  Oito seria um
//! VERTICE, que este termo nao usa.  A clausula, o teste e a mensagem
//! concordavam entre si e estavam os tres errados: concordancia nao e'
//! verificacao quando as tres copias vem da mesma afirmacao nao conferida.
//
//! COBERTA-POR: test-ns3d:a_transicao_2_1_e_o_paralelo_sao_RECUSADOS_na_construcao

//! U5  A PROJECAO ZERA A DIVERGENCIA DE UM CAMPO QUALQUER, e o zero vem com o
//! contraste de que o campo de partida NAO era solenoidal.
//!
//! O oraculo nao precisa de fisica: dado um campo arbitrario, projeta-lo tem
//! de deixar divergencia zero de maquina.  Exercita a cadeia inteira --
//! divergencia, montagem do Poisson, solucao, gradiente -- sem que exista
//! passo de tempo ou termo convectivo.
//!
//! E' O PAR `div`/`grad` QUE FECHA, nao cada um.  Se a projecao usasse um
//! gradiente diferente do que a divergencia do Poisson supos, o campo
//! corrigido nao teria divergencia nula e NADA DIRIA: o solver convergiria
//! igual, e o defeito so' apareceria na fisica, muito depois.  Medido por
//! mutacao -- gradiente da projecao 2% maior: divergencia 1,58 em vez de
//! 1e-11.
//!
//! SEM O CONTRASTE o caso seria vazio: "a divergencia deu 1e-14" nao
//! distingue a projecao de um campo que ja' era solenoidal.
//
//! COBERTA-POR: test-ns3d:a_projecao_zera_a_divergencia_de_um_campo_qualquer

//! U6  A ADJUNCAO `D = -G^T` VALE EM MALHA UNIFORME TRIDIMENSIONAL, e o peso e'
//! `AREA x distancia`.
//!
//! Em duas dimensoes o volume de controle dual era `comprimento x distancia`
//! (clausula O1).  Usar aquele peso aqui nao fecha nem em malha uniforme, e a
//! mutacao mede: tirar a distancia do peso leva o desvio de 1e-13 para 0,649.
//!
//! Afirma sobre o PAR, que nenhum teste de um operador so' ve'.
//
//! COBERTA-POR: test-ns3d:a_divergencia_e_o_gradiente_sao_adjuntos_em_tres_dimensoes

//! U7  A PRESSAO TODA-NEUMANN E' RECUSADA tambem em tres dimensoes.
//!
//! A constante fica no nucleo e o PETSc converge para um elemento qualquer do
//! subespaco -- isso DA' NUMERO, e o numero e' plausivel.  A recusa custou o
//! ajuste de tres casos deste mesmo arquivo, que construiam o solver so' para
//! medir operador e nao tinham pino: a recusa disparou e eles abortaram.  E'
//! o comportamento certo.
//!
//! E A REGIAO DO PINO USA O `h` DA PAREDE, nao um valor digitado: com `h`
//! errado ela pega zero facetas (e volta a recusa) ou duas (e o Poisson fica
//! sobredeterminado).  O mesmo tropeco ja' tinha acontecido no lado 2D.
//
//! COBERTA-POR: test-ns3d:a_projecao_zera_a_divergencia_de_um_campo_qualquer

//! U8  COM `mu` UNIFORME O TERMO VISCOSO E' O LAPLACIANO SOBRE `Re`, EXATO.
//!
//! A diferenca de tres pontos e' exata para quadratica, logo o termo tem de
//! dar `mu lap(u)` em precisao de maquina, em qualquer resolucao.  Isto
//! verifica a estrutura de fluxo nas TRES direcoes.
//!
//! A QUADRATICA TEM AS TRES SEGUNDAS DERIVADAS DIFERENTES ENTRE SI.  Com
//! coeficientes iguais, somar duas direcoes em vez de tres daria um numero
//! PROPORCIONAL ao certo, e o erro passaria por constante de normalizacao.
//! Medido por mutacao (laco sobre duas direcoes): 1,84 em vez de 1e-11.
//!
//! E O TRANSPOSTO E' EXATAMENTE ZERO sem campo de viscosidade -- e' o que faz
//! o caso monofasico reproduzir o laplaciano bit a bit.
//
//! COBERTA-POR: test-ns3d:com_mu_uniforme_o_termo_viscoso_e_o_laplaciano_EXATO

//! U9  `mu` NA ARESTA E' A MEDIA HARMONICA DE QUATRO CELULAS -- quatro, e nao
//! oito.
//!
//! O ponto onde `mu` e' avaliado para a direcao transversal fica sobre dois
//! planos de face com a terceira coordenada num centro de celula: isso e' uma
//! ARESTA, tocada por QUATRO celulas, como a quina em duas dimensoes.  Oito
//! seria um VERTICE -- tres planos de face --, que este termo nao usa.
//!
//! ESTE PROJETO AFIRMOU "OITO" NUMA CLAUSULA, NUM TESTE E NUMA MENSAGEM DE
//! RECUSA, e os tres concordavam entre si porque os tres vieram da mesma
//! afirmacao nao conferida.  Concordancia nao e' verificacao quando as copias
//! tem a mesma origem.
//!
//! HARMONICA E NAO ARITMETICA: e' ela que preserva a continuidade da tensao
//! cisalhante atraves de interface plana.  O caso EXIGE que as duas medias
//! difiram -- senao nao distinguiria uma da outra.  Medido por mutacao
//! (aritmetica no lugar): 5,5 em vez de 1,818.
//
//! COBERTA-POR: test-ns3d:mu_na_aresta_e_a_media_harmonica_de_QUATRO_celulas

//! U10  O TRANSPOSTO DO TENSOR EXISTE, e so' com `mu` variavel.
//!
//! `(d_e mu)(d_d u_e)`.  E' o termo cuja falta, no lado 2D, custou
//! circularidade 0,8361 contra 0,9012 +- 1e-4 na bolha de Hysing -- 650 vezes
//! a dispersao entre as referencias publicadas.
//!
//! O CONTRASTE E' O CASO: um transposto que devolvesse zero sempre passaria na
//! clausula U8, que exige zero com `mu` uniforme.  Aqui se exige que NAO seja
//! zero com `mu` variavel e derivadas cruzadas nao nulas.
//
//! COBERTA-POR: test-ns3d:o_transposto_existe_e_so_com_mu_variavel

//! U11  O CONVECTIVO E' EXATO PARA CAMPO UNIFORME EM QUALQUER ESQUEMA, e para
//! campo LINEAR em todos menos o upwind -- que da' EXATAMENTE METADE.
//!
//! Forma divergente, `sum_e d(u_e u_d)/dx_e`, somada nas TRES direcoes.  Dois
//! oraculos de forma fechada: campo uniforme da' zero; `u_x = A x` da'
//! `2 A^2 x`.
//!
//! O PRIMEIRO PASSA COM O TERMO TROCADO POR ZERO -- verifica estrutura, nao
//! valor.  O segundo nao passa.
//!
//! E O SEGUNDO NAO VALE PARA O UPWIND.  Ele toma `phi` na FACETA a montante
//! (`A x`) e nao no CENTRO da celula (`A(x +- h/2)`), onde fica a face do
//! volume de controle; a conta da' `A^2 x`, metade do certo.  A clausula
//! exige a METADE, que e' oraculo mais forte do que "nao e' exato" -- e a
//! primeira versao deste caso afirmava exatidao para os quatro e falhou com
//! erro 1,47, que e' `A^2 * 0,75` no ponto mais afastado.  A previsao bateu
//! em todos os digitos.
//
//! COBERTA-POR: test-ns3d:o_convectivo_e_EXATO_para_campo_uniforme_e_linear

//! U12  OS QUATRO ESQUEMAS DIFEREM num campo de grau acima do que reproduzem.
//!
//! Se dessem sempre o mesmo, a clausula U11 estaria verificando estrutura e
//! nao esquema -- qualquer interpolacao serviria.  E' a licao do canal de
//! Poiseuille no lado 2D, onde o convectivo se anulava identicamente e os
//! tres esquemas davam o mesmo numero: ali "o CUBISTA esta' correto" nao era
//! afirmacao verificada.
//!
//! CUBISTA e QUICK so' coincidem na faixa central do NVD, e o campo tem de
//! sair dela em algum lugar -- o caso exige que difiram tambem entre si.
//
//! COBERTA-POR: test-ns3d:os_esquemas_DISCRIMINAM_num_campo_de_grau_maior

//! U13  O REBAIXAMENTO DO CONVECTIVO E' CONTADO, e e' ZERO no miolo.
//!
//! Sem o ponto a montante o esquema de segunda ordem nao existe e cai para
//! upwind.  Junto a' parede isso e' inevitavel; no miolo seria defeito.  A
//! clausula exige as DUAS metades: conta nao nula no dominio inteiro -- senao
//! o contador nao esta' contando -- e ZERO no miolo.
//!
//! IMPORTA MAIS EM TRES DIMENSOES: tres direcoes por faceta em vez de duas,
//! logo mais oportunidades de o montante faltar.  Um esquema que se anuncia
//! de segunda ordem e cai para primeira em parte do dominio entrega menos do
//! que diz, e sem o contador isso e' invisivel.
//
//! COBERTA-POR: test-ns3d:o_rebaixamento_do_convectivo_e_CONTADO

//! U14  AS DERIVADAS A MAO DA SOLUCAO MANUFATURADA SAO CONFERIDAS CONTRA
//! DIFERENCAS FINITAS.
//!
//! O termo fonte e' montado de DOZE formulas escritas a mao -- nove
//! `du_i/dx_j` e tres `dp/dx_i`.  Uma delas errada produz um `f` errado, que
//! produz uma solucao numerica que converge BONITO para a coisa errada: a
//! ordem medida sai certa e o resultado nao vale nada.
//!
//! ISTO NAO E' CALCULAR `f` NUMERICAMENTE, e a distincao e' o ponto.  O codigo
//! usa a formula a mao; a clausula compara a formula com uma diferenca
//! central do proprio `u`.  Calcular `f` por diferencas tornaria o oraculo
//! uma comparacao entre dois erros meus -- e e' por isso que o catalogo de
//! manufaturadas deste projeto nao tem diferenciacao automatica em lugar
//! nenhum.
//!
//! Medido por mutacao (um sinal trocado numa das nove): 9,61 contra 1e-06.
//!
//! A CLAUSULA CONFERE TAMBEM que a manufaturada e' SOLENOIDAL.  Ela vem de um
//! potencial vetorial, logo e' por construcao -- mas "por construcao" ja' se
//! provou frase perigosa neste projeto.
//
//! COBERTA-POR: test-ns3d:as_derivadas_a_mao_batem_com_diferencas_finitas

//! U15  CADA OPERADOR REPRODUZ O TERMO ANALITICO DA MANUFATURADA, e e' essa a
//! camada que localiza o defeito quando a ordem do passo inteiro sai errada.
//!
//! Os operadores tem casos com oraculo de EXATIDAO, e o exemplo de MMS mede a
//! ordem do passo COMPLETO.  Entre os dois nao havia nada -- e quando a ordem
//! saiu 0,58 em vez de 2, nao havia como dizer QUAL operador discorda sem
//! instrumentar o exemplo.
//!
//! O oraculo aqui nao e' exatidao: o campo e' trigonometrico e nenhum esquema
//! o reproduz.  E' CONSISTENCIA -- o erro tem de ser da ordem de `h^2` vezes
//! a escala do termo, e nao da ordem do termo.
//!
//! FOI ESTE CASO QUE ACHOU O DEFEITO: viscoso e gradiente passaram, o
//! convectivo discordou 19,6%.  Medindo a ORDEM de cada um -- e nao o valor
//! contra um limiar escolhido -- o viscoso deu 1,88 e 2,00 enquanto o
//! convectivo deu 0,944 e 0,988.  Primeira ordem com diferencas centrais e'
//! defeito, nao discretizacao.
//!
//! A CAUSA: `u_no_ponto` interpolava so' ao longo da propria direcao.  A
//! velocidade que advecta atraves da face transversal do volume de controle
//! vive numa ARESTA, deslocada meia celula em `d` das facetas que a carregam;
//! interpolar so' em `e` devolve o valor no CENTRO e o erro e' `O(h)`.  Hoje
//! a interpolacao e' trilinear sobre a rede de facetas.  E' a mesma classe de
//! defeito que no lado 2D custou 20% em `u` e 37% em `v` na caixa
//! manufaturada -- em tres dimensoes ha' DUAS transversais, logo duas
//! interpolacoes a esquecer.
//
//! COBERTA-POR: test-ns3d:cada_operador_reproduz_o_termo_da_manufaturada

#endif
