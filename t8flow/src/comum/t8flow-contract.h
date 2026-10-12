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

#endif
