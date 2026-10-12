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
// LIMITE DA FAMILIA E, declarado: so' DIRICHLET.  Neumann precisa de derivada na
// faceta, que muda a LINHA do ajuste e nao so' o lado direito -- e' trabalho
// proprio, nao um `else`.  E a selecao de suporte varre a malha para traduzir id
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

#endif
