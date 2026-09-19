#ifndef HIG_MESH_CONTRACT_H
#define HIG_MESH_CONTRACT_H

// =============================================================================
// O CONTRATO DE MESH
//
// Isto e' a especificacao do que uma implementacao de malha hierarquica tem de
// cumprir para servir ao nucleo do HFlow.  Existe para que uma SEGUNDA
// implementacao -- o t8code e' a proxima -- possa ser escrita contra um alvo
// escrito, e verificada sem depender de referencia gravada.
//
// Nao ha' despacho virtual aqui, e isso e' decisao MEDIDA, nao omissao.  Ver
// "A forma do contrato", abaixo.
//
// NENHUM .c INCLUI ESTE CABECALHO, e e' de proposito: ele nao declara funcao, ele
// declara GARANTIAS.  Quem o consome e' o ci/run_higtree_tests.py, que le as
// clausulas daqui e exige que cada uma tenha caso passando -- uma clausula que
// perdeu o teste sai como "SEM TESTE QUE RODE", uma cujo teste falhou sai como
// "REPROVADA".  Fica em src/ e nao em doc/ porque envelhece junto com o codigo que
// restringe, e porque a suite falha quando ele e o codigo divergem.
// =============================================================================
//
//
// ------------------------------------------------------------------ A FORMA
//
// O contrato tem DOIS NIVEIS, porque as duas metades tem frequencias de chamada
// que diferem em uma ordem de grandeza e pedem formas diferentes.
//
// Medido em higflow/src, 19 de setembro de 2026:
//
//     PRODUCAO   ~65 chamadas, nenhuma em laco de celula.  lb_create,
//                lb_add_input_tree, lb_calc_partition, sd_add_higtree,
//                sd_add_fringe_higtree, psd_create, psd_synced_mapper.
//                Executadas uma vez por dominio, na montagem e a cada
//                adaptacao de malha.
//
//     CONSULTA   907 chamadas, 90% DENTRO de laco sobre celulas ou facetas.
//                hig_get_center (228 em laco), sd_get_domain_celliterator
//                (237), hig_get_cid (190), hig_get_delta (122),
//                sfd_get_stencil (120), sd_get_cell_with_point (10).
//                Executadas por celula, por passo de tempo.
//
// Disso decorre a forma:
//
//   O nivel de PRODUCAO e' onde a substituicao acontece, e onde uma interface
//   com despacho dinamico nao custa nada -- sao dezenas de chamadas por
//   simulacao inteira.  E' ali que o t8code entra.
//
//   O nivel de CONSULTA e' contrato COMPORTAMENTAL, nao interface virtual.
//   Torna-lo virtual poria indirecao em 907 sitios de laco quente, e a mesma
//   analise ja' havia rejeitado despacho virtual nos ramos por celula durante a
//   Fase 4 da migracao.  O que amarra este nivel sao as garantias abaixo, e o
//   que as torna reais sao os testes que as verificam -- nao a assinatura.
//
// Uma segunda implementacao, portanto, nao herda de uma classe: ela PRODUZ as
// estruturas que as consultas leem, e e' aceita quando passa nas garantias.
//
//
// ------------------------------------------------------------- AS GARANTIAS
//
// Cada uma e' verificada por um caso da suite (ci/run_higtree_tests.py).  A
// lista abaixo e' a fonte: o driver a compara com os casos que realmente
// rodaram e acusa garantia que ficou sem teste.  Clausula sem teste nao e'
// clausula -- e' comentario.
//
//
// NIVEL DE PRODUCAO
//
//   P1  A particao cobre o dominio exatamente uma vez.  A soma global dos
//       volumes locais e' o volume do dominio, e a soma global das celulas
//       locais e' o total.  Pega celula perdida E celula contada duas vezes,
//       que contagem por rank nao distingue.
//         test-fringe-parallel / particao_cobre_o_dominio_uma_vez
//
//   P2  Havendo vizinho, todo rank recebe franja; nao havendo, nenhum recebe.
//       Reduzido por MIN e por MAX, o que e' mais forte que "alguem tem".
//         test-fringe-parallel / franja_existe_quando_ha_vizinho
//
//   P3  A montagem e' possivel EM SERIE, sem MPI.  Ate' 2026-09-19 nao era: o
//       bloco de faceta so' era preenchido pelo caminho particionado, e metade
//       do contrato nao podia ser exercitada isoladamente.  Interface que so'
//       existe acoplada ao particionamento e' a mais dificil de substituir --
//       que e' justamente o que a segunda implementacao vai fazer.
//         test-facet-domain-serial / montagem_serial_preenche_o_bloco_de_faceta
//
//   P4  O resultado nao depende de COMO o dominio foi dividido.  E' o oraculo
//       que separa "particiona diferente" de "particiona errado", e o t8code
//       vai particionar diferente de proposito.
//         test-partition-independence / uma_arvore_contra_duas
//         test-fringe-parallel       / valor_nao_depende_da_particao
//
//
// NIVEL DE CONSULTA -- geometria e identidade
//
//   C1  centro = (low+high)/2 e delta = high-low, exatos.  Amarra os tres
//       acessores entre si: um deles certo e outro errado nao passa.
//         test-cell-queries / centro_e_delta_coerentes_com_a_caixa
//
//   C2  O centro cai na grade analitica: (k + 1/2) * delta, com k inteiro, em
//       cada direcao e em cada nivel de refino.
//         test-cell-queries / centro_na_grade_analitica
//
//   C3  O identificador e' bijecao sobre [0, n_local): sem repetido, sem
//       buraco, sem id fora da faixa.
//         test-cell-queries / cid_e_bijecao
//
//   C4  O identificador e' estavel: a mesma posicao geometrica devolve o mesmo
//       id, seja alcancada por percurso ou por busca por ponto.
//         test-cell-queries / cid_estavel_entre_percursos
//
//
// NIVEL DE CONSULTA -- iteracao
//
//   C5  O iterador de dominio cobre as celulas locais exatamente uma vez.
//       Afirmado por VOLUME, que pega falta e repeticao; contagem sozinha nao.
//         test-cell-queries / iterador_cobre_o_dominio_uma_vez
//
//   C6  O iterador de dominio NAO inclui a franja.  A franja da' suporte a'
//       interpolacao, nao pertence ao dominio sobre o qual se resolve.
//         test-fringe-support / iterador_de_dominio_ignora_a_franja
//
//
// NIVEL DE CONSULTA -- localizacao por ponto
//
//   C7  A celula devolvida CONTEM o ponto.
//         test-point-location / celula_devolvida_contem_o_ponto
//
//   C8  Ponto sobre face interna vai para a celula de MENOR coordenada: o
//       intervalo e' fechado em cima e aberto embaixo.  Esta convencao era
//       EMERGENTE -- quem a mediu esperava o contrario -- e por isso esta'
//       escrita: convencao que nem quem le o codigo acerta por intuicao nao
//       sobrevive a uma reimplementacao.
//         test-point-location / convencao_de_empate_na_face_esta_fixada
//
//   C9  O desempate nao muda conforme o dominio seja uma arvore ou varias.  Se
//       mudasse, nenhuma comparacao entre decomposicoes valeria.
//         test-point-location / empate_nao_depende_de_como_o_dominio_foi_dividido
//
//
// NIVEL DE CONSULTA -- estencil
//
//   C10 Interpolacao de ordem k reproduz exatamente polinomio de grau <= k, no
//       interior do dominio.
//         test-stencil-value / linear_interior
//
//   C11 O estencil ATRAVESSA salto de nivel maior que 2:1.  A contribuicao de
//       Sousa et al. (2019) e' minimos quadrados moveis em arvore NAO GRADUADA;
//       isto e' requisito da interface, nao detalhe de implementacao.  O
//       balanceamento 2:1 do t8code e' opcional, mas as rotinas de vizinhanca e
//       de ghost dele foram construidas em torno do caso balanceado -- esta
//       clausula e' a que a segunda implementacao tem maior risco de nao
//       cumprir.
//         test-level-jump / reproduz_campo_linear_atravessando_salto_4_para_1
//       MEDIDO no t8code v4.0.0-26.08 (higtree/tests/t8code/probe-c11.c): ele
//       REPRESENTA o salto 4:1 sem balanceamento forcado, e REPORTA os 4
//       vizinhos de face atraves dele.  A premissa sobrevive ao seu teste de
//       maior risco -- o que resta e' produzir, a partir dessa floresta, as
//       estruturas que as consultas leem.
//
//   C12 Fora do dominio, o fechamento usa a parede ATRAVESSADA pela projecao, e
//       nao a mais proxima.
//         test-stencil-selection / fecha_pela_parede_atravessada_e_nao_pela_mais_proxima
//
//   C13 O suporte do estencil ALCANCA a franja, e ela carrega peso.  E' o
//       criterio do que a franja entrega -- nao do tamanho dela, que amarraria
//       a implementacao.  Vale em serie e sob particionamento.
//         test-fringe-support  / com_franja_o_estencil_atravessa_e_carrega_peso
//         test-fringe-parallel / estencil_alcanca_a_franja
//
//   C14 Sobre o contorno (ramo ON_BOUNDARY), a reproducao polinomial vale.
//       Metade do despacho de fechamento passa por aqui, e foi nesse ramo que
//       sobreviveram tres dos sete sitios do defeito de compactacao, verdes por
//       ausencia de teste e nao por estarem certos.
//         test-boundary-path / reproduz_campo_linear_SOBRE_o_contorno
//
//
// ------------------------------------------------- O QUE NAO E' GARANTIA
//
// Registrado para ninguem construir guarda sobre oraculo que nao discrimina:
//
//   Reproducao polinomial NAO detecta franja ausente.  Campo linear sai exato
//   mesmo com suporte so' de um lado, por extrapolacao -- medido nos dois
//   testes de franja.  Quem detecta ausencia e' C13.
//
//   Reproducao polinomial NAO detecta parede de fechamento errada quando o
//   campo e' linear: Lagrange ao longo de qualquer eixo acerta um linear.  Por
//   isso C12 usa campo de grau acima da ordem.
//
//   A suite de exemplos do HiGFlow NAO cobre nenhuma clausula deste contrato.
//   Demonstrado: um defeito no fechamento de contorno viveu anos na base sem
//   mover um digito das 33 execucoes, porque so' aparece quando o valor de
//   contorno VARIA ao longo da parede, e nao-deslizamento e' constante.
//
// =============================================================================

#endif
