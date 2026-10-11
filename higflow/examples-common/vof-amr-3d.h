// *******************************************************************
//  CRITERIO DE REFINO ADAPTATIVO PARA O VOF, 3D
// *******************************************************************
//
// PORTADO do `example2d_VOF` (caso de Hysing), onde o criterio foi desenhado e
// medido.  A estrutura e' a mesma e o que muda e' a dimensao; a maquinaria
// pesada ja' existia em 3D do trabalho do front-tracking e e' REUSADA, nao
// reescrita:
//
//   amr_limiares                 limiares da regra das celulas minimas
//   amr_tabela                   nivel por celula base, por DILATACAO
//   amr_escreve_arquivo          o .amr multinivel
//   malha_adapt_reune_sementes   a reuniao entre ranks
//
// O QUE E' PROPRIO DO VOF, e e' a unica diferenca em relacao ao FT: de onde as
// sementes vem.  No FT saem da GEOMETRIA da frente -- que e' replicada e
// portanto completa em todo rank.  Aqui saem do CAMPO de fracao, que e'
// DISTRIBUIDO: cada rank ve' so' o pedaco de interface que possui.  Por isso a
// reuniao entre ranks NAO e' opcional aqui, e e' o erro que o 2D documenta --
// sem ela o .amr sai refinado apenas em torno do pedaco de quem escreve, e a
// regra das celulas minimas fica violada exatamente onde ninguem esta' olhando.
//
// O PISO DE CINCO CELULAS POR LADO NAO E' NEGOCIAVEL.  O PLIC reconstroi um
// plano por celula a partir da fracao na vizinhanca, e o transporte da fracao
// precisa da mesma vizinhanca: com menos de cinco celulas finas de cada lado da
// interface a reconstrucao nao fecha e o transporte perde massa.  E' correcao,
// nao ajuste, e por isso o piso VENCE a banda fisica -- quando a banda pedida
// da' menos de cinco, a funcao usa cinco e AVISA, em vez de obedecer em silencio
// a um pedido que quebraria o metodo.
//
// A BANDA FISICA, por sua vez, existe para a serie poder convergir: com `celmin`
// FIXO cada resolucao refina largura fisica DIFERENTE -- a mais fina refina
// metade do espaco --, e em 2D foi isso que impediu o estudo de convergencia
// (nenhuma das cinco series monotona, nao-monotonicidade de 12 a 40 vezes o piso
// de amostragem).  Ancorar a banda na resolucao mais GROSSA da serie, onde
// celmin=5, faz as mais finas so' ganharem celulas.

#ifndef VOF_AMR_3D_H
#define VOF_AMR_3D_H

#include "../src/hig-flow-kernel.h"

//! Converte uma banda FISICA em numero de celulas finas, com o piso de CINCO
//! vencendo e avisando.  `var` e' o nome da variavel de ambiente, so' para a
//! mensagem dizer qual pedido foi recusado.  Sem a variavel, devolve `celmin`.
int vof3_celmin_da_banda(const char *var, int celmin,
                         real lx, int nx, int niveis);

//! Escreve o `.amr` multinivel a partir do CAMPO de fracao do solver.
//! Semeia nas celulas de interface (0,001 < f < 0,999), REUNE entre ranks,
//! dilata e escreve.  So' o rank 0 escreve, e todos esperam.
void vof3_escreve_amr_criterio(higflow_solver *ns, real lx, real ly, real lz,
                               int nx, int ny, int nz,
                               int niveis, int cel_min, const char *caminho);

//! Criterio de ARRANQUE: ainda nao ha' dominio, logo a fracao sai da forma
//! ANALITICA -- a mesma funcao que inicializa o fracvol.  Mesma regra e mesmos
//! limiares do dinamico, para a malha do passo 0 nao ser de outra familia.
void vof3_escreve_amr_inicial(real lx, real ly, real lz, int nx, int ny, int nz,
                              int niveis, int cel_min,
                              real (*fracvol)(const Point, const Point, real),
                              const char *caminho);

#endif
