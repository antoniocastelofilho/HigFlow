// TESTE DE LINK, nao de execucao.  A pergunta e' binaria: o caminho
// multifasico VISCOELASTICO resolve todos os simbolos em DIM=3?  Se resolver, o
// que falta e' caso e verificacao, nao codigo.  Se nao resolver, os simbolos que
// faltam nomeiam exatamente o que nao foi portado para tres dimensoes.
//
// A chamada esta' sob uma condicao que nunca ocorre: o objetivo e' obrigar o
// ligador a resolver o simbolo, e nao executar nada.
#include "hig-flow-kernel.h"
#include "hig-flow-step-multiphase.h"
#include "hig-flow-step-multiphase-viscoelastic.h"
#include <stdio.h>

int main(int argc, char **argv) {
    if (argc > 9999) {
        higflow_solver *ns = NULL;
        higflow_solver_step_multiphase(ns);
        higflow_solver_step_multiphase_viscoelastic(ns);
    }
    printf("link 3D do multifasico viscoelastico: OK\n");
    return 0;
}
