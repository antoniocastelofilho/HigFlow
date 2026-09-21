// Adaptador entre o solver e o modulo da fronteira imersa -- o mesmo papel que
// malha-t8.cxx faz para a fonte de malha, e pelo mesmo motivo.
//
// O `hig-flow-step.c` e' ligado por TODOS os exemplos.  Se ele chamasse
// `fi_interpola` direto, os doze Makefiles teriam de ligar o modulo da
// fronteira imersa para um recurso que um exemplo usa.  Com o gancho, quem nao
// instala corpo nao referencia simbolo nenhum, e quem instala compila ESTE
// arquivo junto.

#include "hig-flow-kernel.h"
#include "hig-flow-fronteira-imersa.h"

// A ordem e' a do projeto, secao 3: interpolar a velocidade do passo ANTERIOR
// (o preditor ainda nao rodou), tirar dela a forca direta, e espalhar no campo
// de forca por faceta que a equacao ja' le'.
static void _aplica(higflow_solver *ns, void *ctx)
{
    fi_corpo *corpo = (fi_corpo *) ctx;
    fi_interpola(corpo, ns->sfdu, ns->dpu);
    fi_forca_corpo_rigido(corpo, ns->par.dt);
    fi_espalha(corpo, ns->sfdF, ns->dpFU);
}

//! Instala o corpo no solver.  Uma linha no exemplo, como o `malha_t8_instala`.
extern "C" void fronteira_imersa_instala(higflow_solver *ns, fi_corpo *corpo)
{
    higflow_set_fronteira_imersa(ns, _aplica, corpo);
}
