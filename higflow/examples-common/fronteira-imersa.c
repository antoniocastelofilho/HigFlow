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
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

// ---------------------------------------------------------------------------
// DUAS FORMULACOES, e a diferenca entre elas E' a medida.
//
//   DEFASADA     forca da velocidade do passo ANTERIOR, somada no campo de
//                fonte antes do preditor.  Impoe o nao escorregamento a O(dt).
//
//   UHLMANN      forca da velocidade PROVISORIA do passo corrente, corrigindo
//                o proprio u* depois do preditor.  Sem o termo O(dt).
//
// A escolha vem de HIGFLOW_FI_MODO (defasado | uhlmann).  A omissao e'
// `defasado` DE PROPOSITO: mudar o padrao alteraria em silencio o resultado dos
// exemplos que ja' existem.
// ---------------------------------------------------------------------------

static void _aplica_defasado(higflow_solver *ns, void *ctx)
{
    fi_corpo *corpo = (fi_corpo *) ctx;
    fi_interpola(corpo, ns->sfdu, ns->dpu);
    fi_forca_corpo_rigido(corpo, ns->par.dt);
    fi_zera_forca_passo(corpo);
    fi_acumula_forca_passo(corpo);
    fi_espalha(corpo, ns->sfdF, ns->dpFU);
}

// Cache do Uhlmann, preso ao psfdu corrente.  Criado sob demanda; destruido
// pelo aviso de remalha ANTES de os dominios morrerem -- depois nem o
// dp_destroy seria seguro.  MEDIDO sem isto: heap corrompido dois passos apos
// a reconstrucao, com a falha longe da causa.
static distributed_property *rascunho[DIM] = {NULL};

static void _remalha_avisa(higflow_solver *ns, void *ctx)
{
    (void) ns; (void) ctx;
    for (int dim = 0; dim < DIM; dim++) {
        if (rascunho[dim] != NULL) dp_destroy(rascunho[dim]);
        rascunho[dim] = NULL;
    }
}

static void _corrige_uhlmann(higflow_solver *ns, void *ctx)
{
    fi_corpo *corpo = (fi_corpo *) ctx;

    const char *sa = getenv("HIGFLOW_FI_ALFA");
    const real  alfa = (sa != NULL) ? atof(sa) : 1.0;
    const char *si = getenv("HIGFLOW_FI_ITER");
    int niter = (si != NULL) ? atoi(si) : 1;
    if (niter < 1) niter = 1;

    // CAMPO PROPRIO, ZERADO.  `fi_espalha` e' dono do campo que recebe: ele
    // acumula franja->dono, o que sobre `dpustar` somaria a velocidade de
    // franja nos donos.  Medido: com alfa = 0 -- forca anulada -- a corrida
    // ainda divergia para 1e96, o que prova que a corrupcao nao vinha da forca.
    // O rascunho e' cache PRESO AO DOMINIO: ver _remalha_avisa abaixo.
    for (int dim = 0; dim < DIM; dim++)
        if (rascunho[dim] == NULL)
            rascunho[dim] = psfd_create_property(ns->psfdu[dim]);

    // ITERAR, e o motivo NAO e' a defasagem temporal -- essa ja' foi resolvida
    // por calcular a forca de u*.  E' que espalhar e interpolar nao comutam
    // (H.A != I), entao uma aplicacao so' remove PARTE do escorregamento.
    // MEDIDO: uma aplicacao corta o residuo pela metade (4,47e-2 -> 2,20e-2).
    // E' a razao pela qual o multi-direct forcing da literatura usa 5 a 20.
    fi_zera_forca_passo(corpo);
    for (int it = 0; it < niter; it++) {
        fi_interpola(corpo, ns->sfdu, ns->dpustar);
        fi_forca_corpo_rigido(corpo, ns->par.dt);
        // SOMA SOBRE AS ITERACOES: cada uma acrescenta uma parcela, e a ultima
        // e' a menor.  Ler so' a ultima daria arrasto uma ordem pequeno demais.
        fi_acumula_forca_passo(corpo);

        for (int dim = 0; dim < DIM; dim++) {
            const int n = rascunho[dim]->pdata->total_count;
            for (int lid = 0; lid < n; lid++) dp_set_value(rascunho[dim], lid, 0.0);
        }
        fi_espalha_com_escala(corpo, ns->sfdu, rascunho, alfa * ns->par.dt);

        // u* += alfa * dt * F
        for (int dim = 0; dim < DIM; dim++) {
            const int n = rascunho[dim]->pdata->total_count;
            for (int lid = 0; lid < n; lid++)
                dp_add_value(ns->dpustar[dim], lid, dp_get_value(rascunho[dim], lid));
        }
        for (int dim = 0; dim < DIM; dim++) dp_sync(ns->dpustar[dim]);
    }

    // REINTERPOLA, e isto e' medicao, nao formulacao.  Sem isto o campo
    // `velocidade` do corpo fica com o valor de ANTES da ultima correcao, e o
    // medidor le' o residuo que ainda nao foi corrigido.  Foi o que aconteceu na
    // primeira comparacao, que deu Uhlmann igual a defasado.
    fi_interpola(corpo, ns->sfdu, ns->dpustar);
}

//! Instala o corpo no solver.  Uma linha no exemplo, como o `malha_t8_instala`.
extern "C" void fronteira_imersa_instala(higflow_solver *ns, fi_corpo *corpo)
{
    const char *modo = getenv("HIGFLOW_FI_MODO");
    // 'desligado' nao instala gancho nenhum.  Isto NAO e' o mesmo que alfa = 0:
    // com alfa = 0 a forca e' nula mas fi_espalha continua rodando, e foi
    // justamente por alfa = 0 ainda divergir que se achou a corrupcao de campo
    // pelo PetscSFReduce.  Para isolar o solver do corpo, o gancho tem que nao
    // existir.
    if (modo != NULL && strcmp(modo, "desligado") == 0) {
        print0f("=+=+=+= Fronteira imersa: DESLIGADA (sem gancho) =+=+=+=\n");
        return;
    }
    if (modo != NULL && strcmp(modo, "uhlmann") == 0) {
        higflow_set_fronteira_imersa_pos_preditor(ns, _corrige_uhlmann, corpo);
        higflow_set_remalha_avisa(ns, _remalha_avisa, NULL);
        const char *sa0 = getenv("HIGFLOW_FI_ALFA");
        const char *si0 = getenv("HIGFLOW_FI_ITER");
        print0f("=+=+=+= Fronteira imersa: UHLMANN (pos-preditor), alfa = %s, "
                "%s iteracao(oes) =+=+=+=\n",
                (sa0 != NULL) ? sa0 : "1.0", (si0 != NULL) ? si0 : "1");
    } else {
        if (modo != NULL && strcmp(modo, "defasado") != 0)
            fprintf(stderr, "HIGFLOW_FI_MODO=%s desconhecido; usando 'defasado'\n", modo);
        higflow_set_fronteira_imersa(ns, _aplica_defasado, corpo);
        print0f("=+=+=+= Fronteira imersa: forca DEFASADA (pre-preditor) =+=+=+=\n");
    }
}
