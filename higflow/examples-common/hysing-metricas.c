// Grandezas do benchmark de Hysing et al. (IJNMF 60:1259-1288, 2009), medidas
// EXATAMENTE como a secao 2.5 do artigo as define, e pela MESMA formula nos dois
// metodos -- que e' a condicao para comparar um com o outro e ambos com o
// publicado.
//
//   circularidade  c = Pa/Pb, o perimetro do circulo de AREA EQUIVALENTE sobre o
//                  perimetro da bolha = 2*sqrt(pi*A)/P.  ATENCAO: a forma
//                  4*pi*A/P^2, comum em codigo, e' o QUADRADO desta.  Eu usei a
//                  forma quadratica ate' aqui e comparei com o valor publicado
//                  sem converter: relatei 0,7987 contra 0,9013 quando o valor
//                  comparavel era sqrt(0,7987) = 0,8937.
//
//   centro de massa  yc = int_bolha y dx / int_bolha dx
//   velocidade de subida  Vc = int_bolha u_y dx / int_bolha dx
//
// Os dois integrais saem do campo de FRACAO, que o VOF advecta e o
// front-tracking escreve da geometria da frente -- entao a medida e' a mesma dos
// dois lados, e a diferenca que sobrar e' do metodo.

#include "hig-flow-kernel.h"
#include "hig-flow-eval.h"
#include "hig-mesh-snapshot.h"
#include <math.h>
#include <stdio.h>

//! Area, centro de massa e velocidade de subida, do campo de fracao.
extern "C" void hysing_medidas(higflow_solver *ns, real *area, real *yc, real *vc)
{
    sim_domain *sdm = psd_get_local_domain(ns->ed.mult.psdmult);
    mp_mapper  *mp  = sd_get_domain_mapper(sdm);
    sim_facet_domain *sfdv = psfd_get_local_domain(ns->psfdu[DIM-1]);

    real A = 0.0, sy = 0.0, su = 0.0;
    higcit_celliterator *it;
    for (it = sd_get_domain_celliterator(sdm); !higcit_isfinished(it);
         higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        int clid = mp_lookup(mp, hig_get_cid(c));
        if (clid < 0) continue;
        Point cc, d;
        hig_get_center(c, cc);
        hig_get_delta(c, d);
        real vol = 1.0;
        for (int k = 0; k < DIM; k++) vol *= d[k];
        real f = dp_get_value(ns->ed.mult.dpfracvol, clid);
        if (f <= 0.0) continue;
        real w = f * vol;
        A  += w;
        sy += w * cc[DIM-1];
        // u_y no centro da celula, interpolado das facetas
        su += w * compute_facet_value_at_point(sfdv, cc, cc, 1.0,
                                               ns->dpu[DIM-1], ns->stn);
    }
    higcit_destroy(it);

    real Ag, syg, sug;
    MPI_Allreduce(&A,  &Ag,  1, MPI_HIGREAL, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&sy, &syg, 1, MPI_HIGREAL, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&su, &sug, 1, MPI_HIGREAL, MPI_SUM, MPI_COMM_WORLD);

    *area = Ag;
    *yc   = (Ag > 0.0) ? syg / Ag : 0.0;
    *vc   = (Ag > 0.0) ? sug / Ag : 0.0;
}

//! Circularidade de Hysing a partir de area e perimetro ja' conhecidos.
//! Separada de proposito: o perimetro vem de onde cada metodo sabe medi-lo --
//! do poligono no front-tracking, da reconstrucao no VOF.
extern "C" real hysing_circularidade(real area, real perimetro)
{
    if (perimetro <= 0.0) return 0.0;
    return 2.0 * sqrt(M_PI * area) / perimetro;
}
