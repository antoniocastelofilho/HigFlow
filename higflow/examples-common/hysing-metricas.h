// Grandezas do benchmark de Hysing et al. (2009).  Ver o .c para as definicoes,
// que sao as da secao 2.5 do artigo -- em particular a circularidade, que e'
// Pa/Pb e NAO 4*pi*A/P^2 (esta e' o quadrado daquela).
#ifndef HYSING_METRICAS_H
#define HYSING_METRICAS_H
#include "../src/hig-flow-kernel.h"
#ifdef __cplusplus
extern "C" {
#endif
void hysing_medidas(higflow_solver *ns, real *area, real *yc, real *vc);
real hysing_circularidade(real area, real perimetro);
#ifdef __cplusplus
}
#endif
#endif
