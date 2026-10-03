// *******************************************************************
//  CRITERIO DE REFINO ADAPTATIVO -- ver o .c para o porque' da separacao
// *******************************************************************
//
// Nada aqui sabe de front-tracking nem de VOF.  As funcoes recebem uma LISTA DE
// SEMENTES (centros de celula base cortadas pela interface) e devolvem a tabela
// de niveis e o arquivo .amr.  O que varia entre os metodos e' so' de onde as
// sementes vem.
//
// Sem dependencia de solver: nem PETSc, nem MPI.  O arreio compila este modulo
// direto.

#ifndef HIG_FLOW_AMR_CRITERIO_H
#define HIG_FLOW_AMR_CRITERIO_H

#include "types.h"
#include <stdio.h>
#include <stdlib.h>

// --- criterio de refino ----------------------------------------------------

//! Limiares da REGRA DAS CELULAS MINIMAS, derivados e nao digitados: o nivel
//! mais fino cobre ao menos `cel_min` celulas finas de cada lado da interface,
//! somada a margem que a discretizacao exige.  Saem DECRESCENTES -- thr[0] e' o
//! maior (nivel 1), thr[niveis-1] o menor (o nivel mais fino).  `thr` tem de
//! caber `niveis` reais.
void amr_limiares(real h_base, int niveis, int cel_min, real *thr);

//! Tabela de nivel por celula base, por DILATACAO das sementes.  O caminho de
//! producao: O(nsem x R^3) em vez do O(ncel x nsem) da busca exaustiva, que em
//! 3D daria ~1e9 operacoes por remalhamento numa malha de 250 mil celulas.
//!
//! `sem` tem 3*nsem reais (x,y,z de cada semente).  `tab` tem nx*ny*nz bytes.
//!
//! EXATA, nao aproximada: o nivel e' monotono decrescente na distancia, logo o
//! maximo sobre as sementes e' o nivel da semente mais proxima.  O contrato e'
//! igualdade bit a bit com `amr_tabela_exaustiva`, e o arreio o exige --
//! mesma disciplina da grade espacial.
void amr_tabela(const real *sem, long nsem, real hx, real hy, real hz,
                    int nx, int ny, int nz, const real *thr, int niveis,
                    signed char *tab);

//! A MESMA tabela por busca exaustiva.  Existe para o arreio; sem a segunda
//! implementacao a dilatacao viraria uma aproximacao silenciosa.  NAO usar em
//! producao.
void amr_tabela_exaustiva(const real *sem, long nsem,
                             real hx, real hy, real hz,
                             int nx, int ny, int nz,
                             const real *thr, int niveis, signed char *tab);

//! Escreve o .amr multinivel a partir da tabela de niveis.  Independente do
//! solver, de proposito: o arreio rele' o arquivo e confere as POSICOES dos
//! retalhos contra a tabela, nao so' as contagens.  Devolve 0 em caso de
//! sucesso.  Indice BASE 1 no arquivo.
int amr_escreve_arquivo(const char *caminho, real lx, real ly, real lz,
                            int nx, int ny, int nz, int niveis,
                            const signed char *tab);


#endif
