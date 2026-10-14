// CRITERIO DE REFINO ADAPTATIVO: limiares, tabela de niveis, nivel por ponto.
//
// A ARQUITETURA E' A DO HIGFLOW, e a separacao em tres camadas nao e' estetica.
// Nada aqui sabe DE ONDE vem o que se quer refinar: as funcoes recebem uma
// LISTA DE SEMENTES e devolvem a tabela de niveis.  O que e' especifico de cada
// problema e' so' quem produz as sementes -- a geometria de uma interface, o
// campo de fracao advectado, a vorticidade de um escoamento monofasico.
//
// A CONSEQUENCIA QUE IMPORTA: dois metodos que usem a MESMA construcao de
// tabela produzem malhas que so' diferem pelas SEMENTES.  Com dois criterios
// escritos em lugares diferentes, eles divergiriam com o tempo e a comparacao
// entre os metodos passaria a medir a divergencia dos criterios.
//
// E NAO DEPENDE DO SOLVER: nem PETSc, nem MPI, nem da classe de malha.  E' o
// que permite ao arreio compila-lo direto e exigir igualdade EXATA entre a
// tabela por dilatacao e a exaustiva.
//
// A SEGUNDA IMPLEMENTACAO NAO E' REDUNDANCIA.  A dilatacao e' `O(nsem . R^2)`
// contra `O(ncel . nsem)` da busca exaustiva, e e' por isso que ela e' o
// caminho de producao.  Mas ela so' e' EXATA porque o nivel e' monotono
// decrescente na distancia -- o maximo sobre as sementes e' o nivel da semente
// mais proxima.  Sem a exaustiva ao lado, essa propriedade seria uma crenca, e
// a dilatacao uma aproximacao silenciosa.

#ifndef T8FLOW_2D_CRITERIO2D_HXX
#define T8FLOW_2D_CRITERIO2D_HXX

#include <vector>

#include "malha2d.hxx"

namespace t8flow2d {

//! Os limiares da REGRA DAS CELULAS MINIMAS, DERIVADOS e nao digitados.
//
// O nivel mais fino tem de cobrir ao menos `cel_min` celulas finas de cada
// lado da semente, mais a margem que a decisao por celula base exige: a borda
// da banda e' incerta em `h_base/2` porque o teste e' no CENTRO da celula, e a
// celula se estende outro `h_base/2` alem dele.
//
// Saem DECRESCENTES: `thr[0]` e' o maior (nivel 1) e `thr[niveis-1]` o menor
// (o nivel mais fino).  Digitar os limiares a mao e' o jeito de a banda ficar
// com menos celulas do que se prometeu, e de ninguem perceber.
//
// `largura_fina` E' A LARGURA FISICA da banda mais fina, e o limiar e' o MAIOR
// entre ela e `cel_min` celulas finas.  As duas exigencias sao de naturezas
// diferentes: a contagem de celulas e' da DISCRETIZACAO e escala com `h`; a
// largura e' do PROBLEMA e nao escala.
//
// SO' A CONTAGEM FAZ O ESTUDO ANTICONVERGIR, e foi medido: na bolha de Hysing
// a banda valia 0,50 R em `h_min = 0,0125` e 0,25 R em 0,00625 -- refinar
// encolheu pela metade a regiao resolvida, e `v_max` passou de -6,5% a -10,5%
// enquanto a malha uniforme na mesma `h_min` dava +0,4%.
//
// O PADRAO E' ZERO e reduz exatamente ao comportamento anterior, de modo que
// quem nao passa largura continua bit a bit.  Quem tem escala fisica a
// declarar e' quem conhece o problema -- este arquivo nao sabe de que forma
// vem a semente.
std::vector<double> amr_limiares(double h_base, int niveis, int cel_min,
                                 double largura_fina = 0.0);

//! A tabela de nivel por celula base, por DILATACAO das sementes.
//
// O CAMINHO DE PRODUCAO.  Cada semente carimba um disco de raio dado pelo
// maior limiar, e dentro dele escreve o nivel que a distancia pede -- `O(nsem
// . R^2)` em vez de `O(ncel . nsem)`.  Numa malha de 65 mil celulas com mil
// sementes, a exaustiva daria 65 milhoes de distancias por remalhamento.
//
// `tab` sai com `nx*ny` valores, em ordem `i + nx*j`.
std::vector<signed char> amr_tabela(const std::vector<Ponto2> &sem,
                                    const Ponto2 &lo, double hx, double hy,
                                    int nx, int ny,
                                    const std::vector<double> &thr);

//! A MESMA tabela por busca exaustiva.
//
// EXISTE PARA O ARREIO, e nao para producao.  O contrato e' igualdade BIT A
// BIT com `amr_tabela`; sem ela a dilatacao seria uma aproximacao que ninguem
// conferiu.
std::vector<signed char> amr_tabela_exaustiva(const std::vector<Ponto2> &sem,
                                              const Ponto2 &lo, double hx,
                                              double hy, int nx, int ny,
                                              const std::vector<double> &thr);

//! O nivel que a tabela pede num ponto.  Fora da grade devolve 0.
int amr_nivel_em(const std::vector<signed char> &tab, const Ponto2 &lo,
                 double hx, double hy, int nx, int ny, const Ponto2 &x);

}  // namespace t8flow2d

#endif  // T8FLOW_2D_CRITERIO2D_HXX
