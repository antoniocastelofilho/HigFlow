// VOF: O NUCLEO GEOMETRICO DA RECONSTRUCAO LINEAR POR PEDACOS (PLIC).
//
// Duas funcoes, e uma e' a inversa da outra:
//
//   `plic_area`       dada a normal e o deslocamento da reta, a AREA que ela
//                     corta da celula;
//   `plic_distancia`  dada a normal e a area desejada, o DESLOCAMENTO que a
//                     corta exatamente.
//
// A segunda e' a que o metodo usa a cada passo -- a fracao volumetrica diz
// QUANTO, a normal diz PARA ONDE, e falta onde por a reta.
//
// A CONVENCAO E' A DO HIGFLOW, e vale escreve-la porque la' ela e' implicita e
// e' facil errar o lado:
//
//   * a celula e' `[-h0/2, h0/2] x [-h1/2, h1/2]`, em coordenadas RELATIVAS AO
//     CENTRO.  `d` e' medido do centro, ao longo da normal unitaria;
//   * a area devolvida e' a da regiao `{x : n.x >= d}`, isto e' o lado para o
//     qual a normal APONTA.
//
// Conferido a mao contra `area_left_line_origin_center` de
// `hig-flow-vof-plic.c` em tres casos: eixo-alinhado, meia celula na diagonal,
// e triangulo de canto de area 0,125.
//
// ===========================================================================
// POR QUE HA' DUAS IMPLEMENTACOES DE CADA UMA
// ===========================================================================
//
// `plic_area` e `plic_distancia` sao FORMA FECHADA: reduzem o problema a um
// quadrado unitario canonico e invertem a quadratica por pedacos.  Exatas e sem
// iteracao -- e, sozinhas, INVERIFICADAS: a analise de casos tem seis ramos, e
// um ramo errado num regime estreito (normal quase alinhada a um eixo, celula
// muito achatada) da' numero plausivel.
//
// `plic_area_por_recorte` e `plic_distancia_por_bisseccao` existem PARA O
// ARREIO.  A primeira recorta o retangulo da celula contra o semiplano por
// Sutherland-Hodgman e integra por formula do sapateiro; a segunda bisseca
// sobre ela.  NAO COMPARTILHAM ALGEBRA NENHUMA com as de cima -- e' por isso que
// a concordancia entre as duas vias afirma alguma coisa.
//
// E' o terceiro par deste tipo no sistema, e os dois anteriores pagaram:
// `amr_tabela` contra `amr_tabela_exaustiva`, e o operador de difusao pelos
// dois caminhos.  O padrao e' sempre o mesmo -- a gemea nao e' redundancia, e'
// o que torna a propriedade VERIFICADA em vez de acreditada.
//
// HA' TAMBEM UMA PROPRIEDADE QUE NAO PRECISA DE GEMEA: a ida e volta,
// `plic_area(n, h, plic_distancia(n, h, A)) == A`.  Ela e' mais forte do que
// parece.  `plic_area` e' ESTRITAMENTE DECRESCENTE em `d` no intervalo ativo,
// logo a inversa e' UNICA -- e dai' sai o argumento de correspondencia com o
// HiGFlow: se `plic_area` reproduz a algebra deles e a inversa e' exata, a
// minha `plic_distancia` TEM de dar o que a deles da', sem eu precisar linkar
// o codigo deles para comparar.

#ifndef T8FLOW_2D_VOF2D_HXX
#define T8FLOW_2D_VOF2D_HXX

#include <functional>

#include "campos2d.hxx"
#include "malha2d.hxx"
#include "mls2d.hxx"

namespace t8flow2d {

//! A area de `{x : n.x >= d}` dentro da celula de lados `h` centrada na
//! origem.  `n` NAO precisa ser unitaria -- e' normalizada aqui --, e `d` e'
//! sempre em unidades de comprimento ao longo da normal UNITARIA.
//!
//! Lanca `std::invalid_argument` se `|n|` for nula ou algum lado nao for
//! positivo.  Normal nula nao tem lado para apontar, e aceitar calado daria
//! area arbitraria numa celula de interface.
double plic_area(const Ponto2 &n, const Ponto2 &h, double d);

//! O `d` tal que `plic_area(n, h, d) == area`.
//!
//! `area` tem de estar em `[0, h0*h1]`.  FORA DISSO LANCA, e isso e'
//! deliberado: a fracao volumetrica sai de `[0,1]` em advecao direcional
//! partida -- e' conhecido e e' tratado onde se trata, por ajuste CONTADO, nao
//! aqui por saturacao silenciosa.  Saturar aqui esconderia a unica evidencia
//! de que o esquema de transporte passou do limite.
double plic_distancia(const Ponto2 &n, const Ponto2 &h, double area);

// ---------------------------------------------------------------------------
// AS GEMEAS.  Ver o cabecalho deste arquivo: existem para o arreio exigir que
// as duas vias concordem, e nao para producao.
// ---------------------------------------------------------------------------

//! A mesma area, por RECORTE do retangulo contra o semiplano e formula do
//! sapateiro.  `O(1)` com constante maior, e sem nenhuma algebra em comum.
double plic_area_por_recorte(const Ponto2 &n, const Ponto2 &h, double d);

//! O mesmo `d`, por BISSECCAO sobre `plic_area_por_recorte`.
//!
//! Converge porque a area e' monotona em `d`; o intervalo inicial e' o
//! suporte exato, `[-(|n0|h0+|n1|h1)/2, +(...)/2]`, onde a area vale o volume
//! inteiro e zero.
double plic_distancia_por_bisseccao(const Ponto2 &n, const Ponto2 &h,
                                    double area);

// ===========================================================================
// A FRACAO VOLUMETRICA INICIAL
// ===========================================================================
//
// E' O MESMO CAMINHO DO HIGFLOW, e a correspondencia aqui nao e' preciosismo:
// se o volume inicial diferir, TODAS as grandezas do benchmark partem de outro
// lugar -- a area de referencia, a circularidade, o centro de massa -- e a
// comparacao entre os dois sistemas mede a diferenca de inicializacao
// misturada com a diferenca de metodo.
//
// O CAMINHO, de `get_fracvolN` em `examples-common/vof-geometry-2d.c` e do
// envoltorio em `example2d_VOF/ns-example-2d.c`:
//
//   1. avalia a funcao de forma nos QUATRO CANTOS da celula.  `forma > 0` e'
//      dentro;
//   2. zero cantos dentro -> fracao 0; quatro -> fracao 1.  Sem conta;
//   3. qualquer outro caso -> SUBDIVIDE a celula em `n x n` e soma a area de
//      cada sub-celula pelo caso analitico (um canto, tres cantos, ou dois a
//      dois), com a interface LINEARIZADA dentro de cada sub-celula.
//
// O PASSO 3 E' O QUE FAZ O VOLUME FICAR CERTO, e e' facil nao notar que esta'
// la': linearizar a interface dentro da celula inteira erra a area de uma
// curva de curvatura `1/R` em `O(h^3/R)` por celula cortada.  Com `n = 16` --
// o valor do HiGFlow -- o erro cai por `n^2 = 256`, e a fracao inicial passa a
// ser exata para o que interessa.
//
// O CASO AMBIGUO E' RECUSADO, nao resolvido por convencao.  Com dois cantos
// dentro na DIAGONAL, a celula pode conter duas quinas separadas ou uma
// garganta, e os quatro cantos nao distinguem.  O HiGFlow aborta dizendo
// "malha pouco refinada"; aqui lanca, com a mesma leitura: nao e' defeito do
// codigo, e' a malha nao resolvendo a forma.

//! A AREA (nao a fracao) da parte da celula onde `forma > 0`, por um unico
//! quadrado -- a interface linearizada entre os quatro cantos.
//!
//! EXATA quando `forma` e' afim; e' o caso em que ela tem de coincidir com
//! `plic_area`, que vem por outro caminho.
double area_da_forma(const std::function<double(const Ponto2 &)> &forma,
                     const Ponto2 &centro, const Ponto2 &h);

//! A FRACAO em `[0,1]`, com a subdivisao `n x n` nas celulas cortadas.
//!
//! `n = 16` e' o valor do HiGFlow.  `n = 1` desliga a subdivisao e deixa a
//! interface linearizada na celula inteira -- serve de CONTRASTE, para o erro
//! que a subdivisao remove ficar medido em vez de alegado.
double fracao_inicial(const std::function<double(const Ponto2 &)> &forma,
                      const Ponto2 &centro, const Ponto2 &h, int n = 16);

//! O campo inteiro.  Numa malha adaptada cada folha usa o PROPRIO lado -- a
//! subdivisao e' relativa a' celula, nao a um `h` global.
void escreve_fracao_inicial(const Malha &m,
                            const std::function<double(const Ponto2 &)> &forma,
                            CampoEscalar *frac, int n = 16);

}  // namespace t8flow2d

#endif  // T8FLOW_2D_VOF2D_HXX
