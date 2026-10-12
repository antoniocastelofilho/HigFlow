// O ACOPLAMENTO DA FRENTE COM A MALHA: fase B2.
//
// EM ARQUIVO SEPARADO DE `frente2d`, de proposito e pela mesma razao que o
// HiGFlow o separa: a frente (B1) e' advecao cinematica e nao sabe que existe
// malha, campo ou solver.  Os oraculos dela sao de forma fechada e valem
// sozinhos.  Misturar as duas faria B1 deixar de ser verificavel sem o resto.
//
// SAO TRES COISAS QUE ATRAVESSAM A FRONTEIRA entre a frente e a malha:
//
//   A FRACAO VOLUMETRICA, escrita da GEOMETRIA da frente -- e e' por ela que a
//   comparacao com o VOF e' possivel.  Os dois metodos sao medidos pelas
//   MESMAS integrais sobre o MESMO campo: o VOF o advecta, o front-tracking o
//   escreve.  Qualquer diferenca que sobreviva pertence ao metodo, e nao a'
//   medicao;
//
//   AS PROPRIEDADES -- densidade e viscosidade -- interpoladas da fracao;
//
//   A TENSAO SUPERFICIAL, distribuida dos marcadores para as facetas.
//
// A FRACAO E' EXATA, por recorte de poligono, e nao aproximada por amostragem.
// Isso nao e' refinamento gratuito: com recorte exato, a SOMA das fracoes
// vezes os volumes da' a area do poligono em precisao de maquina -- um oraculo
// que uma aproximacao por amostragem nao permitiria afirmar.  E' tambem o que
// torna o campo utilizavel como medida: as integrais do benchmark saem dele.

#ifndef T8FLOW_2D_ACOPLAMENTO2D_HXX
#define T8FLOW_2D_ACOPLAMENTO2D_HXX

#include <vector>

#include "campos2d.hxx"
#include "frente2d.hxx"
#include "malha2d.hxx"

namespace t8flow2d {

//! A AREA DA INTERSECAO entre um poligono e uma caixa alinhada aos eixos.
//
// Por recorte de Sutherland-Hodgman: a caixa e' convexa, e contra uma regiao
// de recorte convexa o algoritmo da' a area da intersecao EXATAMENTE, mesmo
// para um poligono nao convexo -- as arestas degeneradas que ele pode produzir
// ao longo da borda nao contribuem para o cadarco.
//
// UM ALGORITMO SO', sem casos.  Celula inteiramente dentro devolve a area da
// caixa; inteiramente fora devolve zero.  Tratar "dentro" e "fora" como casos
// separados exigiria um teste de ponto-em-poligono cujo desempate na fronteira
// teria de concordar com o do recorte -- e dois desempates sao dois lugares
// para divergirem.
double area_recortada(const std::vector<Ponto2> &poligono, const Ponto2 &lo,
                      const Ponto2 &hi);

//! Escreve a FRACAO VOLUMETRICA de cada celula a partir da frente.
//
// A soma de `fracao * volume` sobre as celulas da' a area do poligono, em
// precisao de maquina.  E' o oraculo de particao da unidade, e ele so' e'
// afirmavel porque o recorte e' exato.
void escreve_fracao(const Frente &f, const Malha &m, CampoEscalar *c);

//! Densidade e viscosidade por interpolacao LINEAR da fracao.
//
// `prop = fora + (dentro - fora) * c`.
//
// A ALTERNATIVA PARA A VISCOSIDADE -- media HARMONICA -- nao esta' aqui, e a
// ausencia e' escolha declarada: ela e' defensavel (preserva a tensao
// cisalhante atraves de uma interface plana) e DIFERENTE, e o HiGFlow usa a
// linear.  Trocar o modelo de propriedade junto com tudo o mais faria qualquer
// diferenca medida contra ele ficar sem dono.
void interpola_propriedade(const Malha &m, const CampoEscalar &c,
                           double dentro, double fora, CampoEscalar *prop);

//! A FORCA DE TENSAO SUPERFICIAL nas facetas, distribuida dos marcadores.
//
//     f(x) = soma_k  sigma . kappa_k . n_k . ds_k . D_h(x - X_k)
//
// `D_h` e' o nucleo de Peskin de quatro pontos, que tem duas propriedades que
// nao sao decorativas: soma UM sobre qualquer grade deslocada (logo a forca
// total distribuida e' a forca total dos marcadores, exatamente), e e'
// continuo com derivada continua (logo a forca nao salta quando um marcador
// cruza uma faceta, o que produziria ruido no campo de velocidade).
//
// O NUCLEO USA O `h` DA CELULA QUE CONTEM O MARCADOR.  Numa malha adaptada nao
// ha' um `h` -- e distribuir com o `h` errado espalha a forca por uma largura
// errada, o que nao falha: muda a curvatura efetiva e desloca o salto de
// pressao sem que nada acuse.
void forca_tensao(const Frente &f, const Malha &m, const DominioFacetas &fd,
                  double sigma, CampoVetorial *forca);

//! A soma de `forca * medida` sobre as facetas, por direcao.
//
// TEM DE DAR A SOMA DE `sigma . kappa . n . ds` SOBRE OS MARCADORES, que para
// uma curva FECHADA e' ZERO -- a tensao superficial nao exerce forca liquida.
// Um valor nao nulo significa nucleo que nao soma um, ou marcador cuja forca
// caiu fora do dominio.
Ponto2 forca_total(const DominioFacetas &fd, const CampoVetorial &forca);

}  // namespace t8flow2d

#endif  // T8FLOW_2D_ACOPLAMENTO2D_HXX
