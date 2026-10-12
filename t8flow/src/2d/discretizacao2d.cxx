#include "discretizacao2d.hxx"

#include <stdexcept>

namespace t8flow2d {

namespace {

//! A familia de cada variante.  E' a tabela que impede o par incoerente.
Familia familia_de(Variante v)
{
  switch (v) {
    case Variante::FluxoBalanco:
    case Variante::MlsGauss:
      return Familia::VolumesFinitos;
    case Variante::CompactoMehrstellen:
    case Variante::CompactoLele:
    case Variante::MlsPontual:
      return Familia::DiferencasFinitas;
  }
  // Sem `default:` de proposito: assim o compilador avisa quando uma variante
  // nova entra no `enum` e nao entra nesta tabela.  O aviso na compilacao e' o
  // unico lugar onde esse esquecimento e' barato.
  throw std::invalid_argument("Discretizacao: variante sem familia declarada");
}

}  // namespace

Discretizacao::Discretizacao(Variante v, int ordem_mls)
    : variante_(v), familia_(familia_de(v)), ordem_mls_(0)
{
  const bool usa_mls = (v == Variante::MlsGauss || v == Variante::MlsPontual);
  if (usa_mls) {
    if (ordem_mls < 1 || ordem_mls > 4)
      throw std::invalid_argument(
          "Discretizacao: ordem do MLS " + std::to_string(ordem_mls) +
          " fora de [1, 4] -- 4 e' o teto da base propria do ajuste");
    ordem_mls_ = ordem_mls;
  }
  // Nas variantes que nao usam MLS a ordem fica ZERO, e isso e' informacao: um
  // `ordem_mls()` diferente de zero numa variante compacta significaria que
  // alguem acredita estar configurando algo que nao e' lido.
}

bool Discretizacao::incognita_e_media() const
{
  // NENHUMA VARIANTE DESTE SISTEMA USA INCOGNITA-MEDIA HOJE, e dizer isso e' mais
  // util que uma tabela com uma coluna constante.
  //
  // O `FluxoBalanco` e' de segunda ordem, onde os dois eixos coincidem: a media e
  // o valor central diferem por `O(h^2)`, que e' a ordem do proprio esquema, e
  // nenhuma medida deste sistema os separa.  Classifica-lo como pontual nao e'
  // arbitrario -- e' o que o oraculo mede, porque o erro e' comparado com
  // `exata(centroide)`.
  //
  // O `MlsGauss` integra a equacao e mantem a incognita pontual, e isso E'
  // observavel: ele e' exato para grau 3 contra o valor no centro.
  //
  // Volumes finitos com incognita-media -- o classico de alta ordem -- pede um
  // passo de RECONSTRUCAO media para pontual que este sistema nao tem.  Quando
  // tiver, a variante entra aqui e esta funcao deixa de ser constante.
  return false;
}

bool Discretizacao::trata_malha_adaptada() const
{
  // O ESTENCIL COMPACTO DE NOVE PONTOS PRESSUPOE ESPACAMENTO UNIFORME, e isso
  // nao e' uma limitacao da implementacao: os coeficientes `{-20, 4, 1}/6h^2`
  // saem de uma expansao de Taylor em torno de um ponto com os oito vizinhos a
  // distancias `h` e `h*sqrt(2)`.  Numa transicao 2:1 as distancias mudam e os
  // coeficientes deixam de representar o laplaciano.
  //
  // Aplicar assim mesmo nao falharia: daria um operador consistente com OUTRA
  // equacao, e a solucao convergiria -- para a coisa errada.
  return variante_ != Variante::CompactoMehrstellen &&
         variante_ != Variante::CompactoLele;
}

bool Discretizacao::estavel_em_malha_adaptada() const
{
  // As duas compactas nem chegam la': recusam a malha antes.
  if (!trata_malha_adaptada()) return false;

  // O BALANCO DE FLUXOS E' ESTAVEL POR CONSTRUCAO.  Cada faceta contribui
  // `+coef` ao vizinho e `-coef` a' diagonal, com `coef = |face|/dist > 0`: a
  // linha e' diagonalmente dominante e os sinais sao os certos, em qualquer
  // malha.  O `MlsGauss` herda isso, porque tambem e' balanco de fluxos -- o que
  // muda nele e' como a derivada na faceta e' calculada, nao a estrutura.
  if (familia_ == Familia::VolumesFinitos) return true;

  // O LAPLACIANO PONTUAL NAO GARANTE NADA.  Ele monta a linha que aproxima o
  // laplaciano, e a estrutura resultante depende do suporte.  Em malha uniforme
  // o suporte e' regular e sai bem; em adaptada, medido com ordem 3: 20 mil
  // iteracoes sem convergir.  Com ordem 4 convergiu, mas em 1,4 mil iteracoes
  // contra 59 da ordem 3 na uniforme -- ou seja, resolve, e caro.
  //
  // A ordem 4 nao e' declarada estavel aqui com base numa corrida que convergiu:
  // convergir uma vez nao e' estabilidade.  Enquanto nao houver medida que
  // sustente, a resposta e' NAO, e quem quiser usar assim mesmo o faz sabendo.
  return false;
}

int Discretizacao::grau_exato(bool malha_adaptada) const
{
  switch (variante_) {
    case Variante::FluxoBalanco:
      // Grau 2 na uniforme; na adaptada a quadratura da faceta (ponto medio) e a
      // distancia entre centros (`3h/4`, nao `h/2`) tiram a exatidao, e sobra o
      // campo linear.
      return malha_adaptada ? 1 : 2;

    case Variante::MlsGauss:
      // O GRAU ACOMPANHA A ORDEM DO AJUSTE, nos dois tipos de malha.  E' a
      // propriedade que distingue esta variante das outras, e a razao de existir.
      return ordem_mls_;

    case Variante::MlsPontual:
      // O mesmo: medido exato ate' grau 3 com ordem 3, e ate' grau 4 com ordem 4.
      //
      // UMA SUTILEZA QUE A TABELA NAO CAPTURA, e que vale saber: num suporte
      // SIMETRICO -- que e' o que o centro de uma celula ve' em malha uniforme --
      // os monomios IMPARES nao contribuem para o laplaciano, de modo que a
      // ordem 2 ja' alcanca o grau 3.  Nao esta' codificado aqui porque a ordem 2
      // nao foi medida NO POISSON, so' no ajuste isolado, e tabela de grau de
      // exatidao com entrada nao medida e' exatamente o que esta lista nao e'.
      (void)malha_adaptada;
      return ordem_mls_;

    case Variante::CompactoMehrstellen:
    case Variante::CompactoLele:
      // Grau 3, e o limite e' o FECHAMENTO DE PAREDE (o fantasma de quatro
      // pontos), nao o estencil interior -- que alcanca grau 4.  A ordem de
      // convergencia nao sofre com isso, porque o erro do fechamento e' local e
      // tambem `O(h^4)`; o que se perde e' a exatidao nos graus 4 e 5.
      return 3;
  }
  throw std::invalid_argument("Discretizacao::grau_exato: variante desconhecida");
}

int Discretizacao::ordem() const
{
  switch (variante_) {
    case Variante::FluxoBalanco:        return 2;
    case Variante::CompactoMehrstellen:
    case Variante::CompactoLele:        return 4;
    case Variante::MlsGauss:
    case Variante::MlsPontual:          return ordem_mls_;
  }
  throw std::invalid_argument("Discretizacao::ordem: variante desconhecida");
}

std::string Discretizacao::nome() const
{
  switch (variante_) {
    case Variante::FluxoBalanco:        return "vf-fluxo2";
    case Variante::MlsGauss:            return "vf-mlsgauss" + std::to_string(ordem_mls_);
    case Variante::CompactoMehrstellen: return "df-mehrstellen4";
    case Variante::CompactoLele:        return "df-lele4";
    case Variante::MlsPontual:          return "df-mlspontual" + std::to_string(ordem_mls_);
  }
  throw std::invalid_argument("Discretizacao::nome: variante desconhecida");
}

std::string Discretizacao::nomes_validos()
{
  // O PREFIXO `vf-` OU `df-` E' PARTE DO NOME, e nao enfeite: quem escolhe na
  // linha de comando precisa ver a familia na propria escolha.  "mlsgauss" e
  // "mlspontual" sao o mesmo ajuste em familias diferentes, e sem o prefixo a
  // diferenca -- que e' o que a incognita significa -- ficaria invisivel.
  // AS ORDENS APARECEM NA LISTA, uma a uma, porque "volumes finitos de quarta
  // ordem" e' uma ESCOLHA e nao uma configuracao: quem le' esta lista tem de ver
  // que ela existe sem precisar descobrir um parametro separado.
  //
  // `df-mlspontual` sai da lista com uma marca, e nao e' por capricho: ele nao e'
  // estavel em malha adaptada (clausula W6), e uma lista de escolhas que nao
  // distingue o verificado do experimental convida ao erro mais caro.
  return "vf-fluxo2 (2a ordem), vf-mlsgauss2, vf-mlsgauss3, vf-mlsgauss4 "
         "(volumes finitos de ordem alta), df-mehrstellen4, df-lele4 "
         "(so' malha uniforme), df-mlspontual2..4 (EXPERIMENTAL: instavel em "
         "malha adaptada)";
}

namespace {

//! Separa o sufixo de ordem do nome: "vf-mlsgauss4" -> ("vf-mlsgauss", 4).
//
// Devolve a ordem padrao quando nao ha' sufixo, para que `vf-mlsgauss` continue
// valendo -- tirar a forma curta quebraria quem ja' a usa sem ganhar nada.
std::string separa_ordem(const std::string &nome, int padrao, int *ordem)
{
  *ordem = padrao;
  if (nome.empty()) return nome;
  const char ultimo = nome[nome.size() - 1];
  if (ultimo >= '1' && ultimo <= '9') {
    const std::string base = nome.substr(0, nome.size() - 1);
    // "vf-fluxo2" tem digito no fim e NAO e' sufixo de ordem -- o `2` e' parte do
    // nome.  A base tem de terminar em letra para o digito contar como sufixo.
    if (!base.empty() && base != "vf-fluxo") {
      *ordem = ultimo - '0';
      return base;
    }
  }
  return nome;
}

}  // namespace

Discretizacao Discretizacao::do_nome(const std::string &nome, int ordem_mls)
{
  int ord = ordem_mls;
  const std::string base = separa_ordem(nome, ordem_mls, &ord);

  if (nome == "vf-fluxo2")      return Discretizacao(Variante::FluxoBalanco);
  if (base == "vf-mlsgauss")    return Discretizacao(Variante::MlsGauss, ord);
  if (base == "df-mehrstellen") return Discretizacao(Variante::CompactoMehrstellen);
  if (base == "df-lele")        return Discretizacao(Variante::CompactoLele);
  if (base == "df-mlspontual")  return Discretizacao(Variante::MlsPontual, ord);

  // A LISTA VAI NA MENSAGEM.  Um erro que diz apenas "nome invalido" obriga quem
  // o recebe a ler o codigo; com a lista, a correcao e' imediata.
  throw std::invalid_argument("Discretizacao: \"" + nome +
                              "\" nao existe.  Validos: " + nomes_validos());
}

}  // namespace t8flow2d
