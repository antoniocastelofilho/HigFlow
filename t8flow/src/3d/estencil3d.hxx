// O ESTENCIL COM A ORDEM DESEJADA, como na HiGTree.
//
// La' a interface e' `sd_set_interpolator_order` + `sd_get_stencil(sd, centro, x,
// alpha, stn)`, devolvendo um `sim_stencil` com `ids[]`, `vals[]` (os pesos) e um
// `rhs` SEPARADO para o que nao se associa a nenhuma celula -- o valor de
// Dirichlet.  Aqui a forma e' a mesma, e a separacao do `rhs` tambem: guardar o
// valor de contorno junto com os pesos obriga o chamador a saber qual id e' "a
// parede", e e' assim que o fechamento se perde.
//
// O QUE ESTE MODULO ACRESCENTA ao que ja' existia -- o `Mls3` sabia calcular pesos
// de um suporte DADO, e o `LoteEstenceis` sabia faze-lo em lote -- e' a ESCOLHA DO
// SUPORTE:
//
//   1. acha a celula que contem o ponto (busca geometrica da malha);
//   2. cresce por ANEIS de vizinhanca ate' ter pontos para a ordem pedida;
//   3. PARA AO FECHAR O ANEL, nunca no meio dele.
//
// O PASSO 3 NAO E' DETALHE.  Parar no meio de um anel faz o suporte depender da
// ORDEM em que os vizinhos foram visitados, e a resposta passa a mudar com o
// percurso -- que e' exatamente como a clausula de desempate da HiGTree passava
// por sorte.  E crescer SEM PARADA e' pior: busca em largura sem limite devolve a
// malha inteira (medido na HiGTree: 64 pontos numa malha de 64), o ajuste
// reproduz polinomio mesmo assim, e o teste passa enquanto a topologia deixou de
// importar.  Daqui sai um suporte LOCAL, e a clausula F3 afirma isso.
//
// O FECHAMENTO DE CONTORNO.  Quando o anel alcanca o contorno do dominio, os
// pontos que faltam vem das FACETAS de contorno, nao de celulas: o valor ali e'
// conhecido (Dirichlet), entao o peso daquele ponto multiplica o valor e vai para
// o `rhs`.  O suporte tem dois tipos de ponto, e a diferenca nao esta' na
// geometria -- esta' em quem responde pelo valor.
//
// LIMITE DE HOJE, declarado: so' Dirichlet.  Neumann precisa de derivada na
// faceta, que muda a LINHA do ajuste e nao so' o lado direito -- e' trabalho
// proprio, nao um `else` aqui.
#ifndef T8FLOW_ESTENCIL2D_HXX
#define T8FLOW_ESTENCIL2D_HXX

#include <functional>
#include <vector>

#include "contorno3d.hxx"
#include "malha3d.hxx"
#include "mls3d.hxx"

namespace t8flow3d {

//! O ESTENCIL: os ids, os pesos, e o que nao se associa a celula nenhuma.
//
// `valor_no_ponto = soma(pesos[k] * u[ids[k]]) + rhs`
//
// O `rhs` existe porque o fechamento de contorno produz termos CONHECIDOS -- o
// valor de Dirichlet na parede nao e' incognita de ninguem --, e deixa-los na
// matriz exigiria uma coluna que nao existe.
//
// NO LADO 2D ESTA ESTRUTURA VIVE NO CABECALHO DO LOTE, que e' a otimizacao em
// GPU.  E' o lugar errado: quem usa o estencil nao deveria precisar incluir o
// acelerador, e a prova e' que o 3D -- que nao tem lote -- precisou dela.  Aqui
// ela fica onde pertence.
struct Estencil {
  std::vector<long> ids;
  std::vector<double> pesos;
  double rhs = 0.0;
};

//! O dado de uma condicao, como funcao do ponto -- nao uma constante.
//
// Valor de contorno CONSTANTE e' o caso em que quase todo defeito de fechamento
// fica invisivel: com parede de nao-deslizamento, fechar contra a parede errada
// da' o mesmo numero.  Foi assim que nove defeitos de fechamento sobreviveram
// cinco anos na HiGTree.  Daqui para a frente, todo caso de teste usa dado que
// VARIA nas duas direcoes.
using ValorDeContorno = DadoDeContorno;

//! COMO O SUPORTE E' ESCOLHIDO.
//
// A ESCOLHA NAO E' DE GOSTO, e a medida decidiu: com `Aneis`, 504 de 512 ajustes
// de ordem 3 REBAIXARAM a ordem numa malha uniforme de nivel 3 -- 98,4%.
//
// A razao e' a DISTRIBUICAO, nao a contagem.  O crescimento por faces entrega um
// suporte em CRUZ (media de 24,4 pontos para 20 monomios, numero de sobra), mas
// quase todos alinhados com os eixos; os monomios MISTOS -- `xy`, `xz`, `yz`,
// `xyz` -- ficam mal determinados nele.  Em duas dimensoes isso nao aparece
// porque ha' menos termos mistos: a cruz de dois aneis da' 13 pontos para 10
// monomios e sustenta.
enum class Suporte3 {
  //! Anéis de vizinhos de FACE, como no lado 2D.  Mantido porque e' o que casa
  //! com a franja (que tambem e' por face) e basta para a ordem 2.
  Aneis,

  //! CAIXA geometrica: os pontos em `c + (i,j,k)*h` para `|i|,|j|,|k| <= r`.
  //
  // Da' 27 pontos bem distribuidos num cubo 3x3x3, que e' a configuracao
  // classica para ajuste tridimensional -- e condiciona a matriz normal muito
  // melhor que a cruz de mesmo tamanho.
  //
  // So' e' viavel porque `Malha3::id_do_ponto` e' O(niveis): com a varredura
  // linear que o 2D tinha antes, 27 buscas por ponto de consulta tornariam a
  // montagem inviavel.
  Caixa,
};

class ConstrutorEstencil {
 public:
  //! `ordem` e' o grau reproduzido exatamente.  `max_pontos` e' o teto do
  //! suporte; `max_aneis` limita o crescimento -- sem ele, um ponto numa regiao
  //! grossa cercada de contorno cresceria indefinidamente.
  ConstrutorEstencil(const Malha3 &m, int ordem, int max_pontos, int max_aneis = 6);

  //! Instala as condicoes de fronteira POR REGIAO -- parede norte diferente da
  //! sul, que e' o que o Poisson com solucoes manufaturadas precisa.
  //
  // Sem condicao instalada o contorno NAO contribui, e o suporte perto da parede
  // fica com menos pontos -- o que faz `em()` LANCAR em vez de devolver um ajuste
  // de ordem menor sem avisar.  E com condicoes instaladas, uma faceta de contorno
  // que nenhuma regiao cobrir faz `em()` LANCAR tambem: contorno sem condicao e'
  // pior que contorno ausente, porque a conta roda.
  void contorno(ContornoPorRegiao r) { regioes_ = std::move(r); }

  //! Atalho: UMA condicao para todo o contorno.  E' o caso simples, e ele continua
  //! sendo uma regiao -- a que cobre tudo --, nao um caminho separado: dois
  //! caminhos divergem.
  void contorno(Contorno c)
  {
    ContornoPorRegiao r;
    r.padrao(std::move(c));
    regioes_ = std::move(r);
  }

  //! Atalho do atalho, para o caso mais comum de todos.
  void dirichlet(ValorDeContorno v) { contorno(Contorno::dirichlet(std::move(v))); }

  //! As regioes instaladas, para diagnostico e para o teste poder afirmar QUAL
  //! regiao respondeu em cada parede.
  const ContornoPorRegiao &regioes() const { return regioes_; }

  int ordem() const { return ordem_; }
  int max_pontos() const { return max_pontos_; }

  //! O estencil em `x`: ids e pesos das celulas, e `rhs` com a contribuicao do
  //! contorno.
  //
  // Lanca `std::out_of_range` se `x` estiver fora do dominio, e
  // `std::runtime_error` se o suporte nao alcancar o minimo da ordem dentro de
  // `max_aneis` -- recusar e' deliberado: um ajuste de ordem menor que a pedida
  // devolve numeros plausiveis e errados.
  Estencil em(const Ponto3 &x) const;

  //! O estencil da DERIVADA `du/dx_dir` em `x`: ids, pesos e `rhs`.
  //
  // Mesma selecao de suporte do valor -- inclusive o fechamento de contorno, que
  // continua entrando como `(a, b)` e sendo redistribuido na montagem.  O que muda
  // e' so' o funcional avaliado: `pesos_derivada` em vez de `pesos`.
  //
  // E' O QUE A RECONSTRUCAO DE FLUXO DE ALTA ORDEM PEDE, e a razao de ele viver
  // aqui e nao no `Poisson`: o suporte depende apenas do PONTO, nao de quem
  // pergunta.  Logo as duas celulas que compartilham uma faceta obtem o MESMO
  // estencil para o mesmo ponto -- e o esquema fica conservativo por construcao,
  // sem ninguem ter de combinar nada.
  Estencil derivada_em(const Ponto3 &x, int dir) const;

  //! Quantos suportes foram TRUNCADOS pela franja desde o ultimo `zera()`.
  //
  // Zero e' o que se quer afirmar em qualquer numero de processos.  Diferente de
  // zero significa que os pesos perto da fronteira de particao nao sao os do
  // serial, e portanto que a solucao depende da particao.
  long truncados_pela_franja() const { return truncados_; }
  void zera_contadores_de_franja() const { truncados_ = 0; }

  //! O estencil do LAPLACIANO em `x`: ids, pesos e `rhs`.
  //
  // E' o que a familia de DIFERENCAS FINITAS usa para montar a linha, e nao tem
  // quadratura nenhuma: o laplaciano e' avaliado no ponto, nao integrado.
  //
  // Mesmo suporte e mesmo fechamento de contorno do valor e da derivada --
  // reusar a montagem nao e' economia, e' o que garante que o contorno entre do
  // mesmo jeito nos tres funcionais.  Divergir ali so' aparece perto da parede.
  Estencil laplaciano_em(const Ponto3 &x) const;

  //! Quantas CAMADAS de franja este construtor precisa que a malha tenha.
  //
  // E' `raio + 1`, e o raio sai da ALGEBRA: uma caixa de raio `r` tem `2r+1`
  // valores distintos de coordenada por direcao, e grau `g` precisa de `g+1`.
  //
  // QUEM PEDE E' O CHAMADOR, explicitamente -- `m.exige_camadas(ce.camadas_necessarias())`.
  // O construtor poderia pedir sozinho, e nao pede: a malha e' compartilhada, e
  // um objeto que a modifica na construcao faz a ordem de criacao importar sem
  // que isso apareca em lugar nenhum.
  int camadas_necessarias() const { return (ordem_ + 1) / 2 + 1; }

  //! Troca a estrategia de escolha do suporte.  O padrao e' `Caixa`, por
  //! medida: ver o comentario de `Suporte3`.
  void estrategia(Suporte3 e) { estrategia_ = e; }
  Suporte3 estrategia() const { return estrategia_; }

  //! O ajuste MLS por tras, para ler os contadores de rebaixamento de ordem.
  const Mls3 &mls() const { return mls_; }

  // O LOTE EM GPU ESTA' FORA DO ESCOPO DO LADO 3D, por decisao de escopo: o 3D
  // traz "so' o que o solver vai usar", e o lote e' otimizacao.  Quando entrar,
  // entra com o oraculo que o 2D tem -- comparacao ponto a ponto com a CPU, que
  // ja' pegou um nucleo usando minimos quadrados NAO ponderados.


  //! Diagnostico do suporte escolhido em `x`: quantos pontos de celula, quantos
  //! de contorno, quantos aneis, e o diametro.  Existe para a clausula F3 -- "o
  //! suporte cai na grade e e' local" -- poder ser afirmada de fora.
  //! O suporte, e a forma dele e' o que torna os tres tipos de contorno UM caminho.
  //
  // Cada ponto carrega `(id, a, b)`, com o significado `valor_no_ponto = a *
  // u[id] + b`:
  //
  //     celula comum          id = o id dela,        a = 1,  b = 0
  //     parede de Dirichlet   id = -1,               a = 0,  b = g
  //     fantasma (Neu/Robin)  id = a celula ESPELHADA, a = A,  b = B
  //
  // Na montagem, o peso `w` que o ajuste der ao ponto contribui `a*w` no id e
  // `b*w` no `rhs`.  Um laco, tres tipos -- e nada no ajuste sabe que algum ponto
  // e' especial.
  struct Suporte {
    std::vector<Ponto3> pontos;
    std::vector<long> ids;        //!< -1 quando o valor nao depende de incognita
    std::vector<double> a;        //!< coeficiente da incognita `ids[k]`
    std::vector<double> b;        //!< parte conhecida, vai para o `rhs`
    int aneis = 0;
    double diametro = 0.0;
    long fantasmas = 0;           //!< quantos pontos vieram de fantasma
    //! O suporte parou por falta de FRANJA, nao por ter fechado o anel.
    //
    // Em np=1 e' sempre falso.  Em paralelo fica verdadeiro perto da fronteira
    // de particao quando o suporte pedia mais de um anel alem das folhas locais:
    // a franja tem UMA camada, e o processo nao conhece os vizinhos dela.
    //
    // Nao e' erro -- e' uma propriedade declarada do arranjo --, mas e' uma
    // propriedade que muda os PESOS, logo o resultado passa a depender do numero
    // de processos.  Quem monta esquema de alta ordem tem de poder ver isto.
    bool truncado_pela_franja = false;

    long paredes = 0;             //!< quantos vieram de amostrar a parede
  };
  Suporte suporte_em(const Ponto3 &x) const;

  //! Monta o estencil a partir de um suporte e dos pesos do ajuste.
  //
  // Publica e estatica porque e' a traducao `(a, b, peso)` -> `(ids, pesos, rhs)`,
  // e ela e' a MESMA no caminho um-por-um e no do lote.  Duas copias divergiriam,
  // e divergir aqui significa o `rhs` certo num caminho e errado no outro.
  static Estencil monta(const Suporte &s, const std::vector<double> &w);

 private:
  Suporte3 estrategia_ = Suporte3::Caixa;

  //! As duas estrategias.  `suporte_em` despacha entre elas.
  Suporte suporte_caixa(const Ponto3 &x) const;
  Suporte suporte_aneis(const Ponto3 &x) const;

  mutable long truncados_ = 0;

  const Malha3 *m_;
  int ordem_;
  int max_pontos_;
  int max_aneis_;
  ContornoPorRegiao regioes_;
  mutable Mls3 mls_;
};

}  // namespace t8flow3d

#endif
