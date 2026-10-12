// O LOTE DE ESTENCEIS: N consultas numa chamada, com back-end escolhido em
// EXECUCAO.
//
// Pedido do usuario em duas partes, e a segunda e' a que manda: "vários stencils
// de uma vez" e "uma versao em CPU, caso nao tenha GPU".  Logo:
//
//   * **a CPU e' implementacao de PRIMEIRA CLASSE, nao plano B.**  Ela funciona
//     sem OpenCL instalado, sem OpenCL compilado, e sem placa -- e e' ela que
//     define o resultado correto.  Quem nao tem GPU nao perde funcionalidade
//     nenhuma: perde velocidade, e nada mais;
//   * **a GPU e' opcional nas duas pontas**: ou nao existe em tempo de compilacao
//     (sem `T8FLOW_COM_OPENCL`), ou existe mas sem dispositivo em tempo de
//     execucao.  Os dois casos caem na CPU, e `back_end()` DIZ qual rodou;
//   * **a CPU e' o ORACULO da GPU.**  Lote em placa que erra em silencio e' o pior
//     caso possivel para este projeto, e a unica defesa e' poder recalcular na CPU
//     e comparar.  E' por isso que o back-end e' trocavel em execucao e nao em
//     compilacao: sem isso, comparar exigiria dois binarios.
//
// A FORMA DO RESULTADO segue o `sim_stencil` da HiGTree: ids, pesos, e um `rhs`
// SEPARADO para o que nao se associa a nenhuma celula -- o valor de Dirichlet do
// fechamento de contorno.  Guardar o rhs junto com os pesos e' o que obriga o
// chamador a saber qual id e' "a parede", e e' assim que o fechamento se perde.
//
// O QUE AINDA NAO ESTA' AQUI, e esta' dito para nao ser suposto: a ESCOLHA do
// suporte.  Quem chama entrega os pontos e os ids; crescer o raio ate' satisfazer
// a ordem pedida, e o fechamento de contorno que preenche o `rhs`, sao a frente
// seguinte (secao 13.1 do doc/projeto-t8flow.md).  O lote calcula os PESOS de um
// suporte dado, que e' a parte que custa e a parte que a placa acelera.
#ifndef T8FLOW_LOTE2D_HXX
#define T8FLOW_LOTE2D_HXX

#include <cstddef>
#include <string>
#include <vector>

#include "mls2d.hxx"

namespace t8flow2d {

//! Uma consulta: onde interpolar, com que suporte, e os ids daquelas celulas.
struct ConsultaEstencil {
  Ponto2 alvo;
  std::vector<Ponto2> suporte;
  std::vector<long> ids;        //!< ids globais, na MESMA ordem do suporte
};

//! O estencil de uma consulta, no espirito do `sim_stencil` da HiGTree.
struct Estencil {
  std::vector<long> ids;
  std::vector<double> pesos;
  double rhs = 0.0;             //!< o que nao se associa a celula (Dirichlet)
};

enum class BackEnd { Cpu, Gpu };

class LoteEstenceis {
 public:
  //! `ordem` e' o grau reproduzido exatamente; `max_pontos` o teto do suporte.
  LoteEstenceis(int ordem, int max_pontos);
  ~LoteEstenceis();

  LoteEstenceis(const LoteEstenceis &) = delete;
  LoteEstenceis &operator=(const LoteEstenceis &) = delete;

  //! Ha' GPU utilizavel AGORA: compilado com OpenCL, dispositivo presente, e a
  //! ordem pedida suportada pelo kernel.  Falso em qualquer um dos tres casos.
  bool gpu_disponivel() const;

  //! Nome do dispositivo de GPU, ou por que nao ha'.  Diagnostico legivel: "sem
  //! OpenCL em tempo de compilacao" e "OpenCL presente sem dispositivo" sao
  //! estados diferentes, e confundi-los custa um dia.
  std::string diagnostico_gpu() const;

  //! Back-end em uso.  O PADRAO vem de quem compilou (`make GPU=sim` faz a GPU
  //! ser o padrao); sem isso, comeca em Cpu.
  BackEnd back_end() const { return back_end_; }

  //! A GPU foi PEDIDA em tempo de compilacao (`-DT8FLOW_GPU_PADRAO`, de
  //! `make GPU=sim`).  Distinto de `gpu_disponivel()`: pedir e ter sao coisas
  //! diferentes, e o construtor avisa em stderr quando pediram e nao ha'.
  static bool gpu_pedida_na_compilacao();

  //! Troca o back-end.  Pedir `Gpu` sem disponibilidade LANCA -- nao cai na CPU
  //! em silencio, porque quem pediu GPU explicitamente quer saber que nao ha'.
  void usa_back_end(BackEnd b);

  //! Calcula os estenceis de todas as consultas.
  //
  // `saida` e' redimensionada.  Consulta com suporte abaixo do minimo da ordem ou
  // acima de `max_pontos` faz a chamada LANCAR, identificando o indice -- o lote
  // nao calcula parte e cala sobre o resto.
  void calcula(const std::vector<ConsultaEstencil> &consultas,
               std::vector<Estencil> &saida) const;

  int ordem() const { return ordem_; }
  int max_pontos() const { return max_pontos_; }

  //! Monomios da base para esta ordem em 2D: 3 para ordem 1, 6 para ordem 2.
  static int num_monomios(int ordem);

 private:
  void calcula_cpu(const std::vector<ConsultaEstencil> &, std::vector<Estencil> &) const;
  void calcula_gpu(const std::vector<ConsultaEstencil> &, std::vector<Estencil> &) const;

  int ordem_;
  int max_pontos_;
  BackEnd back_end_ = BackEnd::Cpu;
  mutable Mls mls_;             //!< o back-end de CPU E' o nucleo herdado
  struct Gpu;
  Gpu *gpu_ = nullptr;          //!< nulo quando nao ha' GPU, e isso e' normal
};

}  // namespace t8flow2d

#endif
