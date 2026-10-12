// A SAIDA PARA VISUALIZACAO: um escritor so', dirigido por REGISTRO DE CAMPOS.
//
// POR QUE UM SO' e nao um por fisica: o HiGFlow tem mais de vinte variantes de
// `print_vtk` -- 2D, 3D, multifasico, viscoelastico, eletroosmotico, com e sem
// linhas de interface, serial e paralelo --, e cada fisica nova acrescenta
// mais uma.  O que varia entre elas nao e' o FORMATO, e' QUAIS CAMPOS entram.
// Entao o que varia e' dado, nao codigo: quem chama registra os campos e o
// escritor nao sabe que fisica e' essa.
//
// VTK E' VISUALIZACAO, NAO RESULTADO -- e isto esta' aqui porque o proprio
// HiGFlow o registra no cabecalho do IO dele, depois de medir: `np=4` pode
// diferir no VTK com o estado do solver BIT A BIT IDENTICO (as velocidades
// concordando em 0,0 sobre 15.920 facetas enquanto metricas tiradas do VTK
// mudavam por um fator de 115).  A razao e' que o escritor INTERPOLA -- o
// vetor vive nas facetas e o VTK o quer na celula --, e a interpolacao perto
// da fronteira de particao nao e' a mesma com e sem vizinho local.
//
// POR ISSO O CANAL DE RESULTADO E' OUTRO, e nao este: `Metricas` escreve YAML
// de precisao cheia, reduzido entre processos, no mesmo formato das
// referencias do HiGFlow.  Um teste construido sobre VTK mede o erro da saida,
// nao a fisica.  Quem usar este arquivo como oraculo vai ter o que mereceu.

#ifndef T8FLOW_COMUM_IO_SAIDA_VTK2D_HXX
#define T8FLOW_COMUM_IO_SAIDA_VTK2D_HXX

#include <string>
#include <vector>

#include "campos2d.hxx"
#include "malha2d.hxx"

namespace t8flow {

class SaidaVtk2d {
 public:
  explicit SaidaVtk2d(const t8flow2d::Malha &m);

  //! Escalar no CENTRO da celula.  Sai como esta', sem interpolacao.
  void escalar(const std::string &nome, const t8flow2d::CampoEscalar &c);

  //! Vetorial nas FACETAS.  Sai interpolado ao centro -- e e' esta
  //! interpolacao que torna o arquivo inadequado como oraculo.
  void vetorial(const std::string &nome, const t8flow2d::DominioFacetas &fd,
                const t8flow2d::CampoVetorial &v);

  //! Tensorial no centro, as quatro componentes.
  void tensorial(const std::string &nome, const t8flow2d::CampoTensor &t);

  //! Uma FRONTEIRA DE INTERFACE: a poligonal que a descreve.  Vai para um
  //! arquivo PROPRIO, e nao para o da malha: sao topologias diferentes
  //! (segmentos contra quadrilateros), e o VTK legado nao mistura as duas num
  //! arquivo sem que o leitor tenha de adivinhar qual e' qual.
  void fronteira(const std::string &nome,
                 const std::vector<std::vector<t8flow2d::Ponto2>> &curvas);

  //! Escreve `<prefixo>-<quadro>-<rank>.vtk`, e as fronteiras em
  //! `<prefixo>-<nome>-<quadro>-<rank>.vtk`.
  //
  // UM ARQUIVO POR PROCESSO, como o HiGFlow: juntar exigiria comunicar o campo
  // inteiro para um processo, que nao escala e -- pior -- faria a saida mudar
  // o custo da corrida que ela deveria so' observar.
  void escreve(const std::string &prefixo, int quadro) const;

  //! Esvazia os registros.  O escritor e' reutilizado entre quadros, e os
  //! campos sao registrados de novo -- as referencias podem ter mudado.
  void limpa();

 private:
  struct Campo {
    std::string nome;
    int componentes;                 //!< 1 escalar, 3 vetor, 9 tensor
    std::vector<double> v;           //!< por celula local, intercalado
  };

  const t8flow2d::Malha *m_;
  std::vector<Campo> campos_;
  std::vector<std::pair<std::string,
                        std::vector<std::vector<t8flow2d::Ponto2>>>> fronteiras_;
};

}  // namespace t8flow

#endif  // T8FLOW_COMUM_IO_SAIDA_VTK2D_HXX
