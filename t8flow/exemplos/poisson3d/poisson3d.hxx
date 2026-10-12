// O OPERADOR DE POISSON EM TRES DIMENSOES, em volumes finitos.
//
// ESCOPO DECIDIDO PELO USUARIO -- "so' o que o solver vai usar": as duas
// variantes de VOLUMES FINITOS, e nenhuma de diferencas finitas.  Os dois
// esquemas compactos recusam malha adaptada (sao contraste, nao ferramenta) e o
// laplaciano pontual e' instavel nela.
//
//   vf-fluxo2       balanco de fluxos de segunda ordem
//   vf-mlsgauss     derivada normal por ajuste MLS mais quadratura de Gauss,
//                   exato ate' o grau da ordem do ajuste
//
// ESTE ARQUIVO E' EXEMPLO E NAO BIBLIOTECA.  Nada em `src/` sabe que equacao se
// esta' resolvendo, e este resolve UMA.  O passo de projecao de Navier-Stokes
// sera' um Poisson, e o operador que o t8SolverNS usar nascera' da extracao
// deste -- a parte que monta e resolve, sem as funcoes de erro e sem a fonte
// manufaturada, que sao verificacao.
//
// O QUE MUDA DE DUAS PARA TRES DIMENSOES, e nao e' so' um laco a mais:
//
//   A FACETA TEM AREA, nao comprimento.  O fluxo e' `|face| * du/dn`, e `|face|`
//   e' o produto dos dois lados transversais.  Dimensao errada ali ainda
//   converge -- para a equacao errada;
//
//   NA TRANSICAO 2:1 UMA FACE TEM QUATRO VIZINHOS, nao dois, e a razao de areas
//   e' `1/4` e nao `1/2`;
//
//   A QUADRATURA DE GAUSS NA FACETA E' 2x2, nao 2, e no volume e' 2x2x2.

#ifndef T8FLOW_EX_POISSON3D_HXX
#define T8FLOW_EX_POISSON3D_HXX

#include <functional>
#include <map>
#include <string>
#include <vector>

#include "contorno3d.hxx"
#include "discretizacao3d.hxx"
#include "malha3d.hxx"
#include "solverlin.hxx"

namespace t8flow3d {

using t8flow::Motor;
using t8flow::Solver;

class Poisson {
 public:
  Poisson(const Malha3 &m, ContornoPorRegiao bc);

  //! O termo de fonte `f` em `lap(u) = f`.
  void fonte(std::function<double(const Ponto3 &)> f) { f_ = std::move(f); }

  //! A discretizacao -- familia, variante e ordem de uma vez.
  void discretizacao(const Discretizacao &d);
  Discretizacao discretizacao() const { return disc_; }

  //! Teto do suporte e aneis do MLS.  O numero de aneis decide quantas camadas
  //! de franja a malha precisa, e o `resolve` as pede.
  void mls(int max_pontos) { mls_maxpts_ = max_pontos; }

  void resolve(Motor motor, double tol = 1e-12, int max_it = 20000);

  double u(long id) const;
  const std::vector<double> &solucao() const { return sol_; }
  bool convergiu() const { return convergiu_; }
  int iteracoes() const { return iteracoes_; }
  double residuo() const { return residuo_; }

  //! Os tres tempos de `resolve`, SEPARADOS -- e sao tres, nao dois como no
  //! lado bidimensional.
  //
  // SEPARADOS PORQUE ESCALAM DE FORMAS DIFERENTES, e a soma esconde isso:
  //
  //   a FRANJA e' COLETIVA.  Custa o cubo do raio em tres dimensoes, vale zero
  //   em um processo e cresce com o numero deles -- e' a unica fase que pode
  //   piorar com paralelismo;
  //
  //   a MONTAGEM nao comunica nada: cada processo monta as suas linhas.  Se ela
  //   domina em np=1, o tempo total parece escalar bem sem que a parte dificil
  //   escale;
  //
  //   a RESOLUCAO comunica a cada iteracao e e' a que escala pior.
  //
  // Em duas dimensoes o tempo somado chegou a sugerir ganho de 29,8x indo de um
  // para oito processos -- superlinear, impossivel.  A causa era a montagem ser
  // QUADRATICA, e so' apareceu quando as fases foram medidas separadas.
  //
  // Em segundos; valem zero antes do primeiro `resolve`.
  double tempo_de_franja() const { return t_franja_; }
  double tempo_de_montagem() const { return t_monta_; }
  double tempo_de_resolucao() const { return t_resolve_; }

  //! Erro contra a solucao exata, nas duas normas.  REDUZEM entre processos: sem
  //! isso cada um devolveria o pior erro da sua particao, e o numero sairia
  //! MENOR que o verdadeiro.
  double erro_inf(const std::function<double(const Ponto3 &)> &exata) const;
  double erro_l2(const std::function<double(const Ponto3 &)> &exata) const;

  //! O ORACULO DIRETO DO OPERADOR: `max |A u_exata - b|`, sem passar pelo solver.
  //
  // Separa CONSISTENCIA de CONVERGENCIA: um operador pode ser consistente e o
  // sistema divergir por causa do contorno, ou o contrario.  E a montagem e' a
  // MESMA funcao que o solver usa -- duas montagens divergiriam exatamente onde
  // e' dificil ver, perto da parede.
  double residuo_da_exata(const std::function<double(const Ponto3 &)> &exata) const;

  //! Quantos ajustes MLS rebaixaram a ordem, e quantos suportes a franja
  //! truncou.  Zero e' o que se quer poder afirmar nos dois.
  //! DESPEJA cada linha montada em `<prefixo>.<rank>`, no proximo `resolve`.
  //
  // Diagnostico, e usa o MESMO caminho de montagem que o solver -- uma segunda
  // montagem so' para despejar divergiria exatamente onde e' dificil ver.
  //
  // Serve para uma pergunta que nenhuma norma responde: a matriz montada em
  // paralelo e' a MESMA que a serial?  O oraculo `residuo_da_exata` e' um
  // MAXIMO, e um maximo dominado pelas linhas de parede nao se move quando o
  // interior muda.
  void despeja_operador(const std::string &prefixo) { despejo_ = prefixo; }

  long rebaixamentos_mls() const { return rebaixos_; }
  long truncados_pela_franja() const { return truncados_; }

 private:
  void monta_linha(long b, long i, const void *construtor,
                   std::map<long, double> *linha, double *rhs) const;
  void monta_linha_fluxo(long b, long i, std::map<long, double> *linha,
                         double *rhs) const;
  void monta_linha_mls(long b, long i, const void *construtor,
                       std::map<long, double> *linha, double *rhs) const;

  double reduz_max(double v) const;
  double reduz_soma(double v) const;

  const Malha3 *m_;
  ContornoPorRegiao contorno_;
  std::function<double(const Ponto3 &)> f_;
  Discretizacao disc_ = Discretizacao::padrao();
  int mls_maxpts_ = 200;

  long primeiro_ = 0;
  long n_ = 0;
  std::vector<double> sol_;
  bool convergiu_ = false;
  int iteracoes_ = 0;
  double residuo_ = 0.0;
  std::string despejo_;
  double t_franja_ = 0.0;
  double t_monta_ = 0.0;
  double t_resolve_ = 0.0;
  mutable long rebaixos_ = 0;
  mutable long truncados_ = 0;
};

}  // namespace t8flow3d

#endif  // T8FLOW_EX_POISSON3D_HXX
