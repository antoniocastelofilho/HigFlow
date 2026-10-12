// O SOLVER 2D: resolve um sistema de solucao CONHECIDA, e recusa os usos que
// abortariam tres camadas abaixo.
//
// ORACULO: a solucao, nao o residuo.  O residuo e' calculado pelo proprio solver
// e um motor que pare na iteracao maxima tambem devolve; o criterio tem de ser
// externo.  Aqui o sistema e' a tridiagonal de Poisson 1D com -u''=2 e u(0)=u(1)=0,
// cuja solucao do SISTEMA LINEAR e' exatamente x_i = i*(n+1-i)*h^2 -- identidade
// algebrica, nao a solucao continua truncada, logo a tolerancia pode ser apertada.
#include "../testing.h"

#include <stdexcept>
#include <vector>

#include "arranque.hxx"
#include "solver2d.hxx"

using t8flow2d::Motor;
using t8flow2d::Solver;

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);

  const int n = 32;
  const double h = 1.0 / (n + 1);

  // A aritmetica do DESLOCAMENTO entre as duas listas do herdado, cobrada em
  // execucao: `engine_names[i]` corresponde a `(SolverEngine)(i+1)`, e nada
  // confere isso em tempo de compilacao.  Errar faz o programa anunciar um motor
  // e rodar outro.
  t_case("o_nome_do_motor_pedido_bate_com_a_tabela");
  {
    Solver sp(Motor::Petsc, 0, 4, 3);
    T_CHECK_MSG(sp.nome_do_motor() == "petsc",
                "Motor::Petsc deu \"%s\", esperava \"petsc\"", sp.nome_do_motor().c_str());
    Solver sh(Motor::Hypre, 0, 4, 3);
    T_CHECK_MSG(sh.nome_do_motor() == "hypre",
                "Motor::Hypre deu \"%s\", esperava \"hypre\"", sh.nome_do_motor().c_str());
  }

  Solver s(Motor::Petsc, 0, static_cast<std::size_t>(n), 3);
  s.max_iteracoes(1000);
  s.tolerancia(1e-14);

  for (int i = 0; i < n; i++) {
    std::vector<int> cols;
    std::vector<double> vals;
    if (i > 0)     { cols.push_back(i - 1); vals.push_back(-1.0 / (h * h)); }
                     cols.push_back(i);     vals.push_back( 2.0 / (h * h));
    if (i < n - 1) { cols.push_back(i + 1); vals.push_back(-1.0 / (h * h)); }
    s.linha(i, cols, vals);
    s.lado_direito(i, 2.0);
  }

  t_case("linha_repetida_e_recusada");
  {
    bool lancou = false;
    try {
      s.linha(0, {0}, {1.0});
    } catch (const std::logic_error &) {
      lancou = true;
    }
    T_CHECK_MSG(lancou, "definir a mesma linha duas vezes passou sem reclamar");
  }

  t_case("id_fora_da_faixa_local_e_recusado");
  {
    bool lancou = false;
    try {
      s.linha(n + 10, {0}, {1.0});
    } catch (const std::out_of_range &) {
      lancou = true;
    }
    T_CHECK_MSG(lancou, "id fora de [0, %d] passou sem reclamar", n - 1);
  }

  s.monta();
  s.resolve();

  t_case("solver_resolve_sistema_de_solucao_conhecida");
  T_CHECK_MSG(s.convergiu(), "nao convergiu: %d iteracoes, residuo %.3g",
              s.iteracoes(), s.residuo());
  double pior = 0.0;
  int pior_i = -1;
  for (int i = 0; i < n; i++) {
    const double esperado = (i + 1) * static_cast<double>(n - i) * h * h;
    const double e = std::fabs(s.x(i) - esperado);
    if (e > pior) { pior = e; pior_i = i; }
  }
  T_CHECK_MSG(pior <= 1e-10, "pior desvio %.3g na linha %d (tolerancia 1e-10)", pior, pior_i);

  t_case("solucao_inteira_bate_com_a_consulta_individual");
  {
    const std::vector<double> x = s.solucao();
    T_CHECK_MSG(x.size() == static_cast<std::size_t>(n), "%zu valores para %d incognitas",
                x.size(), n);
    double dif = 0.0;
    for (int i = 0; i < n; i++) dif = std::fmax(dif, std::fabs(x[i] - s.x(i)));
    T_NEAR(dif, 0.0, 0.0, "maior diferenca entre solucao() e x(i)");
  }

  // A GUARDA QUE EVITA UM ABORT DENTRO DO HERDADO: o caminho de pre-alocacao do
  // PETSc afirma `col_ids[i] && vals[i]` para TODA linha local.  Linha nao
  // definida nao e' sistema singular, e' abort tres camadas abaixo.
  t_case("linha_faltando_e_recusada_antes_de_montar");
  {
    Solver incompleto(Motor::Petsc, 0, 3, 3);
    incompleto.linha(0, {0}, {1.0});
    incompleto.linha(2, {2}, {1.0});
    bool lancou = false;
    try {
      incompleto.monta();
    } catch (const std::logic_error &) {
      lancou = true;
    }
    T_CHECK_MSG(lancou, "montar com a linha 1 ausente passou sem reclamar");
  }

  return t_end();
}
