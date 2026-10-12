#include "solver2d.hxx"

#include <stdexcept>
#include <string>

#include "solver.h"

namespace t8flow2d {

namespace {

SolverEngine traduz(Motor m)
{
  switch (m) {
    case Motor::Petsc: return SOLVER_PETSC;
    case Motor::Hypre: return SOLVER_HYPRE;
    case Motor::Qualquer:
    default: return SOLVER_ANY;
  }
}

}  // namespace

Solver::Solver(Motor motor, int primeiro_id, std::size_t n_local, int max_nao_zeros)
    : motor_(motor), primeiro_(primeiro_id), n_(n_local)
{
  if (n_local == 0)
    throw std::invalid_argument("Solver: n_local = 0");
  if (max_nao_zeros < 1)
    throw std::invalid_argument("Solver: max_nao_zeros tem de ser >= 1");

  s_ = slv_create(traduz(motor), primeiro_id, n_local);
  if (s_ == nullptr)
    throw std::runtime_error("Solver: slv_create devolveu nulo");

  // Dimensiona a matriz ANTES de qualquer linha: no PETSc isto e' a
  // pre-alocacao, e no herdado ela e' obrigatoria.
  slv_set_maxnonzeros(s_, max_nao_zeros);

  linha_definida_.assign(n_local, false);
}

Solver::~Solver()
{
  if (s_ != nullptr) slv_destroy(s_);
}

std::string Solver::nome_do_motor() const
{
  if (motor_ == Motor::Qualquer) return "qualquer";

  // A unica aritmetica do deslocamento no t8Flow, e ela esta' aqui para ser
  // cobrada por teste: `engine_names[i]` corresponde a `(SolverEngine)(i+1)`.
  const int i = static_cast<int>(traduz(motor_)) - 1;
  const int n = static_cast<int>(sizeof(engine_names) / sizeof(engine_names[0]));
  if (i < 0 || i >= n) return "(fora da tabela)";
  return engine_names[i];
}

void Solver::linha(int id_global, const std::vector<int> &colunas,
                   const std::vector<double> &valores)
{
  if (montado_)
    throw std::logic_error("Solver::linha depois de monta()");
  if (colunas.size() != valores.size())
    throw std::invalid_argument("Solver::linha: colunas e valores de tamanhos diferentes");
  if (colunas.empty())
    throw std::invalid_argument("Solver::linha: linha vazia -- o herdado aborta nisso");

  const long local = static_cast<long>(id_global) - primeiro_;
  if (local < 0 || local >= static_cast<long>(n_))
    throw std::out_of_range(
        "Solver::linha: id " + std::to_string(id_global) + " fora de [" +
        std::to_string(primeiro_) + ", " + std::to_string(primeiro_ + static_cast<long>(n_) - 1) + "]");
  if (linha_definida_[static_cast<size_t>(local)])
    throw std::logic_error("Solver::linha: id " + std::to_string(id_global) +
                           " definido duas vezes -- o herdado afirma isso com assert");

  // O herdado recebe ponteiros nao-const.
  std::vector<int> cols = colunas;
  std::vector<double> vals = valores;
  slv_set_Ai(s_, id_global, static_cast<int>(cols.size()), cols.data(), vals.data());
  linha_definida_[static_cast<size_t>(local)] = true;
}

void Solver::lado_direito(int id_global, double v)
{
  slv_set_bi(s_, id_global, v);
}

void Solver::monta()
{
  // A GUARDA QUE EVITA UM ABORT TRES CAMADAS ABAIXO: o caminho de pre-alocacao
  // do PETSc afirma `col_ids[i] && vals[i]` para toda linha local.
  for (size_t i = 0; i < linha_definida_.size(); i++) {
    if (!linha_definida_[i])
      throw std::logic_error(
          "Solver::monta: a linha " + std::to_string(primeiro_ + static_cast<long>(i)) +
          " nao foi definida; sem ela o herdado aborta dentro da pre-alocacao");
  }
  slv_assemble(s_);     // rhs -> output -> matriz, nesta ordem
  montado_ = true;
}

void Solver::resolve()
{
  if (!montado_) throw std::logic_error("Solver::resolve antes de monta()");
  slv_solve(s_);
}

double Solver::x(int id_global) const
{
  return slv_get_xi(s_, id_global);
}

std::vector<double> Solver::solucao() const
{
  std::vector<double> v(n_, 0.0);
  for (size_t i = 0; i < n_; i++)
    v[i] = slv_get_xi(s_, primeiro_ + static_cast<int>(i));
  return v;
}

bool Solver::convergiu() const { return s_->converged; }
int Solver::iteracoes() const { return s_->SolverIterations; }
double Solver::residuo() const { return s_->SolverResidual; }

void Solver::tolerancia(double eps) { slv_set_eps(s_, eps); }
void Solver::max_iteracoes(int n) { slv_set_maxiteration(s_, n); }

}  // namespace t8flow2d
