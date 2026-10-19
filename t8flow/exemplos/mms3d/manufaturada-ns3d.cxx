#include "manufaturada-ns3d.hxx"

#include <cmath>

namespace t8flow3d {

namespace {
constexpr double PI = 3.14159265358979323846;
inline double S(double a) { return std::sin(PI * a); }
inline double C(double a) { return std::cos(PI * a); }
}  // namespace

double ManufaturadaNS3::u(const Ponto3 &x, int i, double t) const
{
  const double T = std::exp(-t);
  switch (i) {
    case 0: return PI * S(x[0]) * (C(x[1]) * S(x[2]) - S(x[1]) * C(x[2])) * T;
    case 1: return -PI * C(x[0]) * S(x[1]) * S(x[2]) * T;
    default: return PI * C(x[0]) * S(x[1]) * S(x[2]) * T;
  }
}

double ManufaturadaNS3::dudx(const Ponto3 &x, int i, int j, double t) const
{
  // AS NOVE, UMA A UMA.  Cada linha e' a derivada do termo correspondente de
  // `u`, e o fator `pi` extra vem da regra da cadeia.
  const double T = std::exp(-t);
  const double sx = S(x[0]), cx = C(x[0]);
  const double sy = S(x[1]), cy = C(x[1]);
  const double sz = S(x[2]), cz = C(x[2]);
  if (i == 0) {
    // u_x = pi S_x (C_y S_z - S_y C_z) T
    if (j == 0) return PI * PI * cx * (cy * sz - sy * cz) * T;
    if (j == 1) return PI * PI * sx * (-sy * sz - cy * cz) * T;
    return PI * PI * sx * (cy * cz + sy * sz) * T;
  }
  if (i == 1) {
    // u_y = -pi C_x S_y S_z T
    if (j == 0) return PI * PI * sx * sy * sz * T;
    if (j == 1) return -PI * PI * cx * cy * sz * T;
    return -PI * PI * cx * sy * cz * T;
  }
  // u_z = +pi C_x S_y S_z T  -- simetrico de u_y com o sinal trocado
  if (j == 0) return -PI * PI * sx * sy * sz * T;
  if (j == 1) return PI * PI * cx * cy * sz * T;
  return PI * PI * cx * sy * cz * T;
}

double ManufaturadaNS3::dudt(const Ponto3 &x, int i, double t) const
{
  // `T = exp(-t)` -> `dT/dt = -T`, logo `du/dt = -u`.
  return -u(x, i, t);
}

double ManufaturadaNS3::lap(const Ponto3 &x, int i, double t) const
{
  // CADA TERMO E' UM PRODUTO DE TRES TRIGONOMETRICAS de argumento `pi .`, e o
  // laplaciano de um desses e' `-3 pi^2` vezes ele.  Vale componente a
  // componente, e e' por isso que o potencial foi escolhido assim.
  return -3.0 * PI * PI * u(x, i, t);
}

double ManufaturadaNS3::p(const Ponto3 &x) const
{
  return C(x[0]) * C(x[1]) * C(x[2]);
}

double ManufaturadaNS3::dpdx(const Ponto3 &x, int i) const
{
  const double sx = S(x[0]), cx = C(x[0]);
  const double sy = S(x[1]), cy = C(x[1]);
  const double sz = S(x[2]), cz = C(x[2]);
  if (i == 0) return -PI * sx * cy * cz;
  if (i == 1) return -PI * cx * sy * cz;
  return -PI * cx * cy * sz;
}

double ManufaturadaNS3::f(const Ponto3 &x, int i, double t) const
{
  // `f = d_t u + (u.grad)u - lap(u)/Re + grad p`.
  //
  // O CONVECTIVO E' MONTADO AQUI das derivadas exatas, e isto NAO e'
  // diferenciacao numerica: as `dudx` sao formulas a mao.  Montar a soma por
  // codigo e' algebra, e tira do caminho a chance de eu errar a conta de nove
  // produtos escrevendo-os um a um.
  double conv = 0.0;
  for (int j = 0; j < 3; j++) conv += u(x, j, t) * dudx(x, i, j, t);
  return dudt(x, i, t) + conv - lap(x, i, t) / re + dpdx(x, i);
}

}  // namespace t8flow3d
