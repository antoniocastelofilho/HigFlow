#ifndef HIG_FLOW_PROBLEM
#define HIG_FLOW_PROBLEM

#include "higtree.h"

// ---------------------------------------------------------------------------
// A interface que cada exemplo implementa para descrever o SEU problema.
//
// Ate' aqui isso era feito por oito ponteiros de funcao livres, registrados de
// uma vez por higflow_set_external_functions e guardados em ns->func.  O arranjo
// funciona, mas nao consegue exprimir que as oito pertencem ao mesmo problema:
// nao ha' estado compartilhado entre elas, nao ha' como um exemplo especializar
// so' uma parte, e o compilador nao verifica que o conjunto esta completo.
//
// A migracao e' incremental por construcao.  higflow_set_external_functions
// continua existindo e passa a instalar um HigFlowLegacyProblem, que reenvia as
// chamadas aos mesmos ponteiros -- entao um exemplo nao convertido nao muda uma
// linha.  A biblioteca chama sempre ns->problem->metodo(), sem saber qual dos
// dois caminhos esta por tras.
// ---------------------------------------------------------------------------
class HigFlowProblem {
public:
    virtual ~HigFlowProblem() {}

    // Condicao inicial e termos de fonte no interior
    virtual real pressure(Point center, real t) = 0;
    virtual real velocity(Point center, int dim, real t) = 0;
    virtual real source_term(Point center, real t) = 0;
    virtual real facet_source_term(Point center, int dim, real t) = 0;

    // Os mesmos, nas fronteiras identificadas por id
    virtual real boundary_pressure(int id, Point center, real t) = 0;
    virtual real boundary_velocity(int id, Point center, int dim, real t) = 0;
    virtual real boundary_source_term(int id, Point center, real t) = 0;
    virtual real boundary_facet_source_term(int id, Point center, int dim, real t) = 0;
};

// ---------------------------------------------------------------------------
// Adaptador para os exemplos que ainda registram por ponteiro de funcao.
// Guarda os oito e reenvia.  Enquanto existir um exemplo nao convertido, esta
// classe precisa existir; quando o ultimo for convertido, ela e
// higflow_set_external_functions saem juntas.
// ---------------------------------------------------------------------------
class HigFlowLegacyProblem : public HigFlowProblem {
public:
    real (*fn_pressure)(Point, real);
    real (*fn_velocity)(Point, int, real);
    real (*fn_source_term)(Point, real);
    real (*fn_facet_source_term)(Point, int, real);
    real (*fn_boundary_pressure)(int, Point, real);
    real (*fn_boundary_velocity)(int, Point, int, real);
    real (*fn_boundary_source_term)(int, Point, real);
    real (*fn_boundary_facet_source_term)(int, Point, int, real);

    real pressure(Point c, real t)                  { return fn_pressure(c, t); }
    real velocity(Point c, int d, real t)           { return fn_velocity(c, d, t); }
    real source_term(Point c, real t)               { return fn_source_term(c, t); }
    real facet_source_term(Point c, int d, real t)  { return fn_facet_source_term(c, d, t); }
    real boundary_pressure(int id, Point c, real t)                 { return fn_boundary_pressure(id, c, t); }
    real boundary_velocity(int id, Point c, int d, real t)          { return fn_boundary_velocity(id, c, d, t); }
    real boundary_source_term(int id, Point c, real t)              { return fn_boundary_source_term(id, c, t); }
    real boundary_facet_source_term(int id, Point c, int d, real t) { return fn_boundary_facet_source_term(id, c, d, t); }
};

// ---------------------------------------------------------------------------
// Segunda porta: os callbacks especificos de modelo, registrados pelos
// higflow_create_domain_*.  Ao contrario dos oito da primeira porta, estes NAO
// sao um conjunto universal -- cada exemplo usa so' os do modelo que configura, e
// o mesmo nome tem assinaturas diferentes entre modelos (get_viscosity e'
// (Point,q,t) no generalized-newtonian e (Point,q,t,beta,struct_par) no
// variable-viscosity).  Por isso sao interfaces SEPARADAS por modelo, herdadas so'
// por quem precisa, em vez de uma interface gorda com metodos vazios.
// ---------------------------------------------------------------------------
class HigFlowViscoelasticProblem {
public:
    virtual ~HigFlowViscoelasticProblem() {}
    virtual real tensor(Point center, int i, int j, real t) = 0;
    virtual real kernel(int dim, real lambda, real tol) = 0;
    virtual real kernel_inverse(int dim, real lambda, real tol) = 0;
    virtual real kernel_jacobian(int dim, real lambda, real tol) = 0;
};

// Adaptador para quem ainda registra por ponteiro.
class HigFlowLegacyViscoelastic : public HigFlowViscoelasticProblem {
public:
    real (*fn_tensor)(Point, int, int, real);
    real (*fn_kernel)(int, real, real);
    real (*fn_kernel_inverse)(int, real, real);
    real (*fn_kernel_jacobian)(int, real, real);

    real tensor(Point c, int i, int j, real t)      { return fn_tensor(c, i, j, t); }
    real kernel(int d, real l, real tol)            { return fn_kernel(d, l, tol); }
    real kernel_inverse(int d, real l, real tol)    { return fn_kernel_inverse(d, l, tol); }
    real kernel_jacobian(int d, real l, real tol)   { return fn_kernel_jacobian(d, l, tol); }
};

// --- multifasico: viscosidades e densidades das duas fases, e a fracao volumetrica
class HigFlowMultiphaseProblem {
public:
    virtual ~HigFlowMultiphaseProblem() {}
    virtual real viscosity0(Point center, real t) = 0;
    virtual real viscosity1(Point center, real t) = 0;
    virtual real density0(Point center, real t) = 0;
    virtual real density1(Point center, real t) = 0;
    virtual real fracvol(Point center, Point delta, real t) = 0;
};

class HigFlowLegacyMultiphase : public HigFlowMultiphaseProblem {
public:
    real (*fn_viscosity0)(Point, real);
    real (*fn_viscosity1)(Point, real);
    real (*fn_density0)(Point, real);
    real (*fn_density1)(Point, real);
    real (*fn_fracvol)(Point, Point, real);

    real viscosity0(Point c, real t)            { return fn_viscosity0(c, t); }
    real viscosity1(Point c, real t)            { return fn_viscosity1(c, t); }
    real density0(Point c, real t)              { return fn_density0(c, t); }
    real density1(Point c, real t)              { return fn_density1(c, t); }
    real fracvol(Point c, Point d, real t)      { return fn_fracvol(c, d, t); }
};

// --- multifasico viscoelastico.  kernel/kernel_inverse/kernel_jacobian tem a
// mesma assinatura da interface viscoelastica monofasica: um exemplo que herde as
// duas e defina uma vez sobrepoe as duas, que e' o comportamento desejado quando
// ele passa a mesma funcao para os dois registros.
class HigFlowMultiphaseViscoelasticProblem {
public:
    virtual ~HigFlowMultiphaseViscoelasticProblem() {}
    virtual real tensor_multiphase(real fracvol, Point center, int i, int j, real t) = 0;
    virtual real kernel(int dim, real lambda, real tol) = 0;
    virtual real kernel_inverse(int dim, real lambda, real tol) = 0;
    virtual real kernel_jacobian(int dim, real lambda, real tol) = 0;
};

class HigFlowLegacyMultiphaseViscoelastic : public HigFlowMultiphaseViscoelasticProblem {
public:
    real (*fn_tensor_multiphase)(real, Point, int, int, real);
    real (*fn_kernel)(int, real, real);
    real (*fn_kernel_inverse)(int, real, real);
    real (*fn_kernel_jacobian)(int, real, real);

    real tensor_multiphase(real f, Point c, int i, int j, real t) { return fn_tensor_multiphase(f, c, i, j, t); }
    real kernel(int d, real l, real tol)          { return fn_kernel(d, l, tol); }
    real kernel_inverse(int d, real l, real tol)  { return fn_kernel_inverse(d, l, tol); }
    real kernel_jacobian(int d, real l, real tol) { return fn_kernel_jacobian(d, l, tol); }
};

// --- generalized newtonian: so' a viscosidade
class HigFlowGenNewtonianProblem {
public:
    virtual ~HigFlowGenNewtonianProblem() {}
    virtual real viscosity(Point center, real q, real t) = 0;
};
class HigFlowLegacyGenNewtonian : public HigFlowGenNewtonianProblem {
public:
    real (*fn_viscosity)(Point, real, real);
    real viscosity(Point c, real q, real t) { return fn_viscosity(c, q, t); }
};

// --- viscoelastico integral: so' o tensor inicial
class HigFlowIntegralProblem {
public:
    virtual ~HigFlowIntegralProblem() {}
    virtual real tensor(Point center, int i, int j, real t) = 0;
};
class HigFlowLegacyIntegral : public HigFlowIntegralProblem {
public:
    real (*fn_tensor)(Point, int, int, real);
    real tensor(Point c, int i, int j, real t) { return fn_tensor(c, i, j, t); }
};

// --- viscoelastico com viscosidade variavel.  get_viscosity aqui NAO e' o mesmo
// do generalized-newtonian: leva beta e struct_par a mais.  Foi essa colisao de
// nome com assinaturas diferentes que impos interfaces separadas por modelo.
class HigFlowVariableViscosityProblem {
public:
    virtual ~HigFlowVariableViscosityProblem() {}
    virtual real tensor(Point center, int i, int j, real t) = 0;
    virtual real kernel(int dim, real lambda, real tol) = 0;
    virtual real kernel_inverse(int dim, real lambda, real tol) = 0;
    virtual real kernel_jacobian(int dim, real lambda, real tol) = 0;
    virtual real viscosity(Point center, real q, real t, real beta, real struct_par) = 0;
    virtual real structpar(Point center, real q, real t, real beta,
                           real Phi, real Lambda, real Gamma) = 0;
};

class HigFlowLegacyVariableViscosity : public HigFlowVariableViscosityProblem {
public:
    real (*fn_tensor)(Point, int, int, real);
    real (*fn_kernel)(int, real, real);
    real (*fn_kernel_inverse)(int, real, real);
    real (*fn_kernel_jacobian)(int, real, real);
    real (*fn_viscosity)(Point, real, real, real, real);
    real (*fn_structpar)(Point, real, real, real, real, real, real);

    real tensor(Point c, int i, int j, real t)    { return fn_tensor(c, i, j, t); }
    real kernel(int d, real l, real tol)          { return fn_kernel(d, l, tol); }
    real kernel_inverse(int d, real l, real tol)  { return fn_kernel_inverse(d, l, tol); }
    real kernel_jacobian(int d, real l, real tol) { return fn_kernel_jacobian(d, l, tol); }
    real viscosity(Point c, real q, real t, real b, real sp) { return fn_viscosity(c, q, t, b, sp); }
    real structpar(Point c, real q, real t, real b, real P, real L, real G) {
        return fn_structpar(c, q, t, b, P, L, G);
    }
};

// --- eletro-osmotico monofasico.  Onze metodos: cinco no interior, cinco na
// fronteira, mais a permissividade.
class HigFlowElectroosmoticProblem {
public:
    virtual ~HigFlowElectroosmoticProblem() {}
    virtual real source_term(Point center, int dim, real t) = 0;
    virtual real phi(Point center, real t) = 0;
    virtual real psi(Point center, real t) = 0;
    virtual real nplus(Point center, real t) = 0;
    virtual real nminus(Point center, real t) = 0;
    virtual real boundary_source_term(int id, Point center, int dim, real t) = 0;
    virtual real boundary_phi(int id, Point center, real t) = 0;
    virtual real boundary_psi(int id, Point center, real t) = 0;
    virtual real boundary_nplus(int id, Point center, real t) = 0;
    virtual real boundary_nminus(int id, Point center, real t) = 0;
    virtual real permittivity(Point center, real t) = 0;
};

class HigFlowLegacyElectroosmotic : public HigFlowElectroosmoticProblem {
public:
    real (*fn_source_term)(Point, int, real);
    real (*fn_phi)(Point, real);
    real (*fn_psi)(Point, real);
    real (*fn_nplus)(Point, real);
    real (*fn_nminus)(Point, real);
    real (*fn_boundary_source_term)(int, Point, int, real);
    real (*fn_boundary_phi)(int, Point, real);
    real (*fn_boundary_psi)(int, Point, real);
    real (*fn_boundary_nplus)(int, Point, real);
    real (*fn_boundary_nminus)(int, Point, real);
    real (*fn_permittivity)(Point, real);

    real source_term(Point c, int d, real t)   { return fn_source_term(c, d, t); }
    real phi(Point c, real t)                  { return fn_phi(c, t); }
    real psi(Point c, real t)                  { return fn_psi(c, t); }
    real nplus(Point c, real t)                { return fn_nplus(c, t); }
    real nminus(Point c, real t)               { return fn_nminus(c, t); }
    real boundary_source_term(int id, Point c, int d, real t) { return fn_boundary_source_term(id, c, d, t); }
    real boundary_phi(int id, Point c, real t)     { return fn_boundary_phi(id, c, t); }
    real boundary_psi(int id, Point c, real t)     { return fn_boundary_psi(id, c, t); }
    real boundary_nplus(int id, Point c, real t)   { return fn_boundary_nplus(id, c, t); }
    real boundary_nminus(int id, Point c, real t)  { return fn_boundary_nminus(id, c, t); }
    real permittivity(Point c, real t)             { return fn_permittivity(c, t); }
};

// --- eletro-osmotico multifasico.  Os mesmos onze, com fracvol na frente: o
// valor depende da fase local.  Um exemplo que herde as duas interfaces fica com
// SOBRECARGAS, nao com colisao, porque as assinaturas diferem.
class HigFlowMultiphaseElectroosmoticProblem {
public:
    virtual ~HigFlowMultiphaseElectroosmoticProblem() {}
    virtual real source_term(real fracvol, Point center, int dim, real t) = 0;
    virtual real phi(real fracvol, Point center, real t) = 0;
    virtual real psi(real fracvol, Point center, real t) = 0;
    virtual real nplus(real fracvol, Point center, real t) = 0;
    virtual real nminus(real fracvol, Point center, real t) = 0;
    virtual real boundary_source_term(real fracvol, int id, Point center, int dim, real t) = 0;
    virtual real boundary_phi(real fracvol, int id, Point center, real t) = 0;
    virtual real boundary_psi(real fracvol, int id, Point center, real t) = 0;
    virtual real boundary_nplus(real fracvol, int id, Point center, real t) = 0;
    virtual real boundary_nminus(real fracvol, int id, Point center, real t) = 0;
    virtual real permittivity(real fracvol, Point center, real t) = 0;
};

class HigFlowLegacyMultiphaseElectroosmotic : public HigFlowMultiphaseElectroosmoticProblem {
public:
    real (*fn_source_term)(real, Point, int, real);
    real (*fn_phi)(real, Point, real);
    real (*fn_psi)(real, Point, real);
    real (*fn_nplus)(real, Point, real);
    real (*fn_nminus)(real, Point, real);
    real (*fn_boundary_source_term)(real, int, Point, int, real);
    real (*fn_boundary_phi)(real, int, Point, real);
    real (*fn_boundary_psi)(real, int, Point, real);
    real (*fn_boundary_nplus)(real, int, Point, real);
    real (*fn_boundary_nminus)(real, int, Point, real);
    real (*fn_permittivity)(real, Point, real);

    real source_term(real f, Point c, int d, real t) { return fn_source_term(f, c, d, t); }
    real phi(real f, Point c, real t)                { return fn_phi(f, c, t); }
    real psi(real f, Point c, real t)                { return fn_psi(f, c, t); }
    real nplus(real f, Point c, real t)              { return fn_nplus(f, c, t); }
    real nminus(real f, Point c, real t)             { return fn_nminus(f, c, t); }
    real boundary_source_term(real f, int id, Point c, int d, real t) { return fn_boundary_source_term(f, id, c, d, t); }
    real boundary_phi(real f, int id, Point c, real t)    { return fn_boundary_phi(f, id, c, t); }
    real boundary_psi(real f, int id, Point c, real t)    { return fn_boundary_psi(f, id, c, t); }
    real boundary_nplus(real f, int id, Point c, real t)  { return fn_boundary_nplus(f, id, c, t); }
    real boundary_nminus(real f, int id, Point c, real t) { return fn_boundary_nminus(f, id, c, t); }
    real permittivity(real f, Point c, real t)            { return fn_permittivity(f, c, t); }
};

// --- viscoelastico com bandas de cisalhamento (modelo VCM)
class HigFlowShearBandingProblem {
public:
    virtual ~HigFlowShearBandingProblem() {}
    virtual real tensor(Point center, int i, int j, real t) = 0;
    virtual real tensor_A(Point center, int i, int j, real t) = 0;
    virtual real tensor_B(Point center, int i, int j, real t) = 0;
    virtual real nA(Point center, real t) = 0;
    virtual real nB(Point center, real t) = 0;
    virtual real cA(Point center, real t, real CAeq, real chi, real ANA) = 0;
    virtual real cB(Point center, real t, real CBeq, real chi, real ANA) = 0;
    virtual real boundary_nA(int id, Point center, real t) = 0;
    virtual real boundary_nB(int id, Point center, real t) = 0;
};
class HigFlowLegacyShearBanding : public HigFlowShearBandingProblem {
public:
    real (*fn_tensor)(Point, int, int, real);
    real (*fn_tensor_A)(Point, int, int, real);
    real (*fn_tensor_B)(Point, int, int, real);
    real (*fn_nA)(Point, real);
    real (*fn_nB)(Point, real);
    real (*fn_cA)(Point, real, real, real, real);
    real (*fn_cB)(Point, real, real, real, real);
    real (*fn_boundary_nA)(int, Point, real);
    real (*fn_boundary_nB)(int, Point, real);
    real tensor(Point c, int i, int j, real t)   { return fn_tensor(c, i, j, t); }
    real tensor_A(Point c, int i, int j, real t) { return fn_tensor_A(c, i, j, t); }
    real tensor_B(Point c, int i, int j, real t) { return fn_tensor_B(c, i, j, t); }
    real nA(Point c, real t)                     { return fn_nA(c, t); }
    real nB(Point c, real t)                     { return fn_nB(c, t); }
    real cA(Point c, real t, real e, real x, real a) { return fn_cA(c, t, e, x, a); }
    real cB(Point c, real t, real e, real x, real a) { return fn_cB(c, t, e, x, a); }
    real boundary_nA(int id, Point c, real t)    { return fn_boundary_nA(id, c, t); }
    real boundary_nB(int id, Point c, real t)    { return fn_boundary_nB(id, c, t); }
};

// --- elastoviscoplastico
class HigFlowElastoviscoplasticProblem {
public:
    virtual ~HigFlowElastoviscoplasticProblem() {}
    virtual real tensor(Point center, int i, int j, real t) = 0;
    virtual real kernel(int dim, real lambda, real tol) = 0;
    virtual real kernel_inverse(int dim, real lambda, real tol) = 0;
    virtual real kernel_jacobian(int dim, real lambda, real tol) = 0;
};
class HigFlowLegacyElastoviscoplastic : public HigFlowElastoviscoplasticProblem {
public:
    real (*fn_tensor)(Point, int, int, real);
    real (*fn_kernel)(int, real, real);
    real (*fn_kernel_inverse)(int, real, real);
    real (*fn_kernel_jacobian)(int, real, real);
    real tensor(Point c, int i, int j, real t)    { return fn_tensor(c, i, j, t); }
    real kernel(int d, real l, real tol)          { return fn_kernel(d, l, tol); }
    real kernel_inverse(int d, real l, real tol)  { return fn_kernel_inverse(d, l, tol); }
    real kernel_jacobian(int d, real l, real tol) { return fn_kernel_jacobian(d, l, tol); }
};

// --- suspensao que engrossa sob cisalhamento
class HigFlowSuspensionProblem {
public:
    virtual ~HigFlowSuspensionProblem() {}
    virtual real tensor(Point center, int i, int j, real t) = 0;
    virtual real tensor_A(Point center, int i, int j, real t) = 0;
    virtual real X(Point center, real t, real X0, real chi, real chi_J) = 0;
    virtual real vol_frac(Point center, real t) = 0;
    virtual real alpha(Point center, real t, real alpha, real phi, real phircp) = 0;
};
class HigFlowLegacySuspension : public HigFlowSuspensionProblem {
public:
    real (*fn_tensor)(Point, int, int, real);
    real (*fn_tensor_A)(Point, int, int, real);
    real (*fn_X)(Point, real, real, real, real);
    real (*fn_vol_frac)(Point, real);
    real (*fn_alpha)(Point, real, real, real, real);
    real tensor(Point c, int i, int j, real t)   { return fn_tensor(c, i, j, t); }
    real tensor_A(Point c, int i, int j, real t) { return fn_tensor_A(c, i, j, t); }
    real X(Point c, real t, real x0, real x, real xj) { return fn_X(c, t, x0, x, xj); }
    real vol_frac(Point c, real t)               { return fn_vol_frac(c, t); }
    real alpha(Point c, real t, real a, real p, real pr) { return fn_alpha(c, t, a, p, pr); }
};

#endif