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

#endif