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

#endif
