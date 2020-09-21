#include "higtree.h"
#include "higtree-io.h"
#include "higtree-iterator.h"
#include "pdomain.h"
#include "lbal.h"
#include "solver.h"
#define DEBUG
#include "Debug-c.h"

#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950288)
#endif

DECL_CLOCK(total)
DECL_CLOCK(poisson)
DECL_CLOCK(syncprop)
DECL_CLOCK(firstiter)

real cos2( real x) {

    return ( cos(x)*cos(x) );
}

real sin2( real x) {

    return ( sin(x)*sin(x) );
}


real analyticU(int dim, Point pt, real t) {

    real x = pt[0];
    real y = pt[1];
    real z = pt[2];
    real w1 = 2.0;
    real w2 = 2.0;
    real w3 = 2.0;
    real w4 = 1.0;
    real pi = M_PI;
    real ans;

    switch (dim) {
        case 0:
            ans = sin2(w1*pi*x + w2*pi*y + w3*pi*z + w4*t);
            break;
        case 1:
            ans = -cos2(w1*pi*x + w2*pi*y + w3*pi*z + w4*t);
            break;
        case 2:
            ans  = w1 * cos2(w1 * pi* x + w2 * pi* y + w3 * pi* z + w4 * t) / w3
                 + w2 * cos2(w1 * pi* x + w2 * pi* y + w3 * pi* z + w4 * t) / w3;
            break;
    }
    return (ans);
}

real analyticP(Point pt, real t) {

    real x = pt[0];
    real y = pt[1];
    real z = pt[2];
    real w1 = 2.0;
    real w2 = 2.0;
    real w3 = 2.0;
    real w4 = 1.0;
    real pi = M_PI;
    real ans;

    ans = cos(w1*pi*x + w2*pi*y + w3*pi*z + w4*t);

    return(ans);
}

sim_boundary *make_bc(hig_cell *bcg, int type) {
    mp_mapper *bm = mp_create();
    higcit_celliterator *it;
    it = higcit_create_all_leaves(bcg);
    mp_assign_from_celliterator(bm, it, 0);
    higcit_destroy(it);
    sim_boundary *bc = sb_create(bcg, type, bm);
    return bc;
}

void createBCP(psim_domain *psd, higio_amr_info *mi) {
    sim_domain *sd = psd_get_local_domain(psd);
    higcit_celliterator *it;
    for(int dim = 0; dim < DIM; dim++) {
        for(int dir = 0; dir < 2; dir++) {
            hig_cell *bcg = higio_read_bc_from_amr(mi, dim, dir);
            sim_boundary *bc = make_bc(bcg, NEUMANN);
            sd_add_boundary(sd, bc);
        }
    }
    int numbcs = sd_get_num_bcs(sd, NEUMANN);
    for (int i = 0; i < numbcs; i++) {
        sim_boundary *bc = sd_get_bc(sd, NEUMANN, i);
        mp_mapper *bm = sb_get_mapper(bc);
        for(it = sb_get_celliterator(bc); !higcit_isfinished(it); higcit_nextcell(it)) {
            hig_cell *bcell = higcit_getcell(it);
            Point bccenter;
            hig_get_center(bcell, bccenter);
            int bcgid = mp_lookup(bm, hig_get_cid(bcell));
            sb_set_value(bc, bcgid, 0.0);
        }
        higcit_destroy(it);
    }
//    Point l, h, center;
//    hig_cell *c = sd_get_higtree(sd, 0);
//    hig_get_center(c, center);
//    POINT_SUB_SCALAR(l, center, EPSDELTA);
//    POINT_ADD_SCALAR(h, center, EPSDELTA);
//    hig_cell *bcg = hig_create_root(l, h);
//    sim_boundary *bc = make_bc(bcg, DIRICHLET);
//    sd_add_boundary(sd, bc);
//    mp_mapper *bm = sb_get_mapper(bc);
//    int bcgid = mp_lookup(bm, hig_get_cid(bcg));
//    real valp = analyticP(center, 0.0);
//    sb_set_value(bc, bcgid, valp);
}

void updateBCP(psim_domain *psd, real t) {

//    sim_domain *sd = psd_get_local_domain(psd);
//    Point center;
//
//    sim_boundary *bc = sd_get_bc(sd, DIRICHLET, 0); // A numeracao eh diferente para N e D
//    hig_cell *c = sb_get_higtree(bc);
//    hig_get_center(c, center);
//    mp_mapper *bm = sb_get_mapper(bc);
//    int bcgid = mp_lookup(bm, hig_get_cid(c));
//    real valp = analyticP(center, t);
//    sb_set_value(bc, bcgid, valp);
}


void updateBCP_incremental(psim_domain *psd, real t, real dt) {

//    sim_domain *sd = psd_get_local_domain(psd);
//    Point center;
//
//    sim_boundary *bc = sd_get_bc(sd, DIRICHLET, 0); // A numeracao eh diferente para N e D
//    hig_cell *c = sb_get_higtree(bc);
//    hig_get_center(c, center);
//    mp_mapper *bm = sb_get_mapper(bc);
//    int bcgid = mp_lookup(bm, hig_get_cid(c));
//    real valp = analyticP(center, t+dt) - analyticP(center, t);
//    sb_set_value(bc, bcgid, valp);
}



void updateBCU(psim_facet_domain *psfdu[DIM],real t) {

    higcit_celliterator *it;
    sim_facet_domain *sfdu[DIM];

    for(int dim = 0; dim < DIM; dim++) {

        sfdu[dim] = psfd_get_local_domain(psfdu[dim]);
        sim_domain *sd = sfdu[dim]->cdom;
        int numbcs = sd_get_num_bcs(sd, DIRICHLET);
        //DEBUG_INSPECT(numbcs, %d);
        for (int i = 0; i < numbcs; i++) {
            sim_boundary *bc = sd_get_bc(sd, DIRICHLET, i);
            mp_mapper *bm = sb_get_mapper(bc);
            for(it = sb_get_celliterator(bc); !higcit_isfinished(it); higcit_nextcell(it)) {
                hig_cell *bcell = higcit_getcell(it);
                Point bccenter;
                hig_get_center(bcell, bccenter);
                int bcgid = mp_lookup(bm, hig_get_cid(bcell));

                real val = analyticU(dim, bccenter, t);
                sb_set_value(bc, bcgid, val);
            }
            higcit_destroy(it);
        }
    }

}


void createBCU(psim_facet_domain *psfdu[DIM], higio_amr_info *mi) {
    higcit_celliterator *it;
    sim_facet_domain *sfdu[DIM];
    int hig = 0;
    sim_boundary *bctop = NULL;
    for(int dim = 0; dim < DIM; dim++) {
        sfdu[dim] = psfd_get_local_domain(psfdu[dim]);
        for(int dim2 = 0; dim2 < DIM; dim2++) {
            for(int dir = 0; dir < 2; dir++) {
                hig_cell *bcg = higio_read_bc_from_amr(mi, dim2, dir);
                sim_boundary *bc = make_bc(bcg, DIRICHLET);
                sfd_add_boundary(sfdu[dim], bc);
                if (dim == 0 && dir == 1 && dim2 == 1) {
                    bctop = bc;
                }
            }
        }
//        sim_domain *sd = sfdu[dim]->cdom;
//        int numbcs = sd_get_num_bcs(sd, DIRICHLET);
//        DEBUG_INSPECT(numbcs, %d);
//        for (int i = 0; i < numbcs; i++) {
//            sim_boundary *bc = sd_get_bc(sd, DIRICHLET, i);
//            mp_mapper *bm = sb_get_mapper(bc);
//            for(it = sb_get_celliterator(bc); !higcit_isfinished(it); higcit_nextcell(it)) {
//                hig_cell *bcell = higcit_getcell(it);
//                Point bccenter;
//                hig_get_center(bcell, bccenter);
//                int bcgid = mp_lookup(bm, hig_get_cid(bcell));
//
//                real val = analyticU(dim, bccenter, 0.0);
//                sb_set_value(bc, bcgid, val);
//            }
//            higcit_destroy(it);
//        }
    }

    updateBCU(psfdu, 0.0);
}



void initializeP(psim_domain *psdp, distributed_property *dpp) {
    higcit_celliterator *it;
    sim_domain *sdp = psd_get_local_domain(psdp);
    mp_mapper *m = sd_get_domain_mapper(sdp);
    for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        int cgid = mp_lookup(m, hig_get_cid(c));
        Point ccenter;
        hig_get_center(c, ccenter);
        real value = analyticP(ccenter, 0.0);
        dp_set_value(dpp, cgid, value);
    }
    higcit_destroy(it);
}

real rhsU(int dim, real Re, Point pt, real t) {

    real x = pt[0];
    real y = pt[1];
    real z = pt[2];
    real w1 = 2.0;
    real w2 = 2.0;
    real w3 = 2.0;
    real w4 = 1.0;
    real pi = M_PI;
    real rhs = 0.0;

    switch (dim) {
        case 0:
//            rhs =(-4*pow(pi,2)*(pow(w1,2) + pow(w2,2) + pow(w3,2))*cos(2*(t*w4 + pi*(w1*x + w2*y + w3*z))) +
//                  2*Re*(-(pi*w1) + (pi*(w1 - w2) + 2*w4)*cos(t*w4 + pi*(w1*x + w2*y + w3*z)))
//                  *sin(t*w4 + pi*(w1*x + w2*y + w3*z)))/(2.*Re);
            rhs = (-2*pow(pi,2)*(pow(w1,2) + pow(w2,2) + pow(w3,2))*cos(2*(t*w4 + pi*(w1*x + w2*y + w3*z)))
                   + Re*(-(pi*w1) + 2*(pi*w1 + w4)*cos(t*w4 + pi*(w1*x + w2*y + w3*z)))
                   * sin(t*w4 + pi*(w1*x + w2*y + w3*z)))/Re;
            break;

        case 1:
//            rhs = -(4*pow(pi,2)*(pow(w1,2) + pow(w2,2) + pow(w3,2))*cos(2*(t*w4 + pi*(w1*x + w2*y + w3*z))) +
//                    2*Re*(pi*w2 + (pi*(-w1 + w2) - 2*w4)*cos(t*w4 + pi*(w1*x + w2*y + w3*z)))
//                    *sin(t*w4 + pi*(w1*x + w2*y + w3*z)))/(2.*Re);
            rhs = (-2*pow(pi,2)*(pow(w1,2) + pow(w2,2) + pow(w3,2))*cos(2*(t*w4 + pi*(w1*x + w2*y + w3*z)))
                   + Re*(-(pi*w2) + 2*(pi*w1 + w4)*cos(t*w4 + pi*(w1*x + w2*y + w3*z)))
                   *sin(t*w4 + pi*(w1*x + w2*y + w3*z)))/Re;
            break;

        case 2:
//            rhs = (4*pow(pi,2)*(w1 + w2)*(pow(w1,2) + pow(w2,2) + pow(w3,2))*cos(2*(t*w4 + pi*(w1*x + w2*y + w3*z))) -
//                   2*Re*(pi*pow(w3,2) + (w1 + w2)*(pi*(w1 - w2) + 2*w4)*cos(t*w4 + pi*(w1*x + w2*y + w3*z)))
//                   *sin(t*w4 + pi*(w1*x + w2*y + w3*z)))/(2.*Re*w3);
            rhs = (2*pow(pi,2)*(w1 + w2)*(pow(w1,2) + pow(w2,2) + pow(w3,2))
                   * cos(2*(t*w4 + pi*(w1*x + w2*y + w3*z))) - Re*(pi*pow(w3,2)+ 2*(w1 + w2)*(pi*w1 + w4)
                   * cos(t*w4 + pi*(w1*x + w2*y + w3*z)))*sin(t*w4 + pi*(w1*x + w2*y + w3*z)))/(Re*w3);
            break;
    }

    return( rhs );

}


real rhsU_millena(int dim, real Re, Point pt, real t) {

    // Por algum motivo Re nao eh usado abaixo

    real x = pt[0];
    real y = pt[1];
    real z = pt[2];
    real w1 = 2.0;
    real w2 = 2.0;
    real w3 = 2.0;
    real w4 = 1.0;
    real pi = M_PI;
    real ro = 1.0;
    real k1 = 1.0;
    real k2 = 0.0;
    real rhs = 0.0;


    switch (dim) {
        case 0:
            rhs = - sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w1 * pi
            - 4.0 * w1 * w1 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * pi * pi * k2
            * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
            + 2.0 * w1 * w1 * pi * pi * k2 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
            + 2.0 * w2 * w2 * pi * pi * k2 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
            + 2.0 * pi * pi * k2 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w3 * w3
            + 2.0 * ro * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
            * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w4
            - 4.0 * w1 * w1 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * pi * pi * k1
            + 4.0 * k2 * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
            * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w2 * w2 * pi
            * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
            * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * pi
            - 4.0 * w2 * w2 * pi * pi * k1 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
            - 4.0 * pi * pi * k1 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w3 * w3
            - 4.0 * w2 * w2 * pi * pi * k2 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
            * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
            + 4.0 * k2 * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
            * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * pi
            * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
            * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * pi * w3 * w3
            - 4.0 * pi * pi * k2 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
            * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w3 * w3
            + 2.0 * ro * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
            * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * pi * w1
            + 4.0 * k2 * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
            * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w1 * w1 * pi
            * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
            * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * pi
            + 2.0 * w1 * w1 * pi * pi * k1
            + 2.0 * w2 * w2 * pi * pi * k1
            + 2.0 * pi * pi * k1 * w3 * w3;
            break;
        case 1:
            rhs = 2.0 * ro * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
            * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w4
            + 4.0 * k2 * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
            * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w2 * w2 * pi
            * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
            * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * pi
            - 4.0 * w1 * w1 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * pi * pi * k2
            * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
            + 4.0 * k2 * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
            * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * pi
            * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
            * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * pi * w3 * w3
            - 4.0 * pi * pi * k2 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi*z + w4 * t)
            * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w3 * w3
            + 2.0 * ro * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
            * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * pi * w1
            - 4.0 * w2 * w2 * pi * pi * k1 * cos2(w1 * pi* x + w2 * pi * y + w3 * pi * z + w4 * t)
            - 4.0 * w1 * w1 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * pi * pi * k1
            - 4.0 * pi * pi * k1 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w3 * w3
            - sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w2 * pi
            - 4.0 * w2 * w2 * pi * pi * k2 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
            * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
            + 4.0 * k2 * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
            * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w1 * w1 * pi
            * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
            * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * pi
            + 2.0 * w1 * w1 * pi * pi * k2 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
            + 2.0 * pi * pi * k2 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w3 * w3
            + 2.0 * w2 * w2 * pi * pi * k2 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
            + 2.0 * pi * pi * k1 * w3 * w3
            + 2.0 * w2 * w2 * pi * pi * k1
            + 2.0 * w1 * w1* pi * pi * k1;
            break;
        case 2:
            rhs = (-sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w3 * w3 * pi
                   - 2.0 * pi * pi * w2 * w2 * w2 * k2 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
                   - 2.0 * w3 * w3 * pi * pi * k2 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w1
                   - 2.0 * w3 * w3 * pi * pi * k2 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w2
                   - 2.0 * pi * pi * w1 * w1 * k2 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w2
                   - 2.0 * ro * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
                   * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * pi * w1 * w1
                   - 2.0 * ro * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
                   * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * pi * w1 * w2
                   - 2.0 * pi * pi * w2 * w2 * k2 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w1
                   - 2.0 * pi * pi * w2 * w2 * k1 * w1 - 2.0 * pi * pi * w1 * w1 * w1 * k1
                   + 4.0 * pi * pi * w2 * w2 * w2 * k1 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
                   + 4.0 * pi * pi * w1 * w1 * w1 * k1 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
                   - 4.0 * k2 * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
                   * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w3 * w3 * pi
                   * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
                   * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * pi * w1
                   - 4.0 * k2 * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
                   * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w3 * w3 * pi
                   * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
                   * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * pi * w2
                   + 4.0 * w3 * w3 * pi * pi * k1 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w1
                   + 4.0 * w3 * w3 * pi * pi * k1 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w2
                   - 2.0 * ro * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
                   * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w4 * w1
                   - 2.0 * ro * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
                   * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w4 * w2
                   + 4.0 * pi * pi * w1 * w1 * k1 * w2 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
                   + 4.0 * pi * pi * w1 * w1 * w1 * k2 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
                   * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
                   + 4.0 * w3 * w3 * pi * pi * k2 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
                   * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w2
                   + 4.0 * w3 * w3 * pi * pi * k2 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
                   * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w1
                   - 4.0 * k2 * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
                   * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w2 * w2 * pi
                   * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
                   * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * pi * w1
                   + 4.0 * pi * pi * w1 * w1 * k2 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w2
                   * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
                   - 4.0 * k2 * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
                   * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w1 * w1 * w1 * pi
                   * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
                   * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * pi
                   - 4.0 * k2 * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
                   * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w1 * w1 * pi
                   * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
                   * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * pi * w2
                   + 4.0 * pi * pi * w2 * w2 * k2 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
                   * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w1
                   + 4.0 * pi * pi * w2 * w2 * w2 * k2 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
                   * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
                   - 4.0 * k2 * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
                   * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w2 * w2 * w2 * pi
                   * sin(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
                   * cos(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * pi
                   + 4.0 * pi * pi * w2 * w2 * k1 * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t) * w1
                   - 2.0 * pi * pi * w2 * w2 * w2 * k1
                   - 2.0 * w3 * w3 * pi * pi * k1 * w1
                   - 2.0 * pi * pi * w1 * w1 * k1 * w2
                   - 2.0 * pi * pi * w1 * w1 * w1 * k2
                   * cos2(w1 * pi * x + w2 * pi * y + w3 * pi * z + w4 * t)
                   - 2.0 * w3 * w3 * pi * pi * k1 * w2 ) / w3;
            break;
    }
    return ( rhs );

}

//real analyticU(real x, real y) {
//    return 0.0;
//    return 8.0*(x*x*(1.0+x*(-2.0+x)))*(y*(-2.0+4.0*y*y));
//}
//
//real analyticV(real x, real y) {
//    return 0.0;
//    return -8.0*(x*(2.0+x*(-6.0+4.0*x)))*(y*y*(-1.0+y*y));;
//}

void initializeU(psim_facet_domain *psfdu[DIM], distributed_property *dpu[DIM]) {
    higfit_facetiterator *fit;
    sim_facet_domain *sfdu[DIM];
    for(int dim = 0; dim < DIM; dim++) {
        sfdu[dim] = psfd_get_local_domain(psfdu[dim]);
        mp_mapper *m = sfd_get_domain_mapper(sfdu[dim]);
        Point l, h, delta;
        hig_cell *dom = sfd_get_higtree(sfdu[dim], 0);
        hig_get_lowpoint(dom, l);
        hig_get_highpoint(dom, h);
        hig_get_delta(dom, delta);
        for(fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            hig_facet *f = higfit_getfacet(fit);
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            Point x;
            POINT_SUB(x, fcenter, l);
            POINT_DIV(x, x, delta);
            int fgid = mp_lookup(m, hig_get_fid(f));
            real val;
            val = analyticU(dim, x, 0);
            dp_set_value(dpu[dim], fgid, val);
        }
        higfit_destroy(fit);
    }
}

void compute_facet_dudx_and_du2dx2_at_point(sim_facet_domain *sfdu, Point center, Point delta, int dim, real alpha, distributed_property *dpu, sim_stencil *stn, real valueatpoint, real *dudx, real *du2dx2) {

    //DEBUG_PASS;
    Point p;
    POINT_ASSIGN(p, center);
    p[dim] = center[dim] - alpha * delta[dim];
    //DEBUG_PASS;
    real valuel = sfd_dp_interpolate(sfdu, dpu, center, p, stn);
    p[dim] = center[dim] + alpha * delta[dim];
    //DEBUG_PASS;
    real valueh = sfd_dp_interpolate(sfdu, dpu, center, p, stn);
    //DEBUG_PASS;
    if (dudx != NULL) {
        *dudx = (valueh - valuel) / (2.0 * alpha * delta[dim]);
    }
    //DEBUG_PASS;
    if (du2dx2 != NULL) {
        *du2dx2 = (valuel - 2.0*valueatpoint + valueh)/(delta[dim] * delta[dim] * alpha * alpha);
    }
    //DEBUG_PASS;
}

void compute_dpdx_at_point(sim_domain *sdp, Point center, Point delta, int dim, real alpha, distributed_property *dpp, sim_stencil *stn, real *dpdx) {
    Point p;
    POINT_ASSIGN(p, center);
    p[dim] = center[dim] - alpha * delta[dim];
    real valuel = sd_dp_interpolate(sdp, dpp, center, p, stn);
    p[dim] = center[dim] + alpha * delta[dim];
    real valueh = sd_dp_interpolate(sdp, dpp, center, p, stn);
    if (dpdx != NULL) {
        *dpdx = (valueh - valuel) / (2.0 * alpha * delta[dim]);
    }
}

void compute_ustar(psim_facet_domain *psfdu[DIM], int dim, distributed_property *dpu[DIM], real t, real dt, real Re, distributed_property *dpustar[DIM], sim_stencil *stn) {

    higfit_facetiterator *fit;
    sim_facet_domain *sfdu[DIM];
    // Note que o uso de dim esta errado aqui
    for(int dim = 0; dim < DIM; dim++) {
        sfdu[dim] = psfd_get_local_domain(psfdu[dim]);
    }
    mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);

    for(fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
        hig_facet *f = higfit_getfacet(fit);
        int fgid = mp_lookup(mu, hig_get_fid(f));
        Point fcenter;
        Point fdelta;
        hig_get_facet_center(f, fcenter);
        hig_get_facet_delta(f, fdelta);

        real valuec = dp_get_value(dpu[dim], fgid);
        real tot = rhsU(dim, Re, fcenter, t);
        for (int dim2 = 0; dim2 < DIM; dim2++) {
            real du2dx2 = 0.0;
            real dudx = 0.0;
            compute_facet_dudx_and_du2dx2_at_point(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], stn, valuec, &dudx, &du2dx2);
            real u;
            // Bug when value at center
            u = sfd_dp_interpolate(sfdu[dim2], dpu[dim2], fcenter, fcenter, stn);
            tot += - u * dudx + du2dx2 / Re;
        }
        real ustar = valuec + dt * tot;
        dp_set_value(dpustar[dim], fgid, ustar);
    }
    higfit_destroy(fit);
    START_CLOCK(syncprop);
    dp_sync(dpustar[dim]);
    STOP_CLOCK(syncprop);
}

void compute_ustar_incremental(psim_domain *psdp, psim_facet_domain *psfdu[DIM],  int dim, distributed_property *dpp,distributed_property *dpu[DIM], real t, real dt, real Re, distributed_property *dpustar[DIM], sim_stencil *stn) {

    higfit_facetiterator *fit;

    sim_domain *sdp = psd_get_local_domain(psdp);
    sim_facet_domain *sfdu[DIM];

    for(int dim2 = 0; dim2 < DIM; dim2++) {
        sfdu[dim2] = psfd_get_local_domain(psfdu[dim2]);
    }

    mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);

    for(fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
        hig_facet *f = higfit_getfacet(fit);
        int fgid = mp_lookup(mu, hig_get_fid(f));
        Point fcenter;
        Point fdelta;

        hig_get_facet_center(f, fcenter);
        hig_get_facet_delta(f, fdelta);

        real valuec = dp_get_value(dpu[dim], fgid);
        real tot = 0.0;
        real dpdx;
        compute_dpdx_at_point(sdp, fcenter, fdelta, dim, 0.5, dpp, stn, &dpdx);
        tot += -dpdx;
        tot += rhsU(dim, Re, fcenter, t);

        for (int dim2 = 0; dim2 < DIM; dim2++) {
            //DEBUG_INSPECT(dim2, %d);
            real du2dx2 = 0.0;
            real dudx = 0.0;
            //DEBUG_PASS;
            compute_facet_dudx_and_du2dx2_at_point(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], stn, valuec, &dudx, &du2dx2);
            real u;
            // Bug when value at center
            //DEBUG_PASS;
            u = sfd_dp_interpolate(sfdu[dim2], dpu[dim2], fcenter, fcenter, stn);
            tot += - u * dudx + du2dx2 / Re;
            //DEBUG_PASS;
        }
        //DEBUG_PASS;
        real ustar = valuec + dt * tot;
        dp_set_value(dpustar[dim], fgid, ustar);
    }
    higfit_destroy(fit);
    START_CLOCK(syncprop);
    dp_sync(dpustar[dim]);
    STOP_CLOCK(syncprop);
}

void remove_pressure_singularity(psim_domain *psdp, distributed_property *dpn, real t, real dt, int is_incr, solver *slvp) {

    int cgid;
    real value;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Obtem cgid se for rank == 0
    if (rank == 0) {

        sim_domain *sdp = psd_get_local_domain(psdp);
        mp_mapper *mp = sd_get_domain_mapper(sdp);
        higcit_celliterator *it;
        Point ccenter;

        // Toma a primeira celula que aparecer
        it = sd_get_domain_celliterator(sdp);
        hig_cell *c = higcit_getcell(it);
        hig_get_center(c, ccenter);
        cgid = mp_lookup(mp, hig_get_cid(c));
        higcit_destroy(it);

        // Calcula valor do lado direito
        if (is_incr)
            // valor esta errado
            //value = analyticP(ccenter,t+dt) - analyticP(ccenter,t);
            // o correto eh
            value = analyticP(ccenter,t+dt) - dp_get_value(dpn, cgid);
        else
            value = analyticP(ccenter,t+dt);
    }

    MPI_Bcast(&cgid,  1, MPI_INT,    0, MPI_COMM_WORLD);
    MPI_Bcast(&value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    slv_impose_value(slvp, cgid, value);

}



void compute_p(psim_domain *psdp, psim_facet_domain *psfdu[DIM], distributed_property *dpp, distributed_property *dpn, distributed_property *dpustar[DIM], real t, real dt, real Re, int is_incr, sim_stencil *stn) {
    static solver *slvp = NULL;
    if (slvp == NULL) {
        int localdomainsize = psd_get_local_domain_size(psdp);
        slvp = slv_create(SOLVER_ANY, psd_get_first_id(psdp), localdomainsize);
        slv_set_maxnonzeros(slvp, 300);
    }

    higcit_celliterator *it;
    sim_domain *sdp = psd_get_local_domain(psdp);
    sim_facet_domain *sfdu[DIM];
    for(int dim = 0; dim < DIM; dim++) {
        sfdu[dim] = psfd_get_local_domain(psfdu[dim]);
    }
    mp_mapper *mp = sd_get_domain_mapper(sdp);
    for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        Point ccenter;
        Point cdelta;
        hig_get_center(c, ccenter);
        hig_get_delta(c, cdelta);

        real sumdudx = 0.0;
        for(int dim = 0; dim < DIM; dim++) {
            real dudx;
            compute_facet_dudx_and_du2dx2_at_point(sfdu[dim], ccenter, cdelta, dim, 0.5, dpustar[dim], stn, 0.0, &dudx, NULL);
            sumdudx += dudx;
        }
        stn_reset(stn);
        stn_set_rhs(stn, sumdudx / dt);
        real alpha = 0.0;
        for(int dim = 0; dim < DIM; dim++) {
            real w = 1.0/(cdelta[dim]*cdelta[dim]);
            alpha -= 2.0 * w;
            Point p;
            POINT_ASSIGN(p, ccenter);
            p[dim] = ccenter[dim] + cdelta[dim];
            sd_get_stencil(sdp, ccenter, p, w, stn);
            p[dim] = ccenter[dim] - cdelta[dim];
            sd_get_stencil(sdp, ccenter, p, w, stn);
        }
        sd_get_stencil(sdp, ccenter, ccenter, alpha, stn);

        int *ids   = psd_stn_get_gids(psdp, stn);
        real *vals = stn_get_vals(stn);
        int numelems = stn_get_numelems(stn);

        int cgid = psd_get_global_id(psdp, c);

        slv_set_bi(slvp, cgid, stn_get_rhs(stn));
        slv_set_Ai(slvp, cgid, numelems, ids, vals);
    }

    higcit_destroy(it);

    remove_pressure_singularity(psdp, dpn, t, dt, is_incr, slvp);
    slv_assemble(slvp);

    //DEBUG_DIFF_TIME;
    START_CLOCK(poisson);
    slv_solve(slvp);
    STOP_CLOCK(poisson);
    //DEBUG_DIFF_TIME;

    dp_slv_load_from_solver(dpp, slvp);
}

void compute_ut(psim_domain *psdp, psim_facet_domain *psfdu[DIM], int dim, distributed_property *dpp, distributed_property *dpustar[DIM], distributed_property *dpu[dim], real dt, real Re, sim_stencil *stn) {
    sim_domain *sdp = psd_get_local_domain(psdp);
    sim_facet_domain *sfdu[DIM];
    higfit_facetiterator *fit;
    for(int dim = 0; dim < DIM; dim++) {
        sfdu[dim] = psfd_get_local_domain(psfdu[dim]);
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        for(fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            hig_facet *f = higfit_getfacet(fit);
            int fgid = mp_lookup(mu, hig_get_fid(f));
            Point fcenter;
            Point fdelta;
            hig_get_facet_center(f, fcenter);
            hig_get_facet_delta(f, fdelta);

            real dpdx;
            compute_dpdx_at_point(sdp, fcenter, fdelta, dim, 0.5, dpp, stn, &dpdx);
            real ustar = dp_get_value(dpustar[dim], fgid);
            real utdt = ustar - dt * dpdx;
            dp_set_value(dpu[dim], fgid, utdt);
        }
        higfit_destroy(fit);
        START_CLOCK(syncprop);
        dp_sync(dpu[dim]);
        STOP_CLOCK(syncprop);
    }
}

void compute_pt(psim_domain *psdp, distributed_property *ddeltap, distributed_property *dpp, real dt, sim_stencil *stn) {

    sim_domain *sdp = psd_get_local_domain(psdp);

    higcit_celliterator *it;
    mp_mapper *mp = sd_get_domain_mapper(sdp);

    for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);

        int cgid = mp_lookup(mp, hig_get_cid(c));
        real p      = dp_get_value(dpp, cgid);
        real deltap = dp_get_value(ddeltap, cgid);
        real newp   = p + deltap;
        dp_set_value(dpp, cgid, p + deltap);

    }
    higcit_destroy(it);
    START_CLOCK(syncprop);
    dp_sync(dpp);
    STOP_CLOCK(syncprop);
}

void updatePandU(psim_domain *psdp, psim_facet_domain *psfdu[DIM], distributed_property *dpp, distributed_property *dpu[DIM], distributed_property *dpustar[DIM], real tn, real dt, real Re,sim_stencil *stn) {

    int is_incr = 0;
    higfit_facetiterator *fit;
    //DEBUG_DIFF_TIME;

    updateBCU(psfdu, tn+dt);

    for(int dim = 0; dim < DIM; dim++) {
        compute_ustar(psfdu, dim, dpu, tn, dt, Re, dpustar, stn);
    }
    //updateBCP(psdp, t);
    //DEBUG_DIFF_TIME;
    compute_p(psdp, psfdu, dpp, dpp, dpustar, tn, dt, Re, is_incr, stn);
    //DEBUG_DIFF_TIME;

    for(int dim = 0; dim < DIM; dim++) {
        compute_ut(psdp, psfdu, dim, dpp, dpustar, dpu, dt, Re, stn);
    }
    //DEBUG_DIFF_TIME;
}

void updatePandU_incremental(psim_domain *psdp, psim_facet_domain *psfdu[DIM], distributed_property *dpp, distributed_property *ddeltap, distributed_property *dpu[DIM], distributed_property *dpustar[DIM], real tn, real dt, real Re, sim_stencil *stn) {

    int is_incr = 1;
    higfit_facetiterator *fit;
    //DEBUG_DIFF_TIME;


    // Segundo Guermond tem que ser no tempo n+1
    updateBCU(psfdu, tn+dt);

    for(int dim = 0; dim < DIM; dim++) {
        compute_ustar_incremental(psdp, psfdu, dim, dpp, dpu, tn, dt, Re, dpustar, stn);
    }
    //updateBCP_incremental(psdp, t, dt);
    //DEBUG_DIFF_TIME;
    compute_p(psdp, psfdu, ddeltap, dpp, dpustar, tn, dt, Re, is_incr, stn);
    //DEBUG_DIFF_TIME;

    for(int dim = 0; dim < DIM; dim++) {
        compute_ut(psdp, psfdu, dim, ddeltap, dpustar, dpu, dt, Re, stn);
        // Tem que tomar cuidado com o gradiente de pressao nas faces do contorno e com a condicao de contorno para velocidade
    }
    compute_pt(psdp, ddeltap, dpp, dt, stn);
    //DEBUG_DIFF_TIME;

}

void printResults(psim_domain *psdp, psim_facet_domain *psfdu[DIM], distributed_property *dpp, distributed_property *dpu[DIM]) {
    sim_domain *sdp = psd_get_local_domain(psdp);
    sim_stencil *stn = stn_create();
    higcit_celliterator *it;
    Point l, h;
    int hig = 0;

    sim_facet_domain *sfdu[DIM];
    for(int dim = 0; dim < DIM; dim++) {
        sfdu[dim] = psfd_get_local_domain(psfdu[dim]);
    }

    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        Point ccenter;
        Point cdelta;
        hig_get_center(c, ccenter);
        hig_get_delta(c, cdelta);
        real value[DIM];
        for(int dim = 0; dim < DIM; dim++) {
            value[dim] = sfd_dp_interpolate(sfdu[dim], dpu[dim], ccenter, ccenter, stn);
        }
        printf("v(");
        for(int dim = 0; dim < DIM; dim++) {
            printf("%f ", ccenter[dim]);
        }
        printf(") = ");
        for(int dim = 0; dim < DIM; dim++) {
            printf("%0.15f ", value[dim]);
        }
        printf("\n");
    }
    stn_destroy(stn);
}

void compute_error_norm_velocity(int dim, psim_facet_domain *psfdu[DIM], distributed_property *dpu[DIM], real t, real *norm2, real *norminf) {

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    real local_maxu = -1.0;
    real local_sumu = 0.0;
    real global_maxu, global_sumu;

    sim_facet_domain *sfdu;
    sfdu = psfd_get_local_domain(psfdu[dim]);

    higfit_facetiterator *fit;
    mp_mapper *mu = sfd_get_domain_mapper(sfdu);

    for(fit = sfd_get_domain_facetiterator(sfdu); !higfit_isfinished(fit); higfit_nextfacet(fit)) {

        hig_facet *f = higfit_getfacet(fit);
        int fgid = mp_lookup(mu, hig_get_fid(f));
        Point fcenter;
        Point fdelta;
        hig_get_facet_center(f, fcenter);
        hig_get_facet_delta(f, fdelta);

        real unm = dp_get_value(dpu[dim], fgid);
        real uex = analyticU(dim, fcenter, t);

        real umod_inf = fabs( unm - uex );
        real umod_2   = umod_inf * umod_inf;

        for (int i = 0; i < DIM; i++)
            umod_2 *= fdelta[i]; // definicao de norma 2 de grid function do Leveque

        if (umod_inf > local_maxu) local_maxu = umod_inf;
        local_sumu += umod_2;
    }

    MPI_Reduce(&local_sumu, &global_sumu, 1, MPI_HIGREAL, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_maxu, &global_maxu, 1, MPI_HIGREAL, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        *norminf = global_maxu;
        *norm2   = sqrt( global_sumu );
    }

    higfit_destroy(fit);
}

void compute_error_norm_pressure(psim_domain *psdp, distributed_property *dpp, real t, real *norm2, real *norminf) {

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    real local_maxp = -1.0;
    real local_sump =  0.0;
    real global_sump,global_maxp;

    sim_domain *sdp = psd_get_local_domain(psdp);

    higcit_celliterator *it;
    mp_mapper *mp = sd_get_domain_mapper(sdp);

    for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        Point ccenter;
        Point cdelta;
        hig_get_center(c, ccenter);
        hig_get_delta(c, cdelta);

        int cgid = mp_lookup(mp, hig_get_cid(c));
        real pnm = dp_get_value(dpp, cgid);
        real pex = analyticP(ccenter, t);

        real pmod_inf = fabs( pnm - pex );
        real pmod_2   = pmod_inf * pmod_inf;

        for (int i = 0; i < DIM; i++)
            pmod_2 *= cdelta[i]; // definicao de norma 2 de grid function do Leveque

        if (pmod_inf > local_maxp) local_maxp = pmod_inf;
        local_sump += pmod_2;
    }

    // Agora tem que mandar tudo para o processador mestre e fazer a conta lah

    MPI_Reduce(&local_sump, &global_sump, 1, MPI_HIGREAL, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_maxp, &global_maxp, 1, MPI_HIGREAL, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {

        *norminf = global_maxp;
        *norm2 = sqrt( global_sump );
    }

    higcit_destroy(it);

}


void compute_and_print_error_norm_velocity(psim_facet_domain *psfdu[DIM], distributed_property *dpu[DIM], real t, int myrank) {

    real norm2;
    real norminf;

    if (myrank == 0) {
        printf("  Velocity error norms: --------------\n");
        printf("   L_inf        |   L_2 \n");
    }

    for (int dim = 0; dim < DIM; dim++) {

        compute_error_norm_velocity(dim, psfdu, dpu, t, &norm2, &norminf);
        if (myrank == 0) {
            printf("%d: %e | %e\n", dim, norminf, norm2);
        }
    }

}

void compute_and_print_error_norm_pressure(psim_domain *psdp, distributed_property *dpp, real t, int myrank) {

    real norm2;
    real norminf;

    if (myrank == 0) {
        printf("  Pressure error norms: --------------\n");
        printf("   L_inf        |   L_2 \n");
    }

    compute_error_norm_pressure(psdp, dpp, t, &norm2, &norminf);
    if (myrank == 0) {
        printf("   %e | %e\n", norminf, norm2);
    }
}

//void computeextraload(psim_domain *psdp, psim_facet_domain *psfdu[DIM], distributed_property *dpp, distributed_property *dpu[DIM], int extraload) {
//    sim_domain *sdp = psd_get_local_domain(psdp);
//
//    mp_mapper *mp = sd_get_domain_mapper(sdp);
//    DECL_AND_ALLOC(Point, rndpts, extraload);
//    DECL_AND_ALLOC(int, closestpts, extraload);
//    higcit_celliterator *it;
//
//    for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
//        hig_cell *c = higcit_getcell(it);
//        Point lp, hp;
//        hig_get_lowpoint(c, lp);
//        hig_get_highpoint(c, hp);
//        Point x;
//        for (int i = 0; i < extraload; i++) {
//            randpoint(lp, hp, x);
//            POINT_ASSIGN(rndpts[i], x);
//        }
//        real maxpossibledist = co_distance(lp, hp);
//        for (int i = 0; i < extraload; i++) {
//            real shortestdist = maxpossibledist;
//            int closestpoint = i;
//            for (int j = 0; j < extraload; j++) {
//                if (i != j) {
//                    real dist = co_distance(rndpts[i], rndpts[j]);
//                    if (dist < shortestdist) {
//                        shortestdist = dist;
//                        closestpoint = j;
//                    }
//                }
//            }
//            closestpts[i] = closestpoint;
//        }
//    }
//    free(rndpts);
//    free(closestpts);
//}


int main (int argc, char *argv[]) {
    int myrank;
    int ntasks;
    higtree_initialize(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    START_CLOCK(total);
    partition_graph *pg = pg_create(MPI_COMM_WORLD);
    pg_set_fringe_size(pg, 5);

    //DEBUG_DIFF_TIME;
    char *amrfilename = argv[1];
    real dt = atof(argv[2]);
    real Re = atof(argv[3]);
    int numsteps = atoi(argv[4]);
    int incremental = atoi(argv[5]);
    //int extraload = atoi(argv[5]);
    //DEBUG_INSPECT(numsteps, %d);

    FILE *fd = fopen(amrfilename, "r");
    higio_amr_info *mi = higio_read_amr_info(fd);
    fclose(fd);
    DEBUG_PUSH(myrank == 0);
    load_balancer *lb = higio_partition_from_amr_info(pg, mi);
    DEBUG_POP;

    //DEBUG_DIFF_TIME
    sim_domain *sdp = sd_create(NULL);
    sd_set_interpolator_order(sdp, 3);
    unsigned numtrees = lb_get_num_local_trees(lb);
    for(unsigned i = 0; i < numtrees; ++i) {
	    sd_add_higtree(sdp, lb_get_local_tree(lb, i, NULL));
    }
    psim_domain *psdp = psd_create(sdp, pg);
    psd_synced_mapper(psdp);
    distributed_property * dpp     = psd_create_property(psdp);
    distributed_property * ddeltap = psd_create_property(psdp);

    sim_facet_domain  *sfdu[DIM];
    psim_facet_domain *psfdu[DIM];
    distributed_property *dpu[DIM];
    distributed_property *dpustar[DIM];
    for(int dim = 0; dim < DIM; dim++) {
        sfdu[dim] = sfd_create(NULL, dim);
        sfd_set_interpolator_order(sfdu[dim], 3);
        sfd_copy_higtrees_from_center_domain(sfdu[dim], sdp);
        sfd_adjust_facet_ids(sfdu[dim]);
        psfdu[dim] = psfd_create(sfdu[dim], psdp);
    }

    //DEBUG_DIFF_TIME

    createBCP(psdp, mi);
    createBCU(psfdu, mi);

    for(int dim = 0; dim < DIM; dim++) {
        psfd_compute_sfbi(psfdu[dim]);
        psfd_synced_mapper(psfdu[dim]);
        dpu[dim] = psfd_create_property(psfdu[dim]);
        dpustar[dim] = psfd_create_property(psfdu[dim]);
    }

    initializeP(psdp, dpp);
    initializeP(psdp, ddeltap);
    initializeU(psfdu, dpu);
    sim_stencil *stn = stn_create();

    if (incremental) {
        for(int step = 0; step < numsteps; step++) {

            if (myrank == 0) printf("****** \n  Step: %d | Results at t = %g \n", step, (step+1)*dt);
            if (step == 0) { START_CLOCK(firstiter); }
            updatePandU_incremental(psdp, psfdu, dpp, ddeltap, dpu, dpustar, step*dt, dt, Re, stn);
            compute_and_print_error_norm_velocity(psfdu, dpu, (step+1)*dt, myrank);
            compute_and_print_error_norm_pressure(psdp, dpp, (step+1)*dt, myrank);
            if (step == 0) { STOP_CLOCK(firstiter); }
            if (myrank == 0) printf("****** \n\n\n");
            //DEBUG_INSPECT(step, %d);
        }
    }
    else {
        for(int step = 0; step < numsteps; step++) {

            if (myrank == 0) printf("****** \n  Step: %d | Results at t = %g \n", step, (step+1)*dt);
            if (step == 0) { START_CLOCK(firstiter); }
            updatePandU(psdp, psfdu, dpp, dpu, dpustar, step*dt, dt, Re, stn);
            compute_and_print_error_norm_velocity(psfdu, dpu, (step+1)*dt, myrank);
            compute_and_print_error_norm_pressure(psdp, dpp, (step+1)*dt, myrank);
            if (step == 0) { STOP_CLOCK(firstiter); }
            if (myrank == 0) printf("****** \n\n\n");
        }
    }


//    real norm2;
//    real norminf;
//
//    for (int dim = 0; dim < DIM; dim++) {
//
//        compute_error_norm_velocity(dim, psfdu, dpu, numsteps * dt, &norm2, &norminf);
//        if (myrank == 0) {
//            DEBUG_INSPECT(dim, %d);
//            DEBUG_INSPECT(norm2, %e);
//            DEBUG_INSPECT(norminf, %e);
//        }
//    }
//
//    compute_error_norm_pressure(psdp, dpp, numsteps * dt, &norm2, &norminf);
//    DEBUG_PASS;
//    if (myrank == 0) {
//        DEBUG_INSPECT(norm2, %e);
//        DEBUG_INSPECT(norminf, %e);
//    }

    for(int dim = 0; dim < DIM; dim++) {
        dp_destroy(dpustar[dim]);
    }
    stn_destroy(stn);
    STOP_CLOCK(total);
    //DEBUG_DIFF_TIME;

    if (myrank==0) {
        DEBUG_INSPECT(GET_NSEC_CLOCK(total)/1.0e9, %lf);
        DEBUG_INSPECT(GET_NSEC_CLOCK(poisson)/1.0e9, %lf);
        DEBUG_INSPECT(GET_NSEC_CLOCK(syncprop)/1.0e9, %lf);
        DEBUG_INSPECT(GET_NSEC_CLOCK(firstiter)/1.0e9, %lf);
    }
    return 0;
}
