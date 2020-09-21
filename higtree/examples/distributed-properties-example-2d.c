#include "higtree.h"
#include "higtree-io.h"
#include "higtree-iterator.h"
#include "pdomain.h"
#include "lbal.h"
#define DEBUG
#include "Debug-c.h"

int main (int argc, char *argv[]) {
	higtree_initialize(&argc, &argv);

	int myrank;
	int ntasks;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

	partition_graph *pg = pg_create(MPI_COMM_WORLD);
	pg_set_fringe_size(pg, 2);

	DEBUG_DIFF_TIME;
	FILE *fd = fopen(argv[1], "r");
	higio_amr_info *mi = higio_read_amr_info(fd);
	fclose(fd);
	load_balancer *lb = higio_partition_from_amr_info(pg, mi);
	higcit_celliterator *it;
	higfit_facetiterator *fit;
	mp_mapper *mp = mp_create();
	mp_mapper *mu = mp_create();
	mp_mapper *mv = mp_create();
	sim_domain        *sdp = sd_create(mp);
	sim_facet_domain *sfdu = sfd_create(mu, 0);
	sim_facet_domain *sfdv = sfd_create(mv, 1);
	for(unsigned i = 0; i < lb_get_num_local_trees(lb); ++i) {
		sd_add_higtree(sdp, lb_get_local_tree(lb, i, NULL));
	}
	lb_destroy(lb);
	sfd_copy_higtrees_from_center_domain(sfdv, sdp);
	sfd_copy_higtrees_from_center_domain(sfdu, sdp);
	sfd_adjust_facet_ids(sfdv);
	sfd_adjust_facet_ids(sfdu);
	psim_domain        *psdp = psd_create(sdp, pg);
	psim_facet_domain *psfdu = psfd_create(sfdu, psdp);
	psim_facet_domain *psfdv = psfd_create(sfdv, psdp);
	DEBUG_DIFF_TIME
	psd_synced_mapper(psdp);
	psfd_compute_sfbi(psfdu);
	psfd_synced_mapper(psfdu);
	psfd_compute_sfbi(psfdv);
	psfd_synced_mapper(psfdv);
	int localsizep = psd_get_local_domain_size(psdp);
	int localsizeu = psfd_get_local_domain_size(psfdu);
	int localsizev = psfd_get_local_domain_size(psfdv);
	distributed_property * dpp = psd_create_property(psdp);
	distributed_property * dpu = psfd_create_property(psfdu);
	distributed_property * dpv = psfd_create_property(psfdv);

	real cnt = 1.0;
	for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		int cgid = mp_lookup(mp, hig_get_cid(c));
		dp_set_value(dpp, cgid, cgid + 5000.0);
	}
	higcit_destroy(it);
	for(fit = sfd_get_domain_facetiterator(sfdu); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
		hig_facet *f = higfit_getfacet(fit);
		int fgid = mp_lookup(mu, hig_get_fid(f));
		dp_set_value(dpu, fgid, fgid + 1000.0);
	}
	higfit_destroy(fit);
	for(fit = sfd_get_domain_facetiterator(sfdv); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
		hig_facet *f = higfit_getfacet(fit);
		int fgid = mp_lookup(mv, hig_get_fid(f));
		dp_set_value(dpv, fgid, fgid + 2000.0);
	}
	higfit_destroy(fit);
	for(fit = higfit_create_allfacets(sd_get_domain_celliterator(sdp), sfdu->dimofinterest); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
		hig_facet *f = higfit_getfacet(fit);
		Point fcenter;
		hig_get_facet_center(f, fcenter);
		hig_cell *c = f->c;
		Point ccenter;
		hig_get_center(c, ccenter);
		int fgid = mp_lookup(mu, hig_get_fid(f));
		printf("%d - %d - %d (%f, %f) (%f, %f)\n", myrank, hig_get_fid(f), fgid, ccenter[0], ccenter[1], fcenter[0], fcenter[1]);
	}
	for(fit = sfd_get_domain_facetiterator(sfdu); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
		hig_facet *f = higfit_getfacet(fit);
		Point fcenter;
		hig_get_facet_center(f, fcenter);
		hig_cell *c = f->c;
		Point ccenter;
		hig_get_center(c, ccenter);
		int fgid = mp_lookup(mu, hig_get_fid(f));
		printf("=%d - %d - %d (%f, %f) (%f, %f)\n", myrank, hig_get_fid(f), fgid, ccenter[0], ccenter[1], fcenter[0], fcenter[1]);
	}
	DEBUG_PASS;
	higfit_destroy(fit);
	DEBUG_DIFF_TIME;
	for(int i = 0; i < 5; i++) {
	      dp_sync(dpp);
	      dp_sync(dpu);
	      dp_sync(dpv);
	//	if (myrank == 0) {
	//	}
	}
	dp_print("p", myrank, dpp);
	dp_print("u", myrank, dpu);
	dp_print("v", myrank, dpv);
	DEBUG_DIFF_TIME;
	dp_destroy(dpp);
	dp_destroy(dpu);
	DEBUG_DIFF_TIME;
	dp_destroy(dpv);
	DEBUG_DIFF_TIME;
}
