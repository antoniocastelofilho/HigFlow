
#include<stdio.h>
#include<stdlib.h>

#include "mpi.h"
#include "higtree.h"
#include "higtree-io.h"
#include "higtree-parallel.h"

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	printf("rank = %d\n", rank);
	if (size != 2 && rank == 0) {
		printf("usage: mpirun -n 2 %s amr2d.file\n", argv[0]);
		exit(0);
	}
	if (rank == 0) {
		if (argc <= 1) {
			printf("usage: mpirun -n 2 %s amr2d.file\n", argv[0]);
			exit(0);
		}
		FILE * fdin = fopen(argv[1], "r");
		hig_cell * root = higio_read_from_amr2d(fdin);
		mp_mapper *m = mp_create();
		higcit_celliterator *it = higcit_create_all_leaves(root);
		mp_assign_from_celliterator(m, it, 1000);
		higcit_destroy(it);
		fclose(fdin);
		Point l, h;
		hig_get_lowpoint(root, l);
		hig_get_highpoint(root, h);
		MPI_Send(l, DIM, MPI_HIGREAL, 1, 0, MPI_COMM_WORLD);
		MPI_Send(h, DIM, MPI_HIGREAL, 1, 0, MPI_COMM_WORLD);
		higp_send_uniform_refinement_to_node(root, 1);
		h[0] += 0.0;
		hig_cell * root2 = hig_create_root(l, h);
		higp_receive_uniform_refinement_from_node(root2, 1);
		if (hig_equal(root, root2)) {
			printf("Equal!!\n");
		} else {
			printf("Not Equal!!\n");
		}
		higp_sync_mapper_to_node(root, m, 1);
		mp_destroy(m);
		hig_destroy(root);
		hig_destroy(root2);
	} else {
		Point l, h;
		MPI_Status status;
		MPI_Recv(l, DIM, MPI_HIGREAL, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(h, DIM, MPI_HIGREAL, 0, 0, MPI_COMM_WORLD, &status);
		hig_cell *root = hig_create_root(l, h);
		higp_receive_uniform_refinement_from_node(root, 0);
		higp_send_uniform_refinement_to_node(root, 0);
		mp_mapper *m = mp_create();
		higp_sync_mapper_from_node(root, m, 0);
		higcit_celliterator *it;
		for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
			hig_cell *c = higcit_getcell(it);
			uniqueid cid = hig_get_cid(c);
			int cgid = mp_lookup(m, cid);
			printf(" %d", cgid);
		}
		printf("\n");
		hig_destroy(root);
	}
	MPI_Finalize();
}
