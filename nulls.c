// gcc nulls.c -o nm -lgsl -lgslcblas -O3 -DHAVE_INLINE
// or use clang

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define FNSIZE 128

int checkadj(int *adj, int nrow, int ncol) {
	int *msRow = (int*) malloc(nrow * sizeof(int));
	int *msCol = (int*) malloc(ncol * sizeof(int));
	// Check the presence of interactions
	unsigned int noint = 0;
	int nr, nc;
	for (nr = 0; nr < nrow; ++nr) {
		for (nc = 0; nc < ncol; ++nc) {
			if (adj[nc + nr * ncol] == 1) {
				msRow[nr] = 1;
				msCol[nc] = 1;
			}
		}
	}
	// Count interactions for each trophic level
	for (nr = 0; nr < nrow; ++nr) {
		if (msRow[nr] == 0) {
			++noint;
			break;
		}
	}
	for (nc = 0; nc < ncol; ++nc) {
		if (msCol[nc] == 0) {
			++noint;
			break;
		}
	}
	return noint;
}

void null_1(int *template, int nrow, int ncol, int *nullweb, gsl_rng *rng) {
	int nc, nr;
	// Check the total number of interactions
	unsigned int L = 0;
	int i;
	for (i = 0; i < (nrow * ncol); ++i) {
		// Loop through all interactions
		if (template[i] == 1) {
			++L;
		}
	}
	// Assign global interaction probability
	double Pint;
	Pint = (double) L / (nrow * ncol);
	// Loop through all the interactions
	for (nr = 0; nr < nrow; ++nr) {
		for (nc = 0; nc < ncol; ++nc) {
			if (gsl_rng_uniform(rng) <= Pint) {
				nullweb[nc + nr * ncol] = 1;
			}
		}
	}
}

void null_2(int *template, int nrow, int ncol, int *nullweb, gsl_rng *rng) {
	double *msRow = (double*) malloc(nrow * sizeof(double));
	double *msCol = (double*) malloc(ncol * sizeof(double));
	int nr, nc;
	for (nr = 0; nr < nrow; ++nr) {
		msRow[nr] = 0.0;
	}
	for (nc = 0; nc < ncol; ++nc) {
		msCol[nc] = 0.0;
	}
	for (nr = 0; nr < nrow; ++nr) {
		for (nc = 0; nc < ncol; ++nc) {
			if (template[nc + nr * ncol] == 1) {
				msRow[nr] += 1.0;
				msCol[nc] += 1.0;
			}
		}
	}
	// Calculate the interaction probability per row and column
	for (nr = 0; nr < nrow; ++nr) {
		msRow[nr] = msRow[nr] / (double) ncol;
	}
	for (nc = 0; nc < ncol; ++nc) {
		msCol[nc] = msCol[nc] / (double) nrow;
	}
	// Loop through all the interactions
	for (nr = 0; nr < nrow; ++nr) {
		for (nc = 0; nc < ncol; ++nc) {
			if (gsl_rng_uniform(rng) <= ((msRow[nr] + msCol[nc]) / (double) 2)) {
				nullweb[nc + nr * ncol] = 1;
			}
		}
	}
}

int main(int argc, char *argv[]) {
	clock_t start, stop;
	// Initiate the GSL random number generator
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);
	gsl_rng_set(rng, time(NULL));

	/* Start the main program
	 * Arguments order:
	 * template_file ncol nrow target maxiter - e.g.
	 * nm web.txt 10 10 1000 2000
	 * will run 1000 null on a web contained in web.txt,
	 * with 10 rows and 10 columns, and will stop after
	 * 2000 iterations no matter what
	 */

	// Command line arguments
	unsigned int ncol, nrow, target, maxiter;
	ncol = atoi(argv[2]);
	nrow = atoi(argv[3]);
	target = atoi(argv[4]);
	maxiter = atoi(argv[5]);

	// Number of success and trials
	int iter, dn1, dn2;
	iter = 0;
	dn1 = 0;
	dn2 = 0;

	// Read the template file and convert to adjacency matrix if needed
	int *web = (int*) malloc(ncol * nrow * sizeof(int));
	FILE *template;
	template = fopen(argv[1], "rt");
	int nc, nr;
	for (nr = 0; nr < nrow; ++nr) {
		for (nc = 0; nc < ncol; ++nc) {
			fscanf(template, "%d", &web[nc + nr * ncol]);
			web[nc + nr * ncol] = (web[nc + nr * ncol] > 0) ? 1 : 0;
		}
	}
	fclose(template);

	// Pointers to output files
	FILE *n1out;
	FILE *n2out;

	// Start the iterations
	start = clock();
	int i;
	while ((iter < maxiter) & ((dn1 < target) | (dn2 < target))) {
		++iter;
		if (dn1 < target) {
			int* t_n1 = (int*) malloc(nrow * ncol * sizeof(int));
			for (i = 0; i < (nrow * ncol); ++i) {
				t_n1[i] = 0;
			}
			null_1(web, nrow, ncol, t_n1, rng);
			if (checkadj(t_n1, nrow, ncol) == 0) {
				char tfname[FNSIZE]; // The filename buffer.
				snprintf(tfname, sizeof(char) * FNSIZE, "%s_n1_%i.txt", argv[1], dn1);
				n1out = fopen(tfname, "w");
				for (nr = 0; nr < nrow; ++nr) {
					for (nc = 0; nc < ncol; ++nc) {
						fprintf(n1out, "%d ", t_n1[nc + nr * ncol]);
					}
					fprintf(n1out, "\n");
				}
				fclose(n1out);
				++dn1;
			}
			free(t_n1);
		}

		if (dn2 < target) {
			int* t_n2 = (int*) malloc(nrow * ncol * sizeof(int));
			for (i = 0; i < (nrow * ncol); ++i) {
				t_n2[i] = 0;
			}
			null_2(web, nrow, ncol, t_n2, rng);
			if (checkadj(t_n2, nrow, ncol) == 0) {
				char tfname[FNSIZE]; // The filename buffer.
				snprintf(tfname, sizeof(char) * FNSIZE, "%s_n2_%i.txt", argv[1], dn1);
				n2out = fopen(tfname, "w");
				for (nr = 0; nr < nrow; ++nr) {
					for (nc = 0; nc < ncol; ++nc) {
						fprintf(n2out, "%d ", t_n2[nc + nr * ncol]);
					}
					fprintf(n2out, "\n");
				}
				fclose(n2out);
				++dn2;
			}
			free(t_n2);
		}
	}
	stop = clock();

	printf("%d type I and %d type II webs generated in %f s (%d iterations).\n", dn1, dn2, (stop - start) / (float) CLOCKS_PER_SEC, iter);

	gsl_rng_free(rng);
	free(web);
	return EXIT_SUCCESS;
}
