#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <complex.h>

typedef double complex cplx;

int N;
int P;
double PI;
double *xs; // input data vector
cplx *X; // result vector
FILE* in;
FILE* out;

void getArgs(int argc, char **argv) {
	if (argc < 3) {
		fprintf(stdout, "Usage: %s <in.txt> <out.txt> <numThreads>", argv[0]);
		exit(1);
	}

	if ((in = fopen(argv[1], "rt")) == NULL) {
		fprintf(stdout, "Failed to open input file.\n");
		exit(1);
	}
	if ((out = fopen(argv[2], "wt")) == NULL) {
		fprintf(stdout, "Failed to open the output file.\n");
		exit(1);
	}

	P = atoi(argv[3]);
}

void init() {
	int i;

	if (fscanf(in, "%d", &N) == 0) {
		fprintf(stdout, "Couldn't read from file.\n");
		exit(1);
	}

	if ((xs = malloc(sizeof(double) * N)) == NULL) {
		fprintf(stdout, "malloc failed!");
		exit(1);
	}
	if ((X = malloc(sizeof(cplx) * N)) == NULL) {
		fprintf(stdout, "malloc failed");
		exit(1);
	}

	for (i = 0; i < N; i++) {
		if (fscanf(in, "%lf", &xs[i]) == 0) {
				fprintf(stdout, "Couldn't read from file.\n");
				exit(1);
		}

		X[i] = 0.0 + 0.0 * I;
	}

	PI = atan2(1, 1) * 4;
}

void displayResults() {
	fprintf(out, "%d\n", N);

	for (int i = 0; i < N; i++) {
			fprintf(out, "%lf %lf\n", creal(X[i]), cimag(X[i]));
	}

	free(xs);
	free(X);
}

// using the formula
// Xk = sum(xn * e ^ (-2 * PI * i / N * k * n)), n = 0 .. N - 1
// Parallelize the outer for-loop
void* DFT(void *args) {
	int k, n;
	int thread_id = *(int*) args;

	int start = ceil(1.0f * N / P) * thread_id;
	int end = fmin(N, ceil(1.0f * N / P) * (thread_id + 1));

	for (k = start; k < end; k++) {
		for (n = 0; n < N; n++) {
			X[k] += cexp(-I * 2 * PI * k * n / N) * xs[n];
		}
	}

	return NULL;
}

int main(int argc, char * argv[]) {
	int i;
	getArgs(argc, argv);
	init();

	pthread_t tid[P];
	int thread_id[P];

	for (i = 0; i < P; i++) {
		thread_id[i] = i;
	}

	for (i = 0; i < P; i++) {
		pthread_create(&(tid[i]), NULL, DFT, &(thread_id[i]));
	}

	for (i = 0; i < P; i++) {
		pthread_join(tid[i], NULL);
	}

	displayResults();

	return 0;
}
