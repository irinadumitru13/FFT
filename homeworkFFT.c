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
cplx *buf;	// result vector
cplx *aux;	// auxiliary vector
FILE* in;
FILE* out;

void getArgs(int argc, char **argv)
{
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

void init()
{
	if (fscanf(in, "%d", &N) == 0) {
		fprintf(stdout, "Couldn't read from file.\n");
		exit(1);
	}

	if ((xs = malloc(sizeof(double) * N)) == NULL) {
		fprintf(stdout, "malloc failed!");
		exit(1);
	}
	if ((buf = malloc(sizeof(cplx) * N)) == NULL) {
		fprintf(stdout, "malloc failed!");
		exit(1);
	}
	if ((aux = malloc(sizeof(cplx) * N)) == NULL) {
		fprintf(stdout, "malloc failed!");
		exit(1);
	}

	for (int i = 0; i < N; i++) {
		if (fscanf(in, "%lf", &xs[i]) == 0) {
				fprintf(stdout, "Couldn't read from file.\n");
				exit(1);
		}

		buf[i] = xs[i] + 0.0 * I;
		aux[i] = buf[i];
	}

	PI = atan2(1, 1) * 4;
}

// FFT function
void fft(cplx* buf, cplx* aux, int step)
{
	if (step < N) {
		fft(aux, buf, step * 2);
		fft(aux + step, buf + step, step * 2);

		for (int i = 0; i < N; i += 2 * step) {
			cplx t = cexp(-I * PI * i / N) * aux [i + step];
			buf[i / 2] = aux[i] + t;
			buf[(i + N) / 2] = aux[i] - t;
		}
	}
}

// the number of threads is seen as the number of subtrees that
// can be executed in parallel
void* threadFunction(void *args)
{
	int thread_id = *(int*) args;

	// case for 1 and 4 threads
	if (((int)log2(P) % 2) == 0) {
		fft(buf + thread_id, aux + thread_id, P);
	// case for 2 threads
	} else {
		fft(aux + thread_id, buf + thread_id, P);
	}

	return NULL;
}

// last for-loop from FFT function
void unify(cplx* buf, cplx* aux, int step) {
	for (int i = 0; i < N; i += 2 * step) {
		cplx t = cexp(-I * PI * i / N) * aux [i + step];
		buf[i / 2] = aux[i] + t;
		buf[(i + N) / 2] = aux[i] - t;
	}
}

void getResult() {
	if (P == 4) {
		// the results from level 3 must be unified
		unify(aux, buf, P / 2);
		unify(aux + 1, buf + 1, P / 2);

		// unify results from level 2
		unify(buf, aux, P / 4);
	} else if (P == 2) {
		// unify results from level 2
		unify(buf, aux, P / 2);
	}
}

void displayResults()
{
	fprintf(out, "%d\n", N);

	for (int i = 0; i < N; i++) {
		fprintf(out, "%lf %lf\n", creal(buf[i]), cimag(buf[i]));
	}

	free(xs);
	free(buf);
	free(aux);
}

int main(int argc, char * argv[])
{
	int i;
	getArgs(argc, argv);
	init();

	pthread_t tid[P];
	int thread_id[P];

	for (i = 0; i < P; i++) {
		thread_id[i] = i;
	}

	for (i = 0; i < P; i++) {
		pthread_create(&(tid[i]), NULL, threadFunction, &(thread_id[i]));
	}

	for (i = 0; i < P; i++) {
		pthread_join(tid[i], NULL);
	}

	getResult();
	displayResults();

	return 0;
}
