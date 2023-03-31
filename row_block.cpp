#include "mpi.h"
#include "myheader.h"

void scatterData(int rank, int nprcs, int n, float* A, float* b, float* x,
    float* As, float* bs, float* xs)
{
    int m = n / nprcs;
    for (int k = 0; k < m; k++) {
        MPI_Scatter(A + k * n * nprcs, n, MPI_FLOAT, As + k * n, n, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Scatter(b + k * nprcs, 1, MPI_FLOAT, bs + k, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Scatter(x + k * nprcs, 1, MPI_FLOAT, xs + k, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }
}

void gatherRes(int rank, int nprcs, int n, float* xs, float* x)
{
    int m = n / nprcs;
    for (int k = 0; k < m; k++) {
        MPI_Gather(xs + k, 1, MPI_FLOAT, x + k * nprcs, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }
}

void init_s(int rank, int nprcs, int n, int m, float* As, float* s)
{
    for (int i = 0; i < m; i++) {
        float sum = 0.0f;
        for (int j = rank + i * nprcs + 1; j < n; j++) {
            sum += As[i * n + j];
        }
        s[i] = sum;
    }
}
void debug(int i, int rank, int nprcs, int n, float* xs)
{
    float x[n];
    gatherRes(rank, nprcs, n, xs, x);
    if (rank == 0) {
        printX(i, n, x);
    }
}
void guass_seidel_para(int rank, int nprcs, int iter, int n, float* As, float* bs, float* xs)
{
    int m = n / nprcs;
    float* s = new float[m];

    init_s(rank, nprcs, n, m, As, s);

    for (int it = 0; it < iter; it++) {
        int k = 0;
        for (int i = 0; i < n; i++) {
            float xi;
            if (rank == i % nprcs) {
                xi = (bs[k] - s[k]) / As[k * n + i];

                xs[k] = xi;
                s[k] = 0.0f;

                k++;
            }
            MPI_Bcast(&xi, 1, MPI_FLOAT, i % nprcs, MPI_COMM_WORLD);

            for (int j = 0; j < m; j++) {
                if (rank + j * nprcs == i)
                    continue;
                s[j] = s[j] + As[j * n + i] * xi;
            }
        }
        // debug(it, rank, nprcs, n, xs);
    }

    delete[] s;
}

int main(int argc, char** argv)
{
    int rank, nprcs;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc != 3) {
        if (rank == 0)
            printf("Usage: exe itnum filepath!\n");
        return 0;
    }
    int n, m, itnum = atoi(argv[1]);

    // float A[16] = { 10, -1, 2, 0, -1, 11, -1, 3, 2, -1, 10, -1, 0, 3, -1, 8 };
    // float b[4] = { 6, 25, -11, 15 };
    // float x[4] = { 0, 0, 0, 0 };

    float *A, *b, *x_ref, *x_ans;
    float *As, *bs, *xs;

    if (rank == 0) {
        readFromFile(argv[2], n, A, b, x_ref, x_ans);
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    m = n / nprcs;

    As = new float[n * m];
    bs = new float[m];
    xs = new float[m];

    scatterData(rank, nprcs, n, A, b, x_ans, As, bs, xs);

    guass_seidel_para(rank, nprcs, itnum, n, As, bs, xs);

    gatherRes(rank, nprcs, n, xs, x_ans);

    if (rank == 0) {
        if (check(n, 1e-6, x_ref, x_ans)) {
            printf("The answer is right!\n");
            // for (int i = 0; i < n; i++)
            //     printf("%f ", x_ans[i]);
            // printf("\n");
        } else {
            printf("The answer is wrong !!!!!\n");
        }
        delete[] A, delete[] b, delete[] x_ref, delete[] x_ans;
    }
    delete[] As, delete[] bs, delete[] xs;

    MPI_Finalize();
    return 0;
}