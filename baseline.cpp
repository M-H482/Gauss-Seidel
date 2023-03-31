#include "myheader.h"

int main(int argc, char** argv)
{
    // float A[16] = { 10, -1, 2, 0, -1, 11, -1, 3, 2, -1, 10, -1, 0, 3, -1, 8 };
    // float b[4] = { 6, 25, -11, 15 };
    // float x[4] = { 0, 0, 0, 0 };

    if (argc != 3) {
        printf("Usage: exe itnum filepath!\n");
        return 0;
    }
    int n, itnum = atoi(argv[1]);

    float *A, *b, *x_ref, *x_ans;

    readFromFile(argv[2], n, A, b, x_ref, x_ans);

    guass_seidel_serial(itnum, n, (float*)A, b, x_ans);

    if (check(n, 1e-6, x_ref, x_ans)) {
        printf("The answer is right!\n");
        // for (int i = 0; i < n; i++)
        //     printf("%f ", x_ans[i]);
        // printf("\n");
    } else {
        printf("The answer is wrong !!!!!\n");
    }

    delete[] A, delete[] b, delete[] x_ref, delete[] x_ans;
    return 0;
}