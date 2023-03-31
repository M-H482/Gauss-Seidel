#include "myheader.h"

void matVecMul(int n, float* A, float* x, float* b)
{
    for (int i = 0; i < n; i++) {
        float s = 0.0f;
        for (int j = 0; j < n; j++) {
            s += A[i * n + j] * x[j];
        }
        b[i] = s;
    }
}

bool generate(int n, float* A, float* b, float* x_ref, float* x_ans)
{
    srand((unsigned int)time(NULL));

    for (int i = 0; i < n; i++) {
        int sum = 0;
        for (int j = 0; j < n; j++) {
            if (i == j)
                continue;
            A[i * n + j] = rand() % 15;
            sum += A[i * n + j];
        }
        A[i * n + i] = sum + 1;
    }

    for (int i = 0; i < n; i++)
        x_ref[i] = rand() % 5;

    matVecMul(n, A, x_ref, b);

    guass_seidel_serial(50, n, A, b, x_ans);

    for (int i = 0; i < n; i++) {
        if (isnanf(x_ans[i]) || fabsf(x_ref[i] - x_ans[i]) > 1e-6)
            return false;
    }
    return true;
}

void writeToFile(int n, float* A, float* b, float* x)
{
    char file[128];
    sprintf(file, "./data/data_%d.txt", n);

    printf("write to file:%s\n", file);

    FILE* fp = fopen(file, "w");
    fprintf(fp, "n = %d\n", n);

    for (int i = 0; i < n * n; i++) {
        fprintf(fp, "%f\n", A[i]);
    }
    for (int i = 0; i < n; i++) {
        fprintf(fp, "%f\n", b[i]);
    }
    for (int i = 0; i < n; i++) {
        fprintf(fp, "%f\n", x[i]);
    }
    fclose(fp);
}
int main(int argc, char** argv)
{
    if (argc != 2) {
        printf("Usage: exe N\n");
        return 0;
    }

    int n = atoi(argv[1]);
    float* A = new float[n * n];
    float* b = new float[n];
    float* x_ref = new float[n];
    float* x_ans = new float[n];

    int cnt = 0, maxCnt = 10000;
    while (!generate(n, A, b, x_ref, x_ans) && cnt < maxCnt)
        cnt++;

    if (cnt == maxCnt) {
        printf("generation failed!\n");
    } else {
        printf("Generation succeeded!\n");
        writeToFile(n, A, b, x_ans);
    }
    return 0;
}