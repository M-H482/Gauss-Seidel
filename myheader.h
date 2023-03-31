#include <cmath>
#include <iostream>
using namespace std;

void printX(int it, int n, float* x)
{
    printf("it = %d\n", it);
    for (int i = 0; i < n; i++)
        printf("%f  ", x[i]);
    printf("\n\n");
}

void readFromFile(char* file, int& n, float*& A, float*& b, float*& x_ref, float*& x_ans)
{
    FILE* fp = fopen(file, "r");

    int t, r;
    r = fscanf(fp, "n = %d", &t);
    n = t;

    A = new float[n * n];
    b = new float[n];
    x_ref = new float[n];
    x_ans = new float[n];

    for (int i = 0; i < n * n; i++) {
        r = fscanf(fp, "%f", &A[i]);
    }
    for (int i = 0; i < n; i++) {
        r = fscanf(fp, "%f", &b[i]);
    }
    for (int i = 0; i < n; i++) {
        r = fscanf(fp, "%f", &x_ref[i]);
        x_ans[i] = 0.0f;
    }
    fclose(fp);
}

void guass_seidel_serial(int iter, int n, float* A, float* b, float* x)
{
    float* s = new float[n];
    for (int i = 0; i < n; i++) {
        float t = 0.0f;
        for (int j = i + 1; j < n; j++) {
            t += A[i * n + j];
        }
        s[i] = t;
    }

    for (int it = 0; it < iter; it++) {
        for (int i = 0; i < n; i++) {
            x[i] = (b[i] - s[i]) / A[i * n + i];
            s[i] = 0.0f;
            for (int j = 0; j < n; j++) {
                if (j == i)
                    continue;
                s[j] = s[j] + A[j * n + i] * x[i];
            }
        }
        // printX(it, n, x);
    }
    delete[] s;
}

bool check(int n, float threshold, float* ref, float* ans)
{
    for (int i = 0; i < n; i++) {
        if (isnanf(ans[i]) || fabsf(ref[i] - ans[i]) > threshold)
            return false;
    }
    return true;
}