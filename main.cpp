//Michał Czarnobaj 188816
#include <iostream>
#include <chrono>
#include <math.h>
#include <fstream>

int N = 916;
double a1 = 13, a2 = -1, a3 = -1;
const double maxNorm = pow(10, -9);
using namespace std;

void _free(double** A);
double** initMatrix();
double* bVector(int);
double normalize(double*);
void gauss_seidel(double**, double*, chrono::duration<double>*, int*, double*, bool, bool, int);
void jacobi(double**, double*, chrono::duration<double>*, int*, double*, bool, bool, int);
void LU_Decomposition(double**, double*, chrono::duration<double>*, double*);
double* mulMatVec(double** M, double* v);
double* subVec(double* v1, double* v2);

struct Graph {
    int n;
    chrono::duration<double> time;
};

struct GraphIter {
    int n;
    int iterations;
};

struct GraphNorm {
    int n;
    double resNorm;
};

int main()
{
    double** A = initMatrix();
    double* b = bVector(9);
    chrono::duration<double> time;
    int iterations = 0;
    double normRes = 0;
    printf(" ----------------------------------\n");
    printf(" Zadanie B: N = 916, a1 = 13, f = 9\n");
    printf(" ----------------------------------\n");
    jacobi(A, b, &time, &iterations, &normRes, false, true, 0); // by utworzyc plik z danymi do wykresu, przedostatni argument true
    gauss_seidel(A, b, &time, &iterations, &normRes, false, true, 0); // ostatni arguement = 0 dla zadania B, argument = 1 dla zadania C, w innym przypadku dowolna liczba

    a1 = 3;
    _free(A);
    A = initMatrix();

    printf(" ----------------------------------\n");
    printf(" Zadanie C: N = 916, a1 = 3, f = 9\n");
    printf(" ----------------------------------\n");
    jacobi(A, b, &time, &iterations, &normRes, false, true, 1);
    gauss_seidel(A, b, &time, &iterations, &normRes, false, true, 1);

    /*printf(" ----------------------------------\n");
    printf(" Zadanie D: N = 916, a1 = 3, f = 9\n");
    printf(" ----------------------------------\n");

    LU_Decomposition(A, b, &time, &normRes);*/

    /*printf(" ----------------------------------\n");
    printf(" Zadanie E: a1 = 13, f = 9\n");
    printf(" ----------------------------------\n");

    a1 = 13;

    int tabOfN[] = {100, 500, 1000, 2000, 3000, 4000, 5000};
    Graph setOfGraphs[3][7];
    GraphIter setOfGraphsIter[2][7];
    GraphNorm setOfGraphsNorm[3][7];
    _free(A);

    for (int i = 0; i < 7; i++) {
        N = tabOfN[i];
        printf(" N = %d\n\n", N);
        A = initMatrix();
        jacobi(A, b, &time, &iterations, &normRes, false, false, 2);
        setOfGraphs[0][i].n = N;
        setOfGraphs[0][i].time = time;
        setOfGraphsIter[0][i].n = N;
        setOfGraphsIter[0][i].iterations = iterations;
        setOfGraphsNorm[0][i].n = N;
        setOfGraphsNorm[0][i].resNorm = normRes;
        gauss_seidel(A, b, &time, &iterations, &normRes, false, false, 2);
        setOfGraphs[1][i].n = N;
        setOfGraphs[1][i].time = time;
        setOfGraphsIter[1][i].n = N;
        setOfGraphsIter[1][i].iterations = iterations;
        setOfGraphsNorm[1][i].n = N;
        setOfGraphsNorm[1][i].resNorm = normRes;
        LU_Decomposition(A, b, &time, &normRes, false);
        setOfGraphsNorm[2][i].n = N;
        setOfGraphsNorm[2][i].resNorm = normRes;
        setOfGraphs[2][i].n = N;
        setOfGraphs[2][i].time = time;
        _free(A);
        if (i == 6) {
            ofstream jacobiTime, gaussTime, luTime, jacobiIter, gaussIter, jacobiNorm, gaussNorm, luNorm;
            jacobiTime.open("jacobiTime.csv");
            gaussTime.open("gaussTime.csv");
            luTime.open("luTime.csv");
            jacobiIter.open("jacobiIter.csv");
            gaussIter.open("gaussiIter.csv");
            jacobiNorm.open("jacobiNorm.csv");
            gaussNorm.open("gaussNorm.csv");
            luNorm.open("luNorm.csv");
            for (int i = 0; i < 7; i++) {
                jacobiTime << setOfGraphs[0][i].n << ";";
                gaussTime << setOfGraphs[1][i].n << ";";
                luTime << setOfGraphs[2][i].n << ";";
                jacobiIter << setOfGraphsIter[0][i].n << ";";
                gaussIter << setOfGraphsIter[1][i].n << ";";
                jacobiNorm << setOfGraphsNorm[0][i].n << ";";
                gaussNorm << setOfGraphsNorm[1][i].n << ";";
                luNorm << setOfGraphsNorm[2][i].n << ";";
            }
            jacobiTime << endl;
            gaussTime << endl;
            luTime << endl;
            jacobiIter << endl;
            gaussIter << endl;
            jacobiNorm << endl;
            gaussNorm << endl;
            luNorm << endl;
            for (int i = 0; i < 7; i++) {
                jacobiTime << setOfGraphs[0][i].time.count() << ";";
                gaussTime << setOfGraphs[1][i].time.count() << ";";
                luTime << setOfGraphs[2][i].time.count() << ";";
                jacobiIter << setOfGraphsIter[0][i].iterations << ";";
                gaussIter << setOfGraphsIter[1][i].iterations << ";";
                jacobiNorm << setOfGraphsNorm[0][i].resNorm << ";";
                gaussNorm << setOfGraphsNorm[1][i].resNorm << ";";
                luNorm << setOfGraphsNorm[2][i].resNorm << ";";
            }
            jacobiTime.close();
            gaussTime.close();
            luTime.close();
            jacobiIter.close();
            gaussIter.close();
            jacobiNorm.close();
            gaussNorm.close();
            luNorm.close();
        }
    }*/
    return 0;
}



double** initMatrix() {
    double** matrix = new double* [N];
    for (int i = 0; i < N; i++) {
        matrix[i] = new double[N];
        for (int j = 0; j < N; j++) {
            if (j == i) matrix[i][j] = a1;
            else if (j == i - 1 || j == i + 1) matrix[i][j] = a2;
            else if (j == i - 2 || j == i + 2) matrix[i][j] = a3;
            else matrix[i][j] = 0;
        }
    }
    return matrix;

}

double* bVector(int f) {
    double* b = new double[N];
    for (int i = 0; i < N; i++) {
        b[i] = sin(i * f);
    }
    return b;
}

double normalize(double* vector) {
    double result = 0;
    for (int i = 0; i < N; i++) {
        result += pow(vector[i], 2);
    }
    return sqrt(result);
}

void gauss_seidel(double** A, double* b, chrono::duration<double>* time, int* iterations, double* normRes, bool statistics, bool graph, int task) {
    int iteration = 0;
    static int call = 0;
    double* n{ new double[N] {} };
    double* r = new double[N];
    double* v;
    bool flag = false;
    for (int i = 0; i < N; i++) {
        r[i] = 1;
    }
    GraphNorm normGraph[16];
    GraphNorm nGraph[509];
    v = mulMatVec(A, r);
    double* res = subVec(v, b);
    double norm = normalize(res);
    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();
    while (norm > maxNorm) {
        for (int i = 0; i < N; i++) {
            n[i] = b[i] / A[i][i];
            for (int j = 0; j < N; j++) {
                if (j == i) continue;
                n[i] -= ((A[i][j] / A[i][i]) * r[j]);
                r[i] = n[i];
            }
        }
        delete[] res;
        delete[] v;
        v = mulMatVec(A, r);
        res = subVec(v, b);
        norm = normalize(res);
        if (statistics) cout << norm << endl;
        iteration++;
        if (graph) {
            if (task == 0) {
                normGraph[iteration - 1].n = iteration;
                normGraph[iteration - 1].resNorm = norm;
            }
            else if (task == 1) {
                nGraph[iteration - 1].n = iteration;
                nGraph[iteration - 1].resNorm = norm;
            }
        }
        if (norm >= INFINITY) {
            flag = true;
            break;
        }
    }

    end = chrono::system_clock::now();
    chrono::duration<double>duration = end - start;
    *time = duration;
    printf(" \n Gauss-Seidel\n\n");
    cout << " Time elapsed: " << time->count() << "s\n";
    printf(" Iterations: %d\n", iteration);
    if (flag) printf(" Could not find a solution, the function is divirgent\n");
    cout << " Residual norm: " << norm << "\n\n";
    *iterations = iteration;
    *normRes = norm;
    if (graph) {
        ofstream file;
        if (task == 0) {
            file.open("gaussNormB.csv");
            for (int i = 0; i < 16; i++)
                file << normGraph[i].n << ";";
            file << endl;
            for (int i = 0; i < 16; i++)
                file << normGraph[i].resNorm << ";";
            file.close();
        }
        else if (task == 1) {
            file.open("gaussNormC.csv");
            for (int i = 0; i < 509; i++)
                file << nGraph[i].n << ";";
            file << endl;
            for (int i = 0; i < 509; i++)
                file << nGraph[i].resNorm << ";";
            file.close();
        }

    }
    delete[] res;
    delete[] v;
    delete[] r;
    delete[] n;
}

void jacobi(double** A, double* b, chrono::duration<double>* time, int* iterations, double* normRes, bool statistics, bool graph, int task) {
    int iteration = 0;
    static int call = 0;
    double* n = new double[N];
    double* r = new double[N];
    double* r1 = new double[N];
    double* v;
    bool flag = false;
    for (int i = 0; i < N; i++) {
        r[i] = 1;
        r1[i] = 1;
        n[i] = 0;
    }
    GraphNorm normGraph[23];
    GraphNorm nGraph[1223];
    v = mulMatVec(A, r);
    double* res = subVec(v, b);
    double norm = normalize(res);
    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();
    while (norm > maxNorm) {
        for (int i = 0; i < N; i++) {
            r[i] = b[i];
            for (int j = 0; j < N; j++) {
                if (j != i) r[i] -= A[i][j] * r1[j];
            }
            r[i] /= A[i][i];
        }
        for (int i = 0; i < N; i++) r1[i] = r[i];
        delete[] res;
        delete[] v;
        v = mulMatVec(A, r);
        res = subVec(v, b);
        norm = normalize(res);
        if (statistics) cout << norm << endl;
        iteration++;
        if (graph) {
            if (task == 0) {
                normGraph[iteration - 1].n = iteration;
                normGraph[iteration - 1].resNorm = norm;
            }
            else if (task == 1) {
                nGraph[iteration - 1].n = iteration;
                nGraph[iteration - 1].resNorm = norm;
            }
        }
        if (norm >= INFINITY) {
            flag = true;
            break;
        }
    }

    end = chrono::system_clock::now();
    chrono::duration<double>duration = end - start;
    *time = duration;
    printf(" \n Jacobi\n\n");
    cout << " Time elapsed: " << time->count() << "s\n";
    printf(" Iterations: %d\n", iteration);
    if (flag) printf(" Could not find a solution, the function is divirgent\n");
    cout << " Residual norm: " << norm << "\n\n";
    *iterations = iteration;
    *normRes = norm;
    if (graph) {
        ofstream file;
        if (task == 0) {
            file.open("jacobiNormB.csv");
            for (int i = 0; i < 23; i++)
                file << normGraph[i].n << ";";
            file << endl;
            for (int i = 0; i < 23; i++)
                file << normGraph[i].resNorm << ";";
            file.close();
        }
        else if (task == 1) {
            file.open("jacobiNormC.csv");
            for (int i = 0; i < 1223; i++)
                file << nGraph[i].n << ";";
            file << endl;
            for (int i = 0; i < 1223; i++)
                file << nGraph[i].resNorm << ";";
            file.close();
        }
    }
    delete[] res;
    delete[] v;
    delete[] r;
    delete[] r1;
    delete[] n;
}

void LU_Decomposition(double** A, double* b, chrono::duration<double>* time, double* normRes) {
    double** L = new double* [N];
    double** U = new double* [N];
    double* z = new double[N];
    double* x = new double[N];
    for (int i = 0; i < N; i++) {
        L[i] = new double[N];
        U[i] = new double[N];
        z[i] = 0;
        x[i] = 0;
        for (int j = 0; j < N; j++) {
            L[i][j] = 0;
            U[i][j] = 0;
        }
    }
    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();
    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            double sum = 0;
            for (int k = 0; k < i; k++) {
                sum += L[i][k] * U[k][j];
            }
            U[i][j] = A[i][j] - sum;
        }
        for (int j = i; j < N; j++) {
            if (i == j) L[i][i] = 1;
            else {
                double sum = 0;
                for (int k = 0; k < i; k++) {
                    sum += L[j][k] * U[k][i];
                }
                L[j][i] = (A[j][i] - sum) / U[i][i];
            }
        }
    }
    for (int i = 0; i < N; i++) {
        double sum = 0;
        for (int j = 0; j < i; j++) {
            sum += L[i][j] * z[j];
        }
        z[i] = (b[i] - sum) / L[i][i];
    }
    for (int i = N - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = N - 1; j > i; j--) {
            sum += U[i][j] * x[j];
        }
        x[i] = (z[i] - sum) / U[i][i];
    }
    end = chrono::system_clock::now();
    chrono::duration<double>duration = end - start;
    *time = duration;
    double* res = subVec(mulMatVec(A, x), b);
    double norm = normalize(res);
    printf(" \n LU Decomposition\n\n");
    cout << " Time elapsed: " << time->count() << "s\n";
    cout << " Residual norm: " << norm << "\n\n";
    *normRes = norm;
    delete[] x;
    delete[] z;
    delete[] res;
    for (int i = 0; i < N; i++) {
        delete[] L[i];
        delete[] U[i];
    }
    delete[] U;
    delete[] L;
}

double* mulMatVec(double** M, double* v) {
    double* vector = new double[N];
    for (int i = 0; i < N; i++) {
        vector[i] = 0;
        for (int j = 0; j < N; j++) {
            vector[i] += M[i][j] * v[j];
        }
    }
    return vector;
}

double* subVec(double* v1, double* v2) {
    double* vector = new double[N];
    for (int i = 0; i < N; i++) {
        vector[i] = v1[i] - v2[i];
    }
    return vector;
}

void _free(double** A) {
    for (int i = 0; i < N; i++) {
        delete[] A[i];
    }
    delete[] A;
}
