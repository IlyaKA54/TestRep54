//#include <iostream>
//#include <fstream>
//#include <cmath>
//#include <cstring>
//#include <omp.h>
//#include "LAB1_DLL.h"
//
//using namespace std;
//
//const int n = 2;
//const int m = 4;
//
//void PrintMethodResults(double* y, double deltaT, string methodName)
//{
//	cout << "-=" << methodName << "=-" << endl << " Время выполнения = " << deltaT << endl << " Результаты:" << endl;
//
//	for (int i = 0; i < n; i++)
//		cout << '\t' << "y[" << i << "] = " << y[i] << endl;
//
//	cout << endl;
//}
//
//double F1(double* y, double time) {
//	return y[1];
//}
//
//double F2(double* y, double time) {
//	return -0.1 * exp(-0.3 * time);
//}
//
//int main()
//{
//	setlocale(LC_ALL, "Russian");
//	double t0 = 0.0, t = t0, tMax = 10, tau = 0.01, y[n] = { 0.0, 1.0 }, tStart, tEnd, deltaT;
//
//	double (*funcs[])(double*, double) = { F1, F2 };
//
//	tStart = omp_get_wtime();
//	ExplicitEuler(y, n, funcs, t, tMax, tau);
//	tEnd = omp_get_wtime();
//	
//	deltaT = tEnd - tStart;
//
//	PrintMethodResults(y, deltaT, "Метод Эйлера (Явный)");
//
//	t = t0;
//	
//	y[0] = 0.0; y[1] = 1.0;
//	tStart = omp_get_wtime();
//
//	RungeKutta2(y, n, funcs, t, tMax, tau);
//
//	tEnd = omp_get_wtime();
//
//	deltaT = tEnd - tStart;
//
//	PrintMethodResults(y, deltaT, "Метод Рунге-Кутта 2");
//	
//	t = t0;
//	y[0] = 0.0; y[1] = 1.0;
//
//	tStart = omp_get_wtime();
//
//	PredictorCorrector(y, n, funcs, t, tMax, tau);
//
//	tEnd = omp_get_wtime();
//
//	deltaT = tEnd - tStart;
//
//	PrintMethodResults(y, deltaT, "Метод Предиктор-Корретор");
//
//	t = t0;
//	y[0] = 0.0; y[1] = 1.0;
//
//	tStart = omp_get_wtime();
//
//	RungeKutta4(y, n, funcs, t, tMax, tau);
//
//	tEnd = omp_get_wtime();
//
//	deltaT = tEnd - tStart;
//
//	PrintMethodResults(y, deltaT, "Метод Рунге-Кутта 4");
//
//	t = t0;
//	y[0] = 0.0; y[1] = 1.0;
//
//	tStart = omp_get_wtime();
//
//	ImplicitEuler(y, n, funcs, t, tMax, tau);
//
//	tEnd = omp_get_wtime();
//
//	deltaT = tEnd - tStart;
//
//	PrintMethodResults(y, deltaT, "Метод Эйлера (Неявный)"); 
//
//	system("Pause");
//
//	return 0;
//}
//
//
//#include <iostream>
//#include <fstream>
//#include <cmath>
//#include <cstring>
//#include <omp.h>
//#include <stdio.h>
//#include <complex>
//
//
//using namespace std;
//
//char filenameA[256];
//char filename[256];
//int matrixsizes[] = { 100, 250, 500 };
//const int N = 100;
//double tStart, tEnd, deltaT;
//
//void Jacobi(complex <double>** matrixA, complex <double>* matrixB, int n, complex <double>* matrixX) {
//	complex <double>* matrixXX{ new complex <double> [n] {} };
//	for (int i = 0; i < n; i++) matrixXX[i] = 1.0;
//	complex <double> sum = 0.0;
//	double xMax = -1.0, eps = 0.0001;
//	do {
//		xMax = -1.0;
//		for (int i = 0; i < n; i++)
//		{
//			sum = 0.0;
//			for (int j = 0; j < n; j++) {
//				if (j != i) sum += matrixA[i][j] * matrixX[j];
//			}
//			matrixXX[i] = (1.0 / matrixA[i][i]) * (matrixB[i] - sum);
//			complex <double> tmp = matrixX[i] - matrixXX[i];
//			if (xMax < abs(tmp)) xMax = abs(tmp);
//		}
//		for (int i = 0; i < n; i++) matrixX[i] = matrixXX[i];
//
//	} while (xMax > eps);
//	delete[]matrixXX;
//}
//
//void GaussSeidel(complex <double>** matrixA, complex <double>* matrixB, int n, complex <double>* matrixX) {
//	complex <double> sum = 0.0;
//	double eps = 0.0001, xMax = -1.0;
//	complex <double> count = 0.0;
//	do {
//		xMax = -1.0;
//		for (int i = 0; i < n; i++)
//		{
//			count = matrixX[i];
//			sum = 0.0;
//			for (int j = 0; j < n; j++) {
//				if (j != i) sum += matrixA[i][j] * matrixX[j];
//			}
//			matrixX[i] = (1.0 / matrixA[i][i]) * (matrixB[i] - sum);
//			complex <double> tmp = matrixX[i] - count;
//			if (xMax < abs(tmp)) xMax = abs(tmp);
//		}
//	} while (xMax > eps);
//}
//
//void Gauss(complex <double>** matrixA, complex <double>* matrixB, int n, complex <double>* matrixX) {
//	double eps = 0.0001, xMax = -1.0;
//	complex <double> w = 0.0, sum = 0.0;
//	for (int i = 0; i <= n - 1; i++)
//	{
//		for (int j = i + 1; j <= n - 1; j++)
//		{
//			w = matrixA[j][i] / matrixA[i][i];
//			for (int k = i; k <= n - 1; k++) matrixA[j][k] -= w * matrixA[i][k];
//			matrixB[j] -= w * matrixB[i];
//		}
//	}
//	for (int i = n - 1; i >= 0; i--)
//	{
//		sum = 0.0;
//		for (int j = i; j <= n - 1; j++) sum += matrixA[i][j + 1] * matrixX[j + 1];
//		matrixX[i] = (matrixB[i] - sum) / matrixA[i][i];
//	}
//}
//
//
//int main()
//{
//	ifstream matrA("C:/Users/ilya-/Desktop/ConsoleApplication3/matrixA(100).txt");
//	ifstream matrB("C:/Users/ilya-/Desktop/ConsoleApplication3/matrixB(100).txt");
//
//	if (!matrA.is_open() || !matrB.is_open()) {
//		cout << "error\n";
//		return 1;
//	}
//	complex <double>** matrixA = new complex <double>*[N];
//	for (int i = 0; i < N; i++) matrixA[i] = new complex <double>[N];
//	complex <double>* matrixB = new complex <double>[N];
//	complex <double>* matrixX = new complex <double>[N];
//	for (int i = 0; i < N; i++)
//	{
//		for (int j = 0; j < N; j++) matrA >> matrixA[i][j];
//		matrB >> matrixB[i];
//		matrixX[i] = 1.0;
//	}
//
//	GaussSeidel(matrixA, matrixB, N, matrixX);
//
//	char fileX[256];
//	sprintf_s(fileX, "%dmatrixX(%d).txt", 3, N);
//	ofstream matrX(fileX);
//	for (int i = 0; i < N; i++)
//	{
//		matrX << matrixX[i] << endl;
//	}
//
//	matrA.close();
//	matrB.close();
//	delete[]matrixA;
//	delete[]matrixB;
//	delete[]matrixX;





	/*srand(time(0));
	char fileA[256], fileB[256];
	int matrixSizes[] = { 100,250,500 };
	for (int N : matrixSizes)
	{
		sprintf_s(fileA, "matrixA(%d).txt", N);
		sprintf_s(fileB, "matrixB(%d).txt", N);
		ofstream matrixA(fileA);
		ofstream matrixB(fileB);

		if (!matrixA.is_open() || !matrixB.is_open()) {
			cout << "error" << endl;
			return 1;
		}

		complex <double>* matrA = new complex <double>[N];
		complex <double>* matrB = new complex <double>[N];

		for (int i = 0; i < N; i++) {
			complex <double> sum = 0.0;
			for (int j = 0; j < N; j++)
			{
				matrA[j].real(rand() % 100);
				matrA[j].imag(rand() % 100);
				sum += matrA[j];
			}
			matrA[i] = sum + 1.0;
			for (int k = 0; k < N; k++) matrixA << matrA[k] << "\t";
			matrixA << endl;

			matrB[i].real(rand() % 100);
			matrB[i].imag(rand() % 100);
			matrixB << matrB[i] << "\t";
		}

		matrixB.close();
		matrixA.close();
		delete[]matrB;
		delete[]matrA;
	}
	return 0;*/


//}


#include <iostream>
#include <omp.h>
#include <math.h>
#include <corecrt_math.h>

#define MAGIC_AUTOFIX_ENABLE

using namespace std;

double f(double x, double y)
{

    return 5.2 * x - 1.0 * y + exp(0.25 * x * x + 0.15 * y * y);

}

double Derivative(double (*func)(double, double), double x, double y, int i) {
    double dx = 0.0001;
    switch (i) {
    case 0:
        return (func(x + dx, y) - func(x, y)) / dx;
    case 1:
        return (func(x, y + dx) - func(x, y)) / dx;
    }

}

double* GradientDescent(double (*func)(double, double), double x, double y, double xmin, double xmax, double ymin, double ymax)
{
    double coordinates[2] = { x, y };
    double eps = 0.00001;
    double sum = 0;
    double delta = 0.00001;
    double grad[2] = { 0.0 };
    int iteration = 0;

    for (int i = 0; i < 2; i++)
        grad[i] = Derivative(func, coordinates[0], coordinates[1], i);


    do {
        

        sum = 0;
        for (int i = 0; i < 2; i++)
        {
            coordinates[i] = coordinates[i] - delta * grad[i];
        }

        if (coordinates[0] < xmin)
        {
            coordinates[0] = xmin;
            break;
        }

        if (coordinates[0] > xmax)
        {
            coordinates[0] = xmax;
            break;
        }

        if (coordinates[1] < ymin)
        {
            coordinates[1] = ymin;
            break;
        }

        if (coordinates[1] > ymax)
        {
            coordinates[1] = ymax;
            break;
        }

        for (int i = 0; i < 2; i++)
        {
            grad[i] = Derivative(func, coordinates[0], coordinates[1], i);
            sum += grad[i] * grad[i];
        }
        iteration++;

    } while (sum >= eps);

    return new double[2] {coordinates[0], coordinates[1]};



}


double* ScanMethod(double (*func)(double, double), double xmin, double xmax, double ymin, double ymax, double delta1)
{
    double xmax1 = xmax;
    double xmin1 = xmin;
    double ymax1 = ymax;
    double ymin1 = ymin;

    double x = xmin1, y = ymin1;
    double delta2 = delta1 / 100;
    double Xisc = x, Yisc = y, w, Xmin = 1000.0;

    for (int i = 0; i < (xmax1 - xmin1) / delta1 + 1; i++)
    {
        y = ymin1;
        for (int j = 0; j < (ymax1 - ymin1) / delta1 + 1; j++)
        {
            w = func(x, y);
            if (w < Xmin)
            {
                Xmin = w;
                Xisc = x;
                Yisc = y;
				
            }
            y += delta1;
        }
        x += delta1;
    }

    double xmin2, xmax2, ymin2, ymax2;

    xmin2 = Xisc - delta1;
    xmax2 = Xisc + delta1;
    ymin2 = Yisc - delta1;
    ymax2 = Yisc + delta1;

    if (xmin2 < xmin1)
        xmin2 = xmin1;
    if (xmax2 > xmax1)
        xmax2 = xmax1;
    if (ymin2 < ymin1)
        ymin2 = ymin1;
    if (ymax2 > ymax1)
        ymax2 = ymax1;

    Xmin = 10000.0;
    x = xmin2;
    y = ymin2;

    for (int i = 0; i < (xmax2 - xmin2) / delta2 + 1; i++)
    {
        y = ymin2;
        for (int j = 0; j < (ymax2 - ymin2) / delta2 + 1; j++)
        {
            w = func(x, y);
            if (w < Xmin)
            {
                Xmin = w;
                Xisc = x;
                Yisc = y;
            }
            y += delta2;
        }
        x += delta2;
    }

    return new double[3] { Xmin, Xisc, Yisc};
}

int main()
{
	MAGIC_AUTOFIX_ENABLE;
    auto resr = ScanMethod(f, -5.0, 5.0, -5.0, 5.0, 0.01);
    cout << "Min of func = " << resr[0] << endl << "Min x = " << resr[1] << endl << "Min y = " << resr[2] << endl;

    printf("\n");

    auto res = GradientDescent(f, 4, 4, -5.0, 5.0, -5.0, 5.0);
    cout << "Min of func = " << f(res[0], res[1]) << endl << "Min x = " << res[0] << endl << "Min y = " << res[1] << endl;
    system("Pause");
}


