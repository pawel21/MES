#include <cmath>
#include <iostream>
#include "RozwRown.h"
using namespace std;

//RozwRown::RozwRown(void)
//{
//}
//RozwRown::~RozwRown(void)
//{
//}

void RozwRown::cholesky(double **macierz, double *wynik, const int ileWierszy, const int ileKolumn)
{
	/*Rozwi�zaywanie r�wna� metod� Choleskiego
	posta� r�wnania: Ax=b
	nale�y poda� jako argumenty funkcji kolejno
	- macierz wsp��czynnik�w (A)
	- wektor niewiadomych (x)
	- wektor wyraz�w wolnych (b)
	- ilo�� wierszy macierzy A
	- ilo�� kolumn macierzy A

	macierz A powinna by� postaci
	{ a[0][0]  a[0][1]  ...  a[0][n]  ;
	a[1][1]  a[1][2]  ...  a[1][n+1];
	.
	.
	.
	a[m-1][m-1]  a[m-1][m]   0  ...  0;
	a[m][m]      0           0  ...  0}

	gdzie n jest ilo�ci� pasm macierzy natomiast m jest ilo�ci� wierzy macierzy A
	*/


	//rozk�ad LL^T
	//tworzeni macierzy L
	double **L;
	L = new double *[ileWierszy];
	for (int i = 0; i<ileWierszy; i++)
		L[i] = new double[ileKolumn];

	//inicjalizacja zerowych wyraz�w
	for (int i = 0; i < ileWierszy; i++)
		for (int j = 0; j < ileKolumn; j++)
		{
			L[i][j] = 0.0;
		}

	double sum = 0;
	for (int j = 0; j < ileWierszy; j++)
	{
		for (int i = 0; i < ileWierszy; i++)
		{
			if (i == j)
			{
				sum = 0;
				for (int k = 0; k < i; k++)
				{
					if (i - k < ileKolumn){
						sum += L[i][i - k] * L[i][i - k];
					}
				}

				L[i][0] = sqrt((macierz[ileWierszy - i - 1][0] - sum));

			}
			else if (j<i && abs(i - j) <= ileKolumn)
			{
				sum = 0;
				for (int k = 0; k < j; k++)
				{
					if (i - k < ileKolumn && j - k < ileKolumn){
						sum += L[j][j - k] * L[i][i - k];
					}
				}

				if (i - j < ileKolumn){
					L[i][i - j] = (macierz[ileWierszy - i - 1][i - j] - sum) / L[j][0];
				}
			}
		}
	}

	//rozwiazanie L
	double *y;
	y = new double[ileWierszy];

	y[0] = macierz[0][2] / L[0][0];
	for (int i = 1; i < ileWierszy; i++)
	{
		double sum = 0;
		for (int k = 0; k < i; k++)
		{
			if (i - k < ileKolumn){
				sum += L[i][i - k] * y[k];
			}
		}
		y[i] = (macierz[i][2] - sum) / L[i][0];
	}


	//rozwi�zanie L^T
	double *Y;
	Y = new double[ileWierszy];
	Y[ileWierszy - 1] = y[ileWierszy - 1] / L[ileWierszy - 1][0];
	for (int i = ileWierszy - 1; i >= 0; i--)
	{
		double sum = 0;
		for (int k = ileWierszy - 1; k > i; k--)
		{
			if (k - i < ileKolumn){
				sum += (L[k][k - i] * wynik[k]);
			}
		}
		wynik[i] = (y[i] - sum) / L[i][0];
	}

	//usuwanie zb�dnych zmiennych
	for (int i = 0; i<ileWierszy; i++)
		delete[] L[i];

	delete[] L;
	delete[] Y;
	delete[] y;
}

void RozwRown::przemiatania(double **macierz, double *wynik, const int ileWierszy, const int ileKolumn)
{
	//tworzenie niezb�dnych zmiennych
	double *beta;
	beta = new double[ileWierszy];

	double *gamma;
	gamma = new double[ileWierszy];

	//wyliczenie wsp��czynnikow gamma i beta
	beta[0] = macierz[0][0];
	gamma[0] = macierz[0][2]/beta[0];

	for (int i = 1; i < ileWierszy; i++)
	{
		beta[i] = macierz[i][0] - macierz[i - 1][1] * macierz[i - 1][1] / beta[i - 1];
		gamma[i] = (macierz[i][2] - macierz[i - 1][1] * gamma[i - 1]) / beta[i];
	}

	//znalezienie rozwi�zania
	wynik[ileWierszy - 1] = gamma[ileWierszy - 1];

	for (int i = ileWierszy - 2; i >= 0; i--)
	{
		wynik[i] = gamma[i] - macierz[i][1] * wynik[i + 1] / beta[i];
	}

	delete[] beta;
	delete[] gamma;
}

void RozwRown::iteracyjna(double **macierz, double *wynik, const int ileWierszy, const int ileKolumn )
{
	double **C;
	C = new double *[ileWierszy];
	for (int i = 0; i<ileWierszy; i++)
		C[i] = new double[ileKolumn];

	//inicjalizacja zerowych wyraz�w
	for (int i = 0; i < ileWierszy; i++)
		for (int j = 0; j < ileKolumn; j++)
		{
			C[i][j] = 0.0;
		}

	for (int i = 0; i < ileWierszy; i++)
		wynik[i] = 0.1;

	//wyrazy macierzy C
	for (int i = 0; i < ileWierszy; i++)
		for (int j = 0; j < ileKolumn; j++)
		{
			if (j == 0) C[i][0] = 0;
			else C[i][j] = - macierz[i][j] / macierz[i][0];
		}

	//wektor F
	double *F;
	F = new double[ileWierszy];

	for (int i = 0; i < ileWierszy; i++)
		F[i] = macierz[i][2] / macierz[i][0];

	double suma = 0;
	//iteracje
	for (int k = 0; k < 10; k++)
	{
		for (int i = 0; i < ileWierszy; i++)
		{
			suma = 0;
			for (int j = 0; j < ileKolumn; j++)
			{
				if (i + j <= ileWierszy) suma += C[i][j] * wynik[i + j];
				if (i - j >= 0) suma += C[i-j][j] * wynik[i - j];

			}
			wynik[i] = suma + F[i];
			//cout <<i<<"\t"<< suma << "\t" << F[i] << endl;

		}
	}

	//usuwanie zb�dnych zmiennych
	for (int i = 0; i<ileWierszy; i++)
		delete[] C[i];
	delete[] F;
}

void RozwRown::gauss(double **A, double *wynik, const int ileWierszy, const int ileKolumn )
{

	int n = ileWierszy;

	for (int i = 0; i<n; i++) {
		// Search for maximum in this column
		double maxEl = abs(A[i][i]);
		int maxRow = i;
		for (int k = i + 1; k<n; k++) {
			if (abs(A[k][i]) > maxEl) {
				maxEl = abs(A[k][i]);
				maxRow = k;
			}
		}

		// Swap maximum row with current row (column by column)
		for (int k = i; k<n + 1; k++) {
			double tmp = A[maxRow][k];
			A[maxRow][k] = A[i][k];
			A[i][k] = tmp;
		}

		// Make all rows below this one 0 in current column
		for (int k = i + 1; k<n; k++) {
			double c = -A[k][i] / A[i][i];
			for (int j = i; j<n + 1; j++) {
				if (i == j) {
					A[k][j] = 0;
				}
				else {
					A[k][j] += c * A[i][j];
				}
			}
		}
	}

	// Solve equation Ax=b for an upper triangular matrix A

	for (int i = n - 1; i >= 0; i--) {
		wynik[i] = A[i][n] / A[i][i];
		for (int k = i - 1; k >= 0; k--) {
			A[k][n] -= A[k][i] * wynik[i];
		}
	}
    for (int i = 0; i<ileWierszy; i++)
		delete[] A[i];

	delete[] A;
}
