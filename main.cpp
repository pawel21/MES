#include"RozwRown.h"
#include<iostream>

void macierz_zerowa_tablica(double **A, int k, int n);
void create_dynamic_matrix(double **A, double Macierz[3][4], const int ileWierszy, const int ileKolumn);

int main()
{

  int wierszy=3;
  int kolumn=4;
  double **A = new double *[wierszy];
  for(int i = 0; i < wierszy; i++ )
  {
    A[i] = new double [kolumn];
  }

  double Macierz[3][4] = {
                          {3, 4, 5, 0} ,
                          {1, -10, 1, 0},
                          {1, 0, 1, 42.5}
                        };
  create_dynamic_matrix(A, Macierz, wierszy, kolumn);
  double *wynik;

  RozwRown rownanie;

  rownanie.gauss(A, wynik, 3, 4);
  std::cout<<wynik[0]<<"\n";
  std::cout<<wynik[1]<<"\n";
  std::cout<<wynik[2]<<"\n";
  return 0;
}

void create_dynamic_matrix(double **A, double Macierz[3][4], const int ileWierszy, const int ileKolumn)
{
  for(int i = 0; i<ileWierszy; i++)
  {
    for (int j = 0; j<ileKolumn; j++)
      {
        A[i][j] = Macierz[i][j];
        std::cout<<A[i][j]<<" ";
      }
    std::cout<<"\n";
  }
}
