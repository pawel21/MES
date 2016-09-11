#ifndef RozwRown_h
#define RozwRown_h

class RozwRown
{
private:


public:
	//zmienne
	//RozwRown(void);
	//~RozwRown(void);
	void cholesky(double **macierz, double *wynik, const int ileWierszy, const int ileKolumn = 2);
	void przemiatania(double **macierz, double *wynik, const int ileWierszy, const int ileKolumn = 2);
	void iteracyjna(double **macierz, double *wynik, const int ileWierszy, const int ileKolumn = 2);
	void gauss(double **A, double *wynik, const int ileWierszy, const int ileKolumn );
};
#endif // RozwRown_h
