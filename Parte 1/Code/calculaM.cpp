
#include <iostream>
using namespace std;



int main()
{
	// Declaracao de variaveis fixas
	int N = 10;
	double M ;
	double dx, lambda;
	double T = 1;
	
		// Interface do usuario com entradas de dados
	cout << "Digite o valor desejado para 'N'" << endl;
	cin  >> N;
	cout << "Digite o valor desejado para 'lambda'" << endl;
	cin  >> lambda;
	
	// ------------------------------------------
	// Tratamento de dados
		
	dx=(1.0)/N;
	M = T/(dx*dx*lambda);	
	int out = (int) (M+0.5);
	cout << "O M a ser escolhido e': " << out;
	return 0;
}
