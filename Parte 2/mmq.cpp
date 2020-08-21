
#include <iostream>



int main
{
	int nf;
	std::cin >> nf;
	
	//matriz de u, nfXnf
	double A[nf+1][nf+1];

	//matriz coluna de a, nfX1
	double x[nf];
	
	//matriz coluna de <uT,unf> , nfX1
	double b[nf];
	double uT; 
	
	//Montando a matriz de minimos quadrados
	for (int k=1;k<nf+1;k++){//k para linha
		//Como a matriz e'simetrica vamos calcular os elementos partindo da diagonal principal para cada linha 
		for (int j=k; j<nf+1; j++){	//j para coluna
			//Produto interno
			A[k][j] = 0.0;
			for (int i=1;i<N;i++){
				A[k][j] +=  u[i] * u[j];
			}
		}	
	}	
	//Considerando a simetria, vamos reconstruir a matriz inteira
	for(int k=2; k<nf+1; k++){
		for(int i=1; i<k; i++){
			A[k][i] = A[i][k];
		}
	}		
	//Construindo o vetor "b" (<uT,uk>)
	for (int k=1;k<nf+1;k++){
		b[i] = 7*(A[k][1]);
	}
	
	
	
	
	
	return 0;
}

