


#include <iostream>

int main(int argc, char **argv)
{

    int N=5;

	double x[N-1];
    double A[N][N] = {
		{1, 2, 3, 4} ,
        {2, 5, 6, 7} ,
        {3, 6, 8, 9} ,
        {5, 7, 9, 10}
};

    double b[4]  = { 8, 18, 21, 21};
	//x para esses valores: -2, 3, 0, 1



	double L[N][N];
	double D[N-1];




    // Decomposicao em LDL*

    //Preparando as matrizes
    for (int i=0; i!=N ; i++){
        D[i]=0.0; // Matriz D como diagonal principal, o resto da matriz tem 0 como valor
        for (int k=0; k!=N ; k++){
            if(i==k){L[i][k]=1.0;}//Diagonal principal de L
            else{L[i][k]=0.0;}//Limpando o resto
        }
    }

    //Executando

    D[0]=A[0][0]; //Primeiro valor de D
    for(int j=1;j<N-1;j++){//Itera ate acabar a decomposicao
        for (int i=j;i<N-1;i++){//Roda cada linha

            L[i][j-1]=A[i][j-1]/A[j-1][j-1];//Atualiza Cada coluna abaixo de A[i][j-1]


            for(int k=j;k<N-1;k++){//Atualizacao em cada posicao (andando nas diferentes colunas numa linha)
                A[i][k] = A[i][k] - (A[j-1][k])*(L[i][j-1]);
            }

        }

        D[j]=A[j][j]; //Atualizando os valores de D apos cada grande iteracao

    }



	//Resolvendo Ly=b
	double y[N-1];


	for(int i=0; i<N-1; i++){
        double soma=0.0;//Variavel auxiliar que e' a soma de da multiplicacao dos elementos de L na linha i com Y
        for(int j=0;j<i;j++){//Calculo de soma, vai ate a posicao anterior `a diagonal principal de L
            soma = soma + y[j]*(L[i][j]);
        }
		y[i] = b[i] - soma;//Calculo de y
	}



	//Resolvendo DLt*x=y
	double DLt[N-1][N-1];

    // Calculo da multiplicacao D*Lt
	for(int i=0; i<N-1; i++){
        for(int j=0;j<N-1;j++){
            DLt[i][j]=0.0; //Zerando posicoes antes de calcular
            for(int k=0;k<N-1;k++){//iteracao para cada multiplicacao a ser efetuada para um elemento final
                if(k==i){//Como D e' um vetor de diagonal principal, so' vamos multiplicar quando for pelas posicoes k=i
                  DLt[i][j] += D[i]*L[j][k];//Numa matriz normal usariamos D[i][k]*L[k][j], mas como D e' diagonal e queremos a transposta de L fica diferente
                }
            }
         }
    }

  //Solucionando DLt*x=y
   for(int i=N-1; i>-1; i--){
        double soma=0.0;//Variavel auxiliar que e' a soma de da multiplicacao dos elementos de L na linha i com Y
        for(int j=N-1;j>i;j--){//Calculo de soma, vai ate a posicao anterior `a diagonal principal de L
            soma = soma + x[j]*(DLt[i][j]);
        }
		x[i] = (y[i]-soma)/(DLt[i][i]);//Calculo de x
	}



	//Agora temos os valores que solucionam Ax=b;

	std::cout << x[0] <<" "<< x[1] <<" "<< x[2] <<" "<< x[3] <<std::endl;







	return 0;
}

