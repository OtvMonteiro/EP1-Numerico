


#include <iostream>

int main(int argc, char **argv)
{
	int N=5;
	//cin >> N;
		
	double Ad[ ] = { 0, 1, 3, 5, 7};
	double As[ ] = { 0, 0, 2, 4, 6};
	double b[ ]  = { 0, 0, 9, 12, 19};
	//x para esses valores: 2,-1,2,1 
	double x[N-1];
	
	
	
	double D[N-1];
	double L[N-1];
	
	for (int i=0; i!=N ; i++){
        D[i]=0.0;
		L[i]=0.0;
    }
	
	for(int i=1;i!=N;i++){
		L[i+1] = As[i+1]/Ad[i];
		Ad[i+1]= Ad[i+1] - (As[i+1]*L[i+1]);
		D[i]=Ad[i];
	}

	//std::cout << D[1] <<" "<< D[2] <<" "<< D[3] <<" "<< D[4] <<std::endl;
	//std::cout << L[2] <<" "<< L[3] <<" "<< L[4] <<std::endl;
	
	
	//Resolvendo Ly=b
	double y[N-1];
	
	y[0]=0.0;//Para que y[1]=b[1]
	for(int i=1; i!=N; i++){
		y[i] = b[i] - (y[i-1]*L[i]);
	}
	
	//Resolvendo DL*x=y
	double DLconj_d[N-1];
	double DLconj_s[N-1];
	
			//DL*
	for(int i=1; i!=N; i++){
		DLconj_d[i]=D[i];
		DLconj_s[i]=D[i-1]*L[i];
	}
	
		//x=...
	x[N-1]=y[N-1]/DLconj_d[N-1];//Define o valor de Xn-1, para ser capaz de calcular (de modo decrescente) os coeficientes	
	for(int i=N-2; i!=0; i--){
		x[i] = (y[i] - (x[i+1]*DLconj_s[i+1]))/DLconj_d[i];
	}
	//Agora temos os valores que solucionam Ax=b;
	
	std::cout << x[1] <<" "<< x[2] <<" "<< x[3] <<" "<< x[4] <<std::endl;
	
	
	
	return 0;
}

