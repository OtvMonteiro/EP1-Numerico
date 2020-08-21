
#include <iostream>
#include <cmath>
using namespace std;

double u0(double x);
double f(double x,double t, double dx, int nf, int k);
double g1(double t);
double g2(double t);
double gh(double x, double dx, double pk);
double r(double t);
double maximo(double a, double b);


int main()
{
	// Declaracao de variaveis
	int N ;
	int M ;
	int nf;
	double dt, dx, lambda;
	double T = 1.0;

    // Interface do usuario e entradas de dados
	cout << "Digite o valor desejado para 'N' e 'M'" << endl;
	cin  >> N;
	M=N;
	cout << "Digite o valor desejado para 'nf'" << endl;
	cin  >> nf;
	  
	  //Definicoes de dt e dx
	dt=T/M;
	dx=(1.0)/N;

		// Imprime para o usuario dx, dt e lambda
	cout << "dx = " << dx << endl << "dt = " << dt << endl ;
	lambda = dt/(dx*dx);
	cout << "O lambda para as condicoes escolhidas e' : " << lambda << endl;

	double u[N-1][nf+1]; //Vetor que contem valores da temperatura ao longo do tempo (para o atual k)
	
	
//Diferentes valores de k, indo de 1 a nf
for(int k=1; k<nf+1; k++){
	// ------------------------------------------
	
    //Vetores utilizados na decomposicao de Ax=b, LDL*x=b
	double Ad[N-1];
	double As[N-1];
	double b[N-1] ;
	//Vetor que contem as respostas dos calculos de Ax=b
	double x[N-1];		
	//Matrizes de decomposicao LDL*
	double D[N-1];
	double L[N-1];
	
    
    //Valores iniciais para u
    for (int i=0; i!=N ; i++){
        u[i][k] = u0(i*dx);
    }
    
    

// ---- Execucao
	for(int l=0; l!=M ; l++){	
		
		//Atualizando valores nas fronteiras	
		b[1]=u[1][k] + (lambda*(u[0][k]-2*u[1][k]+u[2][k])/2) + (dt*(f(dx*1,dt*l,dx,nf,k)+f(dx*1,dt*(l+1),dx,nf,k))/2) + (lambda*g1(dt*(l+1))/2);
		b[N-1]=u[N-1][k] + (lambda*(u[N-2][k]-2*u[N-1][k])/2) + (dt*(f(dx*(N-1),dt*l,dx,nf,k)+f(dx*(N-1),dt*(l+1),dx,nf,k))/2) + (lambda*g2(dt*(l+1))/2);
		//Atualizando os valores de b para esse loop
		for (int i=2; i!=N-1; i++){
			b[i]=u[i][k] + (lambda*(u[i-1][k]-2*u[i][k]+u[i+1][k])/2) + (dt*(f(dx*i,dt*l,dx,nf,k)+f(dx*i,dt*(l+1),dx,nf,k))/2);
		}
		
		//Limpando matrizes L e D
		for (int i=0; i!=N ; i++){
			D[i]=0.0;
			L[i]=0.0;
		}
		
		//Setando valores de A para Cranl-Nicolson (metade do lambda utilizado no Metodo de Newton)
		for (int i=1; i!=N ; i++){
			Ad[i]=1+lambda;//diagonal principal
		}
		for (int i=2; i!=N ; i++){
			As[i]=-lambda/2;//diagonais secundarias
		}	
		//Decompondo matriz A 
		for(int i=1;i!=N;i++){
			L[i+1] = As[i+1]/Ad[i];
			Ad[i+1]= Ad[i+1] - (As[i+1]*L[i+1]);
			D[i]=Ad[i];
		}

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
		////Valores de x encontrados!

	
		//Setando o valor encontrado como proxima entrada (ou valor final)
		for (int i=1; i!=N; i++){
			u[i][k] = x[i];
		}
	}		

}
	//interativas

       // Imprimir valores em comparacao (para t=l*dt) ao longo de x
    cout<<endl<< "Deseja imprimir valores em T ao longo de x ?[S/n]:"<<endl;
    char gui_choice1 = 'n';
    cin  >> gui_choice1;

    if(gui_choice1!='n'){
		for(int k=1;k<nf+1;k++){
        for (int i=0; i!=N+1; i++){
           cout << "Encontrado :"<< ":" << u[i][k] <<"      para x=" << dx*i <<"  e k="<<k<< endl;
        }
		}
    }



}



double u0(double x){
	return 0.0;
}

double f(double x,double t, double dx, int nf, int k){
	
	double p[nf];
	
	//Setando valores de p, onde pk=p[k] de 1 a nf
	p[1]=0.35;
	p[2]=0.6;
			
	return r(t)*gh(x, dx, p[k]);
}  

double g1(double t){
	return 0.0;
}


double g2(double t){
	return 0.0;
}


double r(double t){
	return 10*(1 + cos(5*t));
}


double gh(double x, double dx, double pk){
	double h=dx;
	if(x>=pk-h/2 && x<=pk+h/2){return (1/h);} 
	else{return 0.0;}
	return 0.0;
}


//Calcula o maximo entre o absoluto de dois valores
double maximo(double a, double b){
    double c = abs(a);
    double d = abs(b);
    if(c<d){return d;}
    else{return c;}
}
