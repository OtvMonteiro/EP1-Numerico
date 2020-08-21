
#include <iostream>
#include <cmath>
using namespace std;

double u0(double x);
double f(double x,double t, double dx, int nf);
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

//Ao inves de matrizes usar algumas variaveis que vao sendo atualizadas...
	// ------------------------------------------
	// Tratamento de dados
	double u[N-1]; //Vetor que contem valores da temperatura ao longo do tempo
	
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
        u[i] = u0(i*dx);
    }
    
    

// ---- Execucao
	for(int k=0; k!=M ; k++){	
		
		//Atualizando valores nas fronteiras	
		b[1]=u[1] + (lambda*(u[0]-2*u[1]+u[2])/2) + (dt*(f(dx*1,dt*k,dx,nf)+f(dx*1,dt*(k+1),dx,nf))/2) + (lambda*g1(dt*(k+1))/2);
		b[N-1]=u[N-1] + (lambda*(u[N-2]-2*u[N-1])/2) + (dt*(f(dx*(N-1),dt*k,dx,nf)+f(dx*(N-1),dt*(k+1),dx,nf))/2) + (lambda*g2(dt*(k+1))/2);
		//Atualizando os valores de b para esse loop
		for (int i=2; i!=N-1; i++){
			b[i]=u[i] + (lambda*(u[i-1]-2*u[i]+u[i+1])/2) + (dt*(f(dx*i,dt*k,dx,nf)+f(dx*i,dt*(k+1),dx,nf))/2);
		}
		
		//Limpando matrizes L e D
		for (int i=0; i!=N ; i++){
			D[i]=0.0;
			L[i]=0.0;
		}
		
		//Setando valores de A para Crank-Nicolson (metade do lambda utilizado no Metodo de Newton)
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
			u[i] = x[i];
		}
	}		

	
	//interativas

       // Imprimir valores em comparacao (para t=k*dt) ao longo de x
    cout<<endl<< "Deseja imprimir valores em T ao longo de x ?[S/n]:"<<endl;
    char gui_choice1 = 'n';
    cin  >> gui_choice1;

    if(gui_choice1!='n'){
        for (int i=0; i!=N+1; i++){
           cout << "Encontrado :"<< ":" << u[i] <<"      para x=" << dx*i << endl;
        }
    }



}



double u0(double x){
	return 0.0;
}

double f(double x,double t, double dx, int nf){
	double soma=0.0;
	double p[nf];
	
	//Setando valores de p, onde pk=p[k] (descartamos o zero) de 1 a nf
	p[1]=0.35;
	
	for(int k=1;k<nf+1;k++){
		soma += gh(x, dx, p[k]);		
	}
	return r(t)*soma;
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
