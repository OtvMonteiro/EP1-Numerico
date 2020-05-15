
#include <iostream>
#include <cmath>
using namespace std;

double u0(double x);
double funcaoA(double x,double t);
double g1(double t);
double g2(double t);
double u_esperado(double x, double t);

int main()
{
	// Declaracao de variaveis fixas
	int N ;
	int M ;
	double dt, dx, lambda;
	double T = 1.0;

		// Interface do usuario com entradas de dados
	cout << "Digite o valor desejado para 'N'" << endl;
	cin  >> N;
	cout << "Digite o valor desejado para 'M'" << endl;
	cin  >> M;


	// ------------------------------------------
	// Tratamento de dados
	double u[N+1][M+1]; //Ma pratica, mas funciona

    // Limpa matriz a ser utilizada
	for(int i=0; i!=N+1; i++){
        for(int k=0; k!=M+1; k++){
        u[i][k]=0.0;
        }
    }

	  //Definicoes de dt e dx
	dt=T/M;
	dx=(1.0)/N;

		// Imprime para o usuario dx, dt e lambda
	cout << "dx = " << dx << endl << "dt = " << dt << endl ;
	lambda = dt/(dx*dx);
	cout << "O lambda para as condicoes escolhidas e' : " << lambda << endl;

		//Condicao inicial
		for (int i=0; i!=N+1 ; i++){
			u[i][0] = u0(i*dx);
		}


	// ---- Execucao
	 // Faz a conta ao longo de x para cada t
	for (int k=0; k!=M; k++){

		//Evolucao temporal das fronteiras
		u[0][k] = g1(k*dt); //u0k
		u[N][k] = g2(k*dt); //uNk

		// Calculo ao longo de x
		for (int i=1; i!=N; i++){
				//Equacao 11
			u[i][k+1] =( u[i][k] + dt*( (( u[i-1][k] - (2.0*(u[i][k])) + u[i+1][k] )/(dx*dx))  +  funcaoA(dx*i, dt*k) )   );

		}
	}



	//teste
	cout  << "funcaoA(dx*8, dt*300) : " << funcaoA(dx*8, dt*300) << endl;

	// Imprimir valores em comparacao (para t=k*dt) ao longo de x
	for (int i=0; i!=N+1; i++){
		cout << "Encontrado:" << u[i][M] << "   Esperado:"<< u_esperado(i*dx, M*dt) << endl;
	}


	// Escolher um output para ser mostrado em comparacao com o esperado
	int i_out, k_out;

	cout << "Escolha um valor de 0 a N para mostrar:"<<endl;
	cin  >> i_out;
	cout << "Escolha um valor de 0 a M para mostrar:"<<endl;
	cin  >> k_out;
	cout << "Os resultados para u(x=" << i_out*dx <<" , t="<< k_out*dt<<"):"<<endl;
	cout << "Encontrado:" << u[i_out][k_out] << "   Esperado:"<< u_esperado(i_out*dx, k_out*dt) << endl;
	return 0;
}



double u0(double x){
	//return 0.0;//c.i.n.
	return x*x*(1.0-x);
}

double funcaoA(double x,double t){
    //return 10*(x*x)*(x-1) - 60*x*t + 20*t ;// Primeira funcao
	return (10.0*cos(10.0*t)*x*x*(1.0-x)*(1.0-x)) - ((1.0 + sin(10.0*t))*((12.0*x*x)-(12.0*x)+2.0));
}

double g1(double t){
	return 0.0*t;
}

double g2(double t){
	return 0.0*t;
}

double u_esperado(double x, double t){
	return (1.0 + sin( 10.0*t))*x*x*(1.0-x)*(1.0-x);
}
