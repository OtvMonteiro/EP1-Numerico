
#include <iostream>
#include <cmath>
using namespace std;

double u0(double x);
double f(double x,double t);
double g1(double t);
double g2(double t);
double u_esperado(double x, double t);
double maximo(double a, double b);


int main()
{
	// Declaracao de variaveis
	int N ;
	int M ;
	double dt, dx, lambda;
	double T = 1.0;
	double tal;


    // Interface do usuario e entradas de dados
	cout << "Digite o valor desejado para 'N'" << endl;
	cin  >> N;
	cout << "Digite o valor desejado para 'M'" << endl;
	cin  >> M;

	  //Definicoes de dt e dx
	dt=T/M;
	dx=(1.0)/N;

		// Imprime para o usuario dx, dt e lambda
	cout << "dx = " << dx << endl << "dt = " << dt << endl ;
	lambda = dt/(dx*dx);
	cout << "O lambda para as condicoes escolhidas e' : " << lambda << endl;


	// ------------------------------------------
	// Tratamento de dados
	double u[N+1][M+1]; //Matriz de temperatura
	double e[N+1][M+1]; // Matriz de erro
    double trunc[N+1][M+1]; //Matriz de erro local de truncamento
    double norma_e[M];// Vetor de normas do erro ao longo do tempo


    // Limpa matrizes a serem utilizadas
	for(int i=0; i!=N+1; i++){
        for(int k=0; k!=M+1; k++){
        u[i][k]=0.0;
        e[i][k]=0.0;
        trunc[i][k]=0.0;
        }
    }

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
			u[i][k+1] = u[i][k] + dt*(  (u[i-1][k] - (2.0*(u[i][k])) + u[i+1][k])/(dx*dx) + f(dx*i, dt*k)  );

                //Equacao 12 //Os valores de u foram calculados em lacos anteriores ou imediatamente antes
            trunc[i][k] = (u[i][k+1] - u[i][k])/dt  - (u[i-1][k] - 2*u[i][k] + u[i+1][k])/(dx*dx)   -  f(i*dx, k*dt);

                //Equacao 18
            e[i][k+1] = e[i][k] + dt*(  (e[i-1][k] - (2.0*(e[i][k])) + e[i+1][k])/(dx*dx) + trunc[i][k]  );
		}
	}


    // Calculo de normas dos erros

    for (int k=0; k!=M; k++){//O erro de truncamento nao foi calculado em T, pois nao e' usado em "e" e precisaria de u em T+dt
        for (int i=0; i!=N; i++){
            // Equacao 15
            tal=maximo(tal, trunc[i][k]);

            //Equacao 19
            norma_e[k+1]=maximo(norma_e[k+1], e[i][k+1]);//Vai de 1 a M em busca do maior valor de erro

        }
    }



  //SAIDAS


	//teste
	//cout  << "f(dx*8, dt*300) : " << f(dx*8, dt*300) << endl;


    cout<<endl<<"A norma do erro para T e': "<< norma_e[M]<<endl;


	//interativas

       // Imprimir valores em comparacao (para t=k*dt) ao longo de x
    cout<<endl<< "Deseja imprimir valores em T ao longo de x para comparacao?[S/n]:"<<endl;
    char gui_choice1 = 'n';
    cin  >> gui_choice1;

    if(gui_choice1!='n'){
        for (int i=0; i!=N+1; i++){
            cout << "Encontrado:" << u[i][M] << "   Esperado:"<< u_esperado(i*dx, M*dt) << endl;
        }
    }


       // Escolher um output para ser mostrado em comparacao com o esperado

    cout<<endl<< "Deseja escolher valores de N e M para avaliar os valores encontrado e esperado?[S/n]:"<<endl;
    char gui_choice2 = 'n';
    cin  >> gui_choice2;

    if(gui_choice2!='n'){

        int i_out, k_out;

        cout << "Escolha um valor de 0 a N para mostrar:"<<endl;
        cin  >> i_out;
        cout << "Escolha um valor de 0 a M para mostrar:"<<endl;
        cin  >> k_out;
        cout << "Os resultados para u(x=" << i_out*dx <<" , t="<< k_out*dt<<"):"<<endl;
        cout << "Encontrado:" << u[i_out][k_out] << "   Esperado:"<< u_esperado(i_out*dx, k_out*dt) << endl;
        return 0;
    }



}



double u0(double x){
	//return 0.0;//c.i.n.
	return x*x*(1.0-x);
}

double f(double x,double t){
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


//Calcula o maximo entre o absoluto de dois valores
double maximo(double a, double b){
    double eps=0.00000001;
    double c = abs(a);
    double d = abs(b);
    if(c-d<eps){return d;}
    else{return c;}
}
