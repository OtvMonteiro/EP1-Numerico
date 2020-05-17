#include <iostream>
#include <cmath>


using namespace std;

double u0(double x);
double f(double x,double t, double dx);
double g1(double t);
double g2(double t);
double u_esperado(double x, double t);
double r(double t);
double gh(double x, double dx);
double maximo(double a, double b);


int main()
{
	// Declaracao de variaveis
	int N ;
	int M ;
	double dt, dx, lambda;
	double T = 1.0;
	double t_decorrido = 0.0; 

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

//Ao inves de matrizes usar algumas variaveis que vao sendo atualizadas...
	// ------------------------------------------
	// Tratamento de dados
	double u_atual[N]; //Vetor de valores ao longo do eixo x no tempo k
	double u_proximo[N];//Vetor de valores ao longo do eixo x no tempo k+1



    //Condicao inVetor de erro ao longo do eixo x no icial
    for (int i=0; i!=N+1 ; i++){
        u_atual[i] = u0(i*dx);
    }


	// ---- Execucao
	 // Faz a conta ao longo de x para cada t
	for (int k=0; k!=M; k++){

		//Evolucao temporal das fronteiras
		u_atual[0] = g1(k*dt); //u0k
		u_atual[N] = g2(k*dt); //uNk

		// Calculo ao longo de x
		for (int i=1; i!=N; i++){
				//Equacao 11
			u_proximo[i] = u_atual[i] + dt*(  (u_atual[i-1] - (2.0*(u_atual[i])) + u_atual[i+1])/(dx*dx) + f(dx*i, dt*k, dx)  );

		}
		
		//Atualiza o valor de tempo decorrido
		t_decorrido += dt;
		//Quando se passam 0.1
		if(t_decorrido>=0.1){
				cout<<endl<< "Tempo: "<<k*dt<<endl<<endl;
				for(int i=0; i!=N+1; i++){
					cout<< i*dx <<" "<< u_proximo[i]<<endl;
				}
				cout<<endl<<endl;
				t_decorrido=0.0; //Zera o valor do tempo para reiniciar a contagem
		}	
            //Atualiza vetores no tempo apos calcular valores para todas as posicoes
        for(int i=1; i!=N ; i++){
            u_atual[i]=u_proximo[i];
        }

	}
	

}



double u0(double x){
	return 0.0;
}

double f(double x,double t, double dx){
	//f=r(t)*gh(x)
	return r(t)*gh(x,dx);
}  

double g1(double t){
	return 0.0;
}


double g2(double t){
	return 0.0;
}

double u_esperado(double x, double t){
	return 0.0;//Deve ser determinado numericamente
}


double r(double t){
	return 10000*(1 - (2*t*t));
}


double gh(double x, double dx){
	double p = 0.25;
	if(x>=p-dx/2 && x<=p+dx/2){return (1/dx);} 
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
