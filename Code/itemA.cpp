
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

//Ao inves de matrizes usar algumas variaveis que vao sendo atualizadas...
	// ------------------------------------------
	// Tratamento de dados
	double u_atual[N]; //Vetor de valores ao longo do eixo x no tempo k
	double u_proximo[N];//Vetor de valores ao longo do eixo x no tempo k+1

	double e_atual[N]; //Vetor de erro ao longo do eixo x no tempo k
	double e_proximo[N];//Vetor de erro ao longo do eixo x no tempo k+1

    double truncIK = 0.0; //Erro local de truncamento, sera atualizado ao longo da execucao
    double norma_e = 0.0; //Norma do erro em t=T-dt
	double norma_ekp = 0.0; //Norma do erro em t=T


    //Condicao inVetor de erro ao longo do eixo x no icial
    for (int i=0; i!=N+1 ; i++){
        u_atual[i] = u0(i*dx);
        e_atual[i] = 0.0;
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
			u_proximo[i] = u_atual[i] + dt*(  (u_atual[i-1] - (2.0*(u_atual[i])) + u_atual[i+1])/(dx*dx) + f(dx*i, dt*k)  );

                //Equacao 12 //Os valores de u foram calculados em lacos anteriores ou imediatamente antes
            truncIK = ((u_proximo[i] - u_atual[i])/dt)  - ((u_atual[i-1] - 2*u_atual[i] + u_atual[i+1])/(dx*dx))   -  f(i*dx, k*dt);

                //Equacao 15
            tal=maximo(tal, truncIK);

                //Equacao 18
            e_proximo[i] = e_atual[i] + dt*(  ((e_atual[i-1] - (2.0*(e_atual[i])) + e_atual[i+1])/(dx*dx)) + truncIK  );

            
		}
            //Atualiza vetores no tempo apos calcular valores para todas as posicoes
        for(int i=0; i!=N+1 ; i++){
            u_atual[i]=u_proximo[i];
			if(k<M-1){//atualiza todos menos o ultimo (usaremos e(k) e e(k+1))
				e_atual[i]=e_proximo[i];
			}
        }

	}
	
	double normaET=0.0;
	for (int i=0; i!=N+1; i++){
	 //Calculo da norma do erro entre a solucao aproximada e a exata, em tk=T
     normaET = maximo(normaET, (u_esperado(i*dx, M*dt) - u_proximo[i]) );
     
     //Equacao 19
     norma_e = maximo(norma_e, e_atual[i]);//Vai em busca do maior valor de erro para t=T-dt
	
	 //Calculo de da norma de e em t=T
     norma_ekp = maximo(norma_e, e_proximo[i]);//Vai em busca do maior valor de erro para t=T
	 
    }


  //SAIDAS
	cout<<endl<<endl;
	cout<<"A norma do erro entre a solucao aproximada e a exata, em tk=T, e': " << normaET <<endl;
    cout<<"O erro de truncamento delimitado, T(dt,dx) e': "<< tal <<endl;
    cout<<"A norma do erro para T encontrada e': "<< norma_ekp <<endl;
	cout<<"A norma do erro para T esperada e': |e(i,k+1)|<="<<((1-2*lambda)+2*lambda)*norma_e  + dt*tal << endl;//Equacao 22
		
	
	
	
	//interativas

       // Imprimir valores em comparacao (para t=k*dt) ao longo de x
    cout<<endl<< "Deseja imprimir valores em T ao longo de x para comparacao?[S/n]:"<<endl;
    char gui_choice1 = 'n';
    cin  >> gui_choice1;

    if(gui_choice1!='n'){
        for (int i=0; i!=N+1; i++){
            cout<<u_proximo[i]<<" & "<<u_esperado(i*dx, M*dt)<<" & "<<abs(u_proximo[i]- u_esperado(i*dx, M*dt)) <<"  \\\\ \\hline"<<endl;//Saida para tabela em latex
			//cout << "Encontrado:" << u_proximo[i] << "   Esperado:"<< u_esperado(i*dx, M*dt) << endl;
        }
    }


       // Escolher um output para ser mostrado em comparacao com o esperado

    cout<<endl<< "Deseja escolher o valor de N para avaliar o valor encontrado e esperado?[S/n]:"<<endl;
    char gui_choice2 = 'n';
    cin  >> gui_choice2;

    if(gui_choice2!='n'){

        int i_out;

        cout << "Escolha um valor de 0 a N para mostrar:"<<endl;
        cin  >> i_out;
        cout << "Os resultados para u(x=" << i_out*dx <<"):"<<endl;
        cout << "Encontrado:" << u_proximo[i_out] << "   Esperado:"<< u_esperado(i_out*dx, M*dt) << endl;
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
    double c = abs(a);
    double d = abs(b);
    if(c<d){return d;}
    else{return c;}
}
