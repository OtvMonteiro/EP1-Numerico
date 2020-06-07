
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
	double tal;//Erro de truncamento delimitado


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
	double u[N-1]; //Vetor que contem valores da temperatura ao longo do tempo
	
	double e_atual[N]; //Vetor de erro ao longo do eixo x no tempo k
	double e_proximo[N];//Vetor de erro ao longo do eixo x no tempo k+1
	double norma_e = 0.0; //Norma do erro em t=T-dt
	double norma_ekp = 0.0; //Norma do erro em t=T

    double truncIK = 0.0; //Erro local de truncamento, sera atualizado ao longo da execucao
    
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
		b[1]=u[1] + (lambda*(u[0]-2*u[1]+u[2])/2) + (dt*(f(dx*1,dt*k)+f(dx*1,dt*(k+1)))/2) + (lambda*g1(dt*(k+1))/2);
		b[N-1]=u[N-1] + (lambda*(u[N-2]-2*u[N-1])/2) + (dt*(f(dx*(N-1),dt*k)+f(dx*(N-1),dt*(k+1)))/2) + (lambda*g2(dt*(k+1))/2);
		//Atualizando os valores de b para esse loop
		for (int i=2; i!=N-1; i++){
			b[i]=u[i] + (lambda*(u[i-1]-2*u[i]+u[i+1])/2) + (dt*(f(dx*i,dt*k)+f(dx*i,dt*(k+1)))/2);
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

		//Calculando erros
		for (int i=1; i!=N; i++){
            e_proximo[i] = e_atual[i] + dt*(  ((e_atual[i-1] - (2.0*(e_atual[i])) + e_atual[i+1])/(dx*dx)) + truncIK  );
            norma_e = maximo(norma_e, e_proximo[i]);//Vai em busca do maior valor de erro		    
		    //Truncamento
		    truncIK = ((x[i] - u[i])/dt)  - ((u[i-1] - 2*u[i] + u[i+1])/(dx*dx))   -  f(i*dx, k*dt);
			tal=maximo(tal, truncIK);
		}

	
		//Setando o valor encontrado como proxima entrada (ou valor final)
		for (int i=1; i!=N; i++){
			u[i] = x[i];
			if(k<M-1){//atualiza todos menos o ultimo (usaremos e(k) e e(k+1))
				e_atual[i]=e_proximo[i];
			}
		}
	}		

	double normaET=0.0;
	for (int i=1; i!=N; i++){
		 //Calculo da norma do erro entre a solucao aproximada e a exata, em tk=T
		 normaET = maximo(normaET, (u_esperado(i*dx, M*dt) - u[i]) );
		 
		 //Equacao 19
		 norma_e = maximo(norma_e, e_atual[i]);//Vai em busca do maior valor de erro para t=T-dt
		
		 //Calculo de da norma de e em t=T
		 norma_ekp = maximo(norma_e, e_proximo[i]);//Vai em busca do maior valor de erro para t=T
	 
    }



   //SAIDAS
	cout<<endl<<endl;
	cout<<N<<" & "<<M<<" & "<<lambda<<" & "<<normaET<<" & "<<tal<<" & "<<norma_ekp<<" & "<<norma_e  + dt*tal<<"  \\\\ \\hline"<<endl<<endl;//saida para tabela latex
	cout<<"A norma do erro entre a solucao aproximada e a exata, em tk=T, e': " << normaET <<endl;
    cout<<"O erro de truncamento delimitado, T(dt,dx) e': "<< tal <<endl;
    cout<<"A norma do erro para T encontrada e': "<< norma_ekp <<endl;
	cout<<"A norma do erro para T esperada e': |e(i,k+1)|<="<<norma_e  + dt*tal << endl;//Equacao 22

	
	
	
	
	//interativas

       // Imprimir valores em comparacao (para t=k*dt) ao longo de x
    cout<<endl<< "Deseja imprimir valores em T ao longo de x para comparacao?[S/n]:"<<endl;
    char gui_choice1 = 'n';
    cin  >> gui_choice1;

    if(gui_choice1!='n'){
        for (int i=0; i!=N+1; i++){
           // cout<<u[i]<<" & "<<u_esperado(i*dx, M*dt)<<" & "<<abs(u[i]- u_esperado(i*dx, M*dt)) <<"  \\\\ \\hline"<<endl;//Saida para tabela em latex
            cout << "Encontrado:" << u[i] << "   Esperado:"<< u_esperado(i*dx, M*dt) << endl;
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
        cout << "Encontrado:" << u[i_out] << "   Esperado:"<< u_esperado(i_out*dx, M*dt) << endl;
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
