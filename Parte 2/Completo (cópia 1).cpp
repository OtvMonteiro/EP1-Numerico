
#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

double u0(double x);
double f(double x,double t, double dx, int k);
double g1(double t);
double g2(double t);
double gh(double x, double dx, double pk);
double r(double t);
double maximo(double a, double b);


char item='a';//Variavel universal para escolher detalhes dos itens

int main()
{
	// Declaracao de variaveis
	int N = 128;
	int M ;
	int nf = 1;
	double dt, dx, lambda;
	double T = 1.0;

    // Interface do usuario e entradas de dados
	cout << "Digite o item desejado ('a','b','c' ou 'd' sem apostrofes)" << endl;
	cin  >> item;

	//N e M
	if(item=='c'||item=='d'){
		cout << "Digite o valor desejado para 'N' e 'M'" << endl;
		cin  >> N;
	}
	//nf
	if(item=='b'){nf=4;}
	else if(item=='c'||item=='d'){nf=10;}
	
	  
	//Definicoes
	M=N;
	dt=T/M;
	dx=(1.0)/N;
	lambda = dt/(dx*dx);
	
	
	double u[N-1][nf]; //Vetor que contem valores da temperatura ao longo do tempo (para o atual k)



//Para os itens a e b usamos o metodo de cranck nicolson para determinar os valores de u em T
if(item=='a'||item=='b') {
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
    
    

// ---- Determinacao de uk(T,xi) por Cranck Nicolson
	for(int l=0; l!=M ; l++){	
		
		//Atualizando valores nas fronteiras	
		b[1]=u[1][k] + (lambda*(u[0][k]-2*u[1][k]+u[2][k])/2) + (dt*(f(dx*1,dt*l,dx,k)+f(dx*1,dt*(l+1),dx,k))/2) + (lambda*g1(dt*(l+1))/2);
		b[N-1]=u[N-1][k] + (lambda*(u[N-2][k]-2*u[N-1][k])/2) + (dt*(f(dx*(N-1),dt*l,dx,k)+f(dx*(N-1),dt*(l+1),dx,k))/2) + (lambda*g2(dt*(l+1))/2);
		//Atualizando os valores de b para esse loop
		for (int i=2; i!=N-1; i++){
			b[i]=u[i][k] + (lambda*(u[i-1][k]-2*u[i][k]+u[i+1][k])/2) + (dt*(f(dx*i,dt*l,dx,k)+f(dx*i,dt*(l+1),dx,k))/2);
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
}
	
	
	//Calculo do MMQ para encontrar as intensidades (dado uk(T, xi))
	
	
	//matriz de u, nfXnf
	double A[nf+1][nf+1];

	//matriz coluna de a, nfX1
	double a[nf];//saida
	double x[nf];//para execucao
	
	//matriz coluna de <uT,unf> , nfX1
	double bmmq[nf];
	
	//Montando a matriz de minimos quadrados
	for (int k=1;k<nf+1;k++){//k para linha
		//Como a matriz e'simetrica vamos calcular os elementos partindo da diagonal principal para cada linha 
		for (int j=k; j<nf+1; j++){	//j para coluna
			//Produto interno
			A[k][j] = 0.0;
			for (int i=1;i<N;i++){
				A[k][j] +=  (u[i][j])*(u[i][k]);
			}
		}	
	}	
	//Considerando a simetria, vamos reconstruir a matriz inteira
	for(int k=2; k<nf+1; k++){
		for(int i=1; i<k; i++){
			A[k][i] = A[i][k];
		}
	}		
	//Imprimindo a matriz
	cout<< "Matriz de MMQ"<<endl;
	for(int l=1;l<nf+1;l++){
		for(int c=1; c<nf+1;c++){
			cout<<A[l][c]<<" ";
		}
		cout<<endl;
	}
	//Construindo o vetor "bmmq" (<uT,uk>) ; bmmq[0] e' bmmq1, para k=1, valendo (<uT,u1>) ///Isso por causa de como e' escrito o metodo seguinte
	for (int k=1;k<nf+1;k++){
			if(item=='a'){bmmq[k-1] = 7*(A[k][1]);}
		else if(item=='b'){bmmq[k-1] = 2.3*(A[k][1]) + 3.7*(A[k][2])  + 0.3*(A[k][3])  +4.2*(A[k][4]);}
		else if(item=='c'){bmmq[k-1] = 7*(A[k][1]);}
		else if(item=='d'){bmmq[k-1] = 7*(A[k][1]);}
		cout<<"k= "<<k<<"   bmmq[k]="<<bmmq[k-1]<<endl;
	}
	
	
	//Resolucao do sistema linear de MMQ por decomposicao LDL*
		
	int Nf = nf+1; //Esse metodo de decomposicao vai de 0 a N-1, e as matrizes sao de 1 a nf 
	//Por isso vamos adaptar a matriz de decomposicao
	double Ammq[Nf-1][Nf-1];
	for(int u=0;u<Nf-1;u++){
		for(int v=0;v<Nf-1;v++){
			Ammq[u][v]=A[u+1][v+1];
		}
	}

    // Decomposicao em LDL*
	
	double Lmmq[Nf][Nf];
	double Dmmq[Nf-1];

    //Preparando as matrizes auxiliares
    for (int i=0; i!=Nf ; i++){
        Dmmq[i]=0.0; // Matriz Dmmq como diagonal principal, o resto da matriz tem 0 como valor
        for (int l=0; l!=Nf ; l++){
            if(i==l){Lmmq[i][l]=1.0;}//Diagonal principal de Lmmq
            else{Lmmq[i][l]=0.0;}//Limpando o resto
        }
    }

    //Executando

    Dmmq[0]=Ammq[0][0]; //Primeiro valor de Dmmq
    for(int j=1;j<Nf-1;j++){//Itera ate acabar a decomposicao
        for (int i=j;i<Nf-1;i++){//Roda cada linha

            Lmmq[i][j-1]=Ammq[i][j-1]/Ammq[j-1][j-1];//Ammqtualiza Cada coluna abaixo de Ammq[i][j-1]


            for(int l=j;l<Nf-1;l++){//Atualizacao em cada posicao (andando nas diferentes colunas numa linha)
                Ammq[i][l] = Ammq[i][l] - (Ammq[j-1][l])*(Lmmq[i][j-1]);
            }

        }

        Dmmq[j]=Ammq[j][j]; //Atualizando os valores de Dmmq apos cada grande iteracao

    }



	//Resolvendo Ly=b
	double y[Nf-1];


	for(int i=0; i<Nf-1; i++){
        double soma=0.0;//Variavel auxiliar que e' a soma de da multiplicacao dos elementos de Lmmq na linha i com Y
        for(int j=0;j<i;j++){//Calculo de soma, vai ate a posicao anterior `a diagonal principal de Lmmq
            soma = soma + y[j]*(Lmmq[i][j]);
        }
		y[i] = bmmq[i] - soma;//Calculo de y
	}



	//Resolvendo DLt*x=y
	double DLt[Nf-1][Nf-1];

    // Calculo da multiplicacao D*Lt
	for(int i=0; i<Nf-1; i++){
        for(int j=0;j<Nf-1;j++){
            DLt[i][j]=0.0; //Zerando posicoes antes de calcular
            for(int l=0;l<Nf-1;l++){//iteracao para cada multiplicacao a ser efetuada para um elemento final
                if(l==i){//Como Dmmq e' um vetor de diagonal principal, so' vamos multiplicar quando for pelas posicoes l=i
                  DLt[i][j] += Dmmq[i]*Lmmq[j][l];//Numa matriz normal usariamos D[i][l]*L[l][j], mas como D e' diagonal e queremos a transposta de L fica diferente
                }
            }
         }
    }

  //Solucionando DLt*x=y
   for(int i=Nf-1; i>-1; i--){
        double soma=0.0;//Variavel auxiliar que e' a soma de da multiplicacao dos elementos de L na linha i com Y
        for(int j=Nf-1;j>i;j--){//Calculo de soma, vai ate a posicao anterior `a diagonal principal de L
            soma = soma + x[j]*(DLt[i][j]);
        }
		x[i] = (y[i]-soma)/(DLt[i][i]);//Calculo de x
		a[i+1]=x[i];//Resultado, partindo de a1=a[1]
	}

		
	
	
	
	//SAIDAS


		//saida de verificacao das intensidades
		for (int i=1; i!=nf+1; i++){
           cout << "a"<< i <<" = " << a[i] << endl;
        }
	
       // Imprimir valores em comparacao (para t=k*dt) ao longo de x
    cout<<endl<< "Deseja imprimir valores em T ao longo de x ?[S/n]:"<<endl;
    char gui_choice1 = 'n';
    cin  >> gui_choice1;

    if(gui_choice1!='n'){
		for(int k=1;k<nf+1;k++){
        for (int i=0; i!=N+1; i++){
           cout << "u :"<< u[i][k] <<"      para x=" << dx*i <<"  e k="<<k<< endl;
        }
		}
    }



}



double u0(double x){
	return 0.0;
}

double f(double x,double t, double dx, int k){
	
	double p[20];//valor arbitrario para criar vetor grande suficiente
	
	//Setando valores de p, onde pk=p[k] de 1 a nf
	if(item=='a'){p[1]=0.35;}
	if(item=='b'){p[1]=0.15;p[2]=0.3;p[3]=0.7;p[4]=0.8;}
	if(item=='c'){ ifstream arquivo;
			arquivo.open("teste.txt");
			for(int j=1; j<10+1;j++){arquivo>>p[j];}
			arquivo.close();
	}	
	
		
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
