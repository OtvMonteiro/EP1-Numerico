


#include <iostream>

int main(int argc, char **argv)
{
	int N=5;
	//cin >> N;
		
	double Ad[ ] = { 0, 1, 3, 5, 7};
	double As[ ] = { 0, 0, 2, 4, 6};
	
	
	
	
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

	std::cout << D[1] <<" "<< D[2] <<" "<< D[3] <<" "<< D[4] <<std::endl;
	std::cout << L[2] <<" "<< L[3] <<" "<< L[4] <<std::endl;
	
	return 0;
}

