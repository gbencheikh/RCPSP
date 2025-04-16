#include"Heuristique.h"
#include"Data.h"
#include<list>
#include <ilcplex/ilocplex.h>
#include<ilconcert/iloexpression.h>

using namespace :: std;

Heuristique :: Heuristique(string filemane){
	Data data(filemane);
	Njobs = data.getNjobs();
	Nresources = data.getNresources();

	nbE = new int[Njobs];
	E = new int*[Njobs];
	P = new int[Njobs];
	B = new int[Nresources];
	b = new int*[Njobs];
	for(int i = 0 ; i < Njobs ; i++){
	    nbE[i] = data.getnbEvelue(i);
	    E[i] = new int[nbE[i] ];
	    b[i] = new int[Nresources];
	    P[i] = data.getValueP(i);
	    for(int j = 0 ; j < nbE[i] ; j++){
	    	E[i][j] = data.getValueE(i,j);
	    }
	    for(int j = 0 ; j < Nresources ; j++){
	        b[i][j] = data.getValueb(i,j);
	    	B[j] = data.getValueB(j);
	    }
	}
	liste = new int[Njobs];
	S = new int[Njobs];
	for(int i = 0 ; i < Njobs ; i++){
		liste[i] = -1;
		S[i] = 0;
	}
	Z = 0 ;
}
void Heuristique :: creatList(){
	bool *T;

	T = new bool[Njobs];
	for(int i = 0 ; i < Njobs ; i++){
		T[i] = false;
		S[i] = -100;
	}

	T[0] = true;
	S[0] = 0;

	list<int> l;
	l.push_back(0);

	while(!l.empty()){
		list<int>::iterator max = l.begin();
		for(list<int>::iterator iter = l.begin() ; iter != l.end() ; ++iter){
			if(S[*iter] > S[*max])
				max = iter;
		}
		for(int i = 0 ; i < Njobs ; i++){
			if(arc(*max,i)==true){
				if(S[i] < S[*max] + P[*max]){
						S[i] = S[*max] + P[*max];
				}
				if(T[i] == false){
					T[i] = true;
					l.push_back(i);
				}
			}
		}
		T[*max] = false;
		l.erase(max);
	}
	for(int i = 0 ; i < Njobs ; i++){
		T[i] = false;
	}
	for(int i = 0 ; i < Njobs ; i++){
		int min = Njobs-1;
		for(int j = 0 ; j < Njobs ; j++){
			if(S[j] < S[min] && T[j] == false)
				min = j;
		}
		liste[i] = min;
		T[min] = true;
	}
}
void Heuristique :: serie(){
	int j = 0;
	int t = 0;
	int **r;
	r = new int*[Nresources];
	for(int k = 0 ; k < Nresources; k++){
		r[k] = new int[getsumP()];
		for(int t = 0 ; t < getsumP() ; t++){
			r[k][t] = 0;
		}
	}
	while(j != Njobs-1){
		int temp = 0;
		int job = liste[j];

		for(int k = 0 ; k < Nresources ; k++){
			for(int ro = t ; ro <= t + P[job] ; ro++)
				if(r[k][ro] + b[job][k] > B[k]){
					temp = 1;
				}
		}
		if(temp == 0){
			if(t > S[job])
				S[job] = t;
			for(int k = 0 ; k < Nresources ; k++){
				for(int ro = S[job] ; ro <= S[job] + P[job] ; ro++){
					r[k][ro] += b[job][k];
				}
			}
			j++;
			temp = 0;
			for(int i = 0 ; i < j ; i++){
				if(arc(liste[i],liste[j]) == true){
					t = S[liste[i] ]+P[liste[i] ];
					S[liste[j] ] = t;
					temp = 1;
				}
			}
			if(temp == 0)
				t = S[job];
		}
		else
			t++;
	}
	Z = S[Njobs-1];
}
bool Heuristique :: canBeThere(int i, int t, int **r){
	for(int to = t ; to < t+P[i]; to++){
		for(int k = 0 ; k < Nresources ; k++){
			if(b[i][k] + r[k][to] > B[k])
				return false;
		}
	}
	return true;
}
void Heuristique :: parallele(){
	int **r;
	r = new int*[Nresources];
	for(int k = 0 ; k < Nresources; k++){
		r[k] = new int[getsumP()];
		for(int t = 0 ; t < getsumP() ; t++){
			r[k][t] = 0;
		}
	}
	bool *T;
	T = new bool[Njobs];
	for(int i = 0 ; i < Njobs ; i++){
		T[i] = false;
	}
	T[0] = true;
	S[0] = 0;
	int t = 0;
	int nbJobs = 1;
	while(nbJobs != Njobs){
		bool found = false;
		int i = 0;
		while(found == false && i < Njobs){
			if(T[i] == false){
				bool temp = true;
				for(int j = 1 ; j < Njobs ; j++){
					if(i != j){
						if(arc(j,i) == true){
							if((T[j] == true && t < S[j]+P[j]) || T[j] == false)
								temp = false;
						}
					}
				}
				if(temp == true && canBeThere(i,t,r) == true){
					for(int k = 0 ; k < Nresources ; k++){
						for(int to = t ; to < t + P[i] ; to++){
							r[k][to] += b[i][k];
						}
					}
					S[i] = t;
					T[i] = true;
					nbJobs++;
					found = true;
				}
				else
					i++;
			}
			else
				i++;
		}
		if(found == false){
			t++;
		}
	}
	Z = S[Njobs-1];
}
bool Heuristique :: arc(int i, int j){ // check if arc (i,j) exists
	for(int k = 0 ; k < nbE[i] ; k++)
		if(E[i][k] == j)
			return true;
	return false;
}
int Heuristique :: getsumP() const {
	int sum = 0 ;
	for(int i = 0 ; i < Njobs ; i++)
		sum += P[i];
	return sum;
}
int Heuristique :: getmaxP() const {
	int max = 0;
	for(int i = 1 ; i < Njobs ; i++)
		if(P[i] > P[max])
			max = i;
	return P[max];
}
int Heuristique :: getminP() const {
	int min = 1;
	for(int i = 2 ; i < Njobs ; i++)
		if(P[min] > P[i])
			min = i;
	return P[min];
}
void Heuristique :: printList(){
	for(int i = 0 ; i < Njobs ; i++){
		cout << liste[i] << " ";
	}
}
void Heuristique :: printSolution(){
	for(int i = 0 ; i < Njobs ; i++){
			cout << S[i] << " ";
	}
}
void Heuristique :: printResource(){
}
void tri(int *L, int n){ // sort the completions times
	for(int i = 0 ; i < n-1; i++){
		int min = i;
		for(int j = i+1; j < n ; j++){
			if(L[j] > L[min]){
				min = j ;
			}
		}
		int temp = L[min];
		L[min] = L[i];
		L[i] = temp;
	}
}
bool Heuristique :: resource(){ // check if resources constraints are respected by the S
	int *L;
	L = new int[Njobs];
	for(int i = 0 ; i < Njobs ; i++){
		L[i] = S[i] + P[i];
	}
	tri(L,Njobs);
	for(int t = 0 ; t < Njobs ; t++){
		for(int k = 0 ; k < Nresources ; k++){
			int o = 0 ;
			for(int j = 0 ; j < Njobs ; j++){
				if(S[j] < L[t] && L[j] >= L[t] && b[j][k]>0){
					o += b[j][k];
				}
				if(o > B[k]){
					cout << k ;
					return false;

				}
			}
		}
	}
	return true;
}
