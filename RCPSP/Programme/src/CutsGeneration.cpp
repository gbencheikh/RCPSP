#include"Data.h"
#include"CutsGeneration.h"
#include<list>
#include<stack>
#include <ilcplex/ilocplex.h>
#include<ilconcert/iloexpression.h>

using namespace :: std;

CutsGeneration :: CutsGeneration(string filemane){
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

}
/*
ILOUSERCUTCALLBACK1(CtCallback, IloExprArray, cuts){
	IloInt n = cuts.getSize();
	for(int i = 0 ; i < n ; i++){
		IloRange c;
		try{
			c = (cuts[i] >= 1);
			add(c).end();

		}catch(...){
			c.end();
			throw;
		}
	}
}

void makearg(list<list<int> > critique, list<list<int> > incompatible, int **a){
	int j = 0;
	for(list<list<int> >::iterator iterListe = critique.begin() ; iterListe != critique.end() ; ++iterListe){
		int i = 0;
		for(list<int> ::iterator iter = (*iterListe).begin() ; iter != (*iterListe).end() ; ++iter){
			a[j][i] = *iter;
			i++;
		}
		j++;
	}
}
*/
void CutsGeneration :: makecuts(const BoolVarMatrix& x, int **a, IloExprArray cuts){

	for(unsigned int i = 0 ; i < critique.size(); i++){
		cuts.add(x[a[i][0]][a[i][1]] + x[a[i][0]][a[i][2]] + x[a[i][1]][a[i][2]] + x[a[i][1]][a[i][0]] + x[a[i][2]][a[i][0]] + x[a[i][2]][a[i][1]] );
	}

	for(list<list<int> >::iterator iterListe = incompatible.begin(); iterListe != incompatible.end() ; ++iterListe){
		for(list<int>::iterator iter1 = (*iterListe).begin() ; iter1 != (*iterListe).end() ; ++iter1){
			for(list<int>::iterator iter2 = iter1 ; iter2 != (*iterListe).end() ; ++iter2){
				if(*iter1 != *iter2){
					cuts.add(x[*iter1][*iter2] + x[*iter2][*iter1]);
				}
			}
		}
	}
}

void CutsGeneration :: model(IloModel& model,BoolVarMatrix& x, IloFloatVarArray& S){

	IloEnv env = model.getEnv();

	/******** Allocation ******/

	x = BoolVarMatrix(env, Njobs);
	for (IloInt j = 0; j < Njobs; j++) {
		x[j] = IloBoolVarArray(env, Njobs);
		for (IloInt t = 0; t < Njobs; t++) {
			stringstream name;
			name << "x_" << j << "_" << t;
			x[j][t].setName(name.str().c_str());
		}
	}

	S = IloFloatVarArray(env,Njobs);
	for(IloInt i = 0 ; i < Njobs ; i++)
	    S[i] = IloFloatVar(env,0,10000);

	/******** Objective function ******/

	IloExpr v(env);
	v+=S[Njobs-1];

	/********* Constraint 1 ***********/
	for (IloInt i = 0; i < Njobs; i++) {
		for(IloInt j = 0 ; j < Njobs ; j++){
			if( arc(i,j) == true){
				IloExpr u(env);
				u += x[i][j];
				model.add(u == 1);
			}
		}
	}
	/********** Constraint 2 ***********/
	int m = M();
	for(IloInt i = 0 ; i < Njobs ; i++){
    	for(IloInt j = 0 ; j < Njobs ; j++){
    		IloExpr u(env);
    		u += S[j] - S[i] + m - (P[i] + m)*x[i][j];
    		model.add(u >= 0);
    	}
    }
	/********** Constraint 3 ***********/
	for(IloInt i = 0 ; i < Njobs ; i++){
		IloExpr u(env);
		u += x[i][i];
		model.add(u == 0);
	}
	/********** Constraint 4 ***********/

	int *a;
	a = new int[3];

	for(list<list<int> >::iterator iterListe = critique.begin() ; iterListe != critique.end() ; ++iterListe){
		int i = 0;
		for(list<int> ::iterator iter = (*iterListe).begin() ; iter != (*iterListe).end() ; ++iter){
			a[i] = *iter;
			i++;
		}
		IloExpr u(env);
		u += x[a[0]][a[1]] + x[a[0]][a[2]] + x[a[1]][a[2]] + x[a[1]][a[0]] + x[a[2]][a[0]] + x[a[2]][a[1]];
		model.add(u >= 1);
	}

	/***********************************/
	model.add(IloMinimize(env, v));
}
bool CutsGeneration :: arc(int i, int j){ // check if arc (i,j) exists
	for(int k = 0 ; k < nbE[i] ; k++)
		if(E[i][k] == j)
			return true;
	return false;
}

int CutsGeneration :: M() {
	return getsumP();
}
int CutsGeneration :: getsumP() const {
	int sum = 0 ;
	for(int i = 0 ; i < Njobs ; i++)
		sum += P[i];
	return sum;
}

double CutsGeneration :: callModel(BoolVarMatrix& x, IloFloatVarArray& S){
	IloEnv env;
	stringstream inf;
	try {
				IloModel mod(env);
				IloTimer Time(env);


				model(mod,x,S);

				IloCplex cplex(env);

				int **a;
				a = new int*[critique.size()];
				for(unsigned int i = 0 ; i < critique.size(); i++)
					a[i] = new int[3];
/*
				IloExprArray cuts(env);
				makearg(critique,incompatible,a);
				makecuts(x,a,cuts);
				cplex.use(CtCallback(env,cuts));
				cuts.end();
*/
				cplex.extract(mod);
				cplex.setParam(IloCplex::MIPInterval, 1);
				cplex.setParam(IloCplex::ClockType, 1);
				cplex.setParam(IloCplex::TiLim , 1800);//30 minutes

				cplex.solve();

				return cplex.getObjValue() ;


		} catch (IloException& e) {
				cout << "ERROR: " << e << endl;
		}
		env.end();
		return -1;
}

string CutsGeneration :: resolvingModel(){

    BoolVarMatrix x;
    IloFloatVarArray S;

    clock_t time = clock();
    double time_init = (double) time/CLOCKS_PER_SEC;

    double res = callModel(x,S);

    time = clock();
    double time_resolution = (double) time/CLOCKS_PER_SEC - time_init;

    stringstream chaine;
    if(res == -1)
    	chaine << "0;0";
    else {
		chaine << res << ";" << time_resolution;
    }
    return chaine.str();
}

void CutsGeneration :: incompatibleJobs(){
	for(int i = 1 ; i < Njobs-2 ; i++){
		list<int> liste;
		liste.push_back(i);
		for(int j = i+1 ; j < Njobs-1 ; j++){
			bool temp = true;
			for(list<int>::iterator iter = liste.begin();iter!=liste.end();++iter){
				if(canBeTogether(*iter,j) == true){
					temp = false;
				}
			}
			if(temp == true){
				liste.push_back(j);
			}
		}
		if(liste.size() > 1 )
			incompatible.push_back(liste);
	}
}

void CutsGeneration :: critiqueJobs(){
	for(int i = 0 ; i < Njobs - 3 ; i++){
		for(int j = i+1 ; j < Njobs-2 ; j++){
			if(canBeTogether(i,j) == true){
				for(int t = j+1 ; t < Njobs-1 ; t++ ){
					if(canBeTogether(i,t) == true && canBeTogether(j,t)==true && tupleCanBeTogether(i,j,t) == false){
						list<int> liste;
						liste.push_back(i);
						liste.push_back(j);
						liste.push_back(t);
						critique.push_back(liste);
					}
				}
			}
		}
	}
}
bool CutsGeneration :: tupleCanBeTogether(int i, int j, int t, int k){
	if(b[i][k] + b[j][k] + b[t][k] <= B[k])
		return true;
	return false;
}
bool CutsGeneration :: tupleCanBeTogether(int i, int j, int t){
	for(int k = 0 ; k < Nresources ; k++){
		if(tupleCanBeTogether(i,j,t,k) == false)
			return false;
	}
	return true;
}
bool CutsGeneration :: canBeTogether(int i, int j){
	for(int k = 0 ; k < Nresources ; k++){
		if(canBeTogether(i,j,k) == false)
			return false;
	}
	return true;
}
bool CutsGeneration :: canBeTogether(int i, int j, int k){
	if(b[i][k] + b[j][k] <= B[k])
		return true;
	return false;
}
void CutsGeneration :: printCritiqueJobs(){
	cout << "set of critical tasks : " << endl;
	for(list<list<int> >::iterator iterListe = critique.begin() ; iterListe != critique.end() ; ++iterListe){
		cout << "critique : ";
		for(list<int> ::iterator iter = (*iterListe).begin() ; iter != (*iterListe).end() ; ++iter){
			cout << *iter << " ";
		}
		cout << endl;
	}
}
void CutsGeneration :: printIncompatibleJobs(){
	cout << "set of pairs tasks that cannot be executed together " << endl;
	for(list<list<int> >::iterator iterListe = incompatible.begin() ; iterListe != critique.end() ; ++iterListe){
		cout << "set : ";
		for(list<int> ::iterator iter = (*iterListe).begin() ; iter != (*iterListe).end() ; ++iter){
			cout << *iter << " ";
		}
		cout << endl;
	}
}
