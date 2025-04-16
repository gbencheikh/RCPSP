#include"Resolution.h"
#include"Heuristique.h"
#include"Data.h"
#include<ctime>
#include<list>

using namespace std;

int min(int a , int b){
	if(a < b) return a;
	return b;
}
int max(int a, int b){
	if(a>b) return a;
	return b;
}

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

void Resolution :: makecuts(const BoolVarMatrix& x, int **a, IloExprArray cuts){
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

Resolution :: Resolution(string filename){
	Data data(filename);
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
	ES = new int[Njobs];
	LS = new int[Njobs];
	for(int i = 0 ; i < Njobs ; i++){
		ES[i] = -1;
		LS[i] = getsumP();
	}
	earliestStarting();
	latestStarting(filename);
}
/*********************************************/
/*                                           */
/*                                           */
/*            Flood's formulation            */
/*                                           */
/*                                           */
/*********************************************/

void Resolution :: ModelFlot(IloModel& model,BoolVarMatrix& x, IloFloatVarArray& S, IntVarMatrix3& f){

	IloEnv env = model.getEnv();

	/******** Allocation ******/

	x = BoolVarMatrix(env, Njobs);

	for (IloInt j = 0; j < Njobs; j++) {
		x[j] = IloBoolVarArray(env, Njobs);
	}

	f = IntVarMatrix3(env, Njobs);

	for (int i = 0; i < Njobs ; i++) {
			f[i] = IloArray<IloIntVarArray> (env, Njobs);
			for (int j = 0; j < Njobs; j++) {
				f[i][j] = IloIntVarArray(env, Nresources);
				for(IloInt k = 0 ; k < Nresources ; k++){
					f[i][j][k] = IloIntVar(env,0,10000);
				}
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
	/********** Constraint 2***********/
	for(IloInt i = 0 ; i < Njobs-1 ; i++){
		for(IloInt j = i+1 ; j < Njobs ; j++){
		    IloExpr u(env);
		    u += x[i][j] + x[j][i];
		    model.add(u <= 1);
		}
	}
	/********** Constraint 3 ***********/

	for(IloInt i = 0 ; i < Njobs ; i++){
		for(IloInt j = 0 ; j < Njobs ; j++){
			for(IloInt k = 0 ; k < Njobs ; k++){
				IloExpr u(env);
				u += -x[i][k] + x[j][k] + x[i][j];
				model.add(u <= 1);
			}
		}
	}
	/********** Constraint 4 ***********/
	for(IloInt i = 0 ; i < Njobs ; i++){
    	for(IloInt j = 0 ; j < Njobs ; j++){
    		IloExpr u(env);
    		int M = calculM(i,j);
    		//int M = getsumP();
    		u += S[j] - S[i] + M - (P[i] + M)*x[i][j];
    		model.add(u >= 0);
    	}
    }
    /********** Constraint 5 ***********/
    for(IloInt i = 0 ; i < Njobs-1 ; i++){
    	for(IloInt k = 0 ; k < Nresources ; k++){
    		IloExpr u(env);
    		for(IloInt j = 1 ; j < Njobs ; j++){
    			u += f[i][j][k];
    		}
    		model.add(u==b[i][k]);
    	}
    }
    /********** Constraint 6 ***********/
    for(IloInt j=1 ; j < Njobs ; j++){
    	for(IloInt k = 0 ; k < Nresources ; k++){
    		IloExpr u(env);
    		for(IloInt i = 0 ; i < Njobs-1 ; i++){
    			u += f[i][j][k];
    		}
    		model.add(u==b[j][k]);
    	}
    }
    /********** Constraint 7 ***********/
    for(IloInt i = 1 ; i < Njobs-1 ; i++){
    	for(IloInt j = 1 ; j < Njobs-1 ; j++){
    		for(IloInt k = 0 ; k < Nresources ; k++){
    			IloExpr u(env);
    			u += f[i][j][k] - min(b[i][k],b[j][k])*x[i][j];
    			model.add(u <= 0);
    		}
    	}
    }
    /***********************************/
	model.add(IloMinimize(env, v));
}
void Resolution :: ModelFlotRelaxed(IloModel& model,FloatVarMatrix& x, IloFloatVarArray& S, IntVarMatrix3& f){

	IloEnv env = model.getEnv();

	/******** Allocation ******/

	x = FloatVarMatrix(env, Njobs);
	for(int i = 0 ; i < Njobs ; i++){
		x[i] = IloFloatVarArray(env,Njobs);
		for(int j = 0 ; j < Njobs ; j++){
			x[i][j] = IloFloatVar(env,0,1);
		}
	}

	f = IntVarMatrix3(env, Njobs);

	for (int i = 0; i < Njobs ; i++) {
			f[i] = IloArray<IloIntVarArray> (env, Njobs);
			for (int j = 0; j < Njobs; j++) {
				f[i][j] = IloIntVarArray(env, Nresources);
				for(IloInt k = 0 ; k < Nresources ; k++){
					f[i][j][k] = IloIntVar(env,0,10000);
				}
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
	/********** Constraint 2***********/
	for(IloInt i = 0 ; i < Njobs-1 ; i++){
		for(IloInt j = i+1 ; j < Njobs ; j++){
		    IloExpr u(env);
		    u += x[i][j] + x[j][i];
		    model.add(u <= 1);
		}
	}
	/********** Constraint 3 ***********/

	for(IloInt i = 0 ; i < Njobs ; i++){
		for(IloInt j = 0 ; j < Njobs ; j++){
			for(IloInt k = 0 ; k < Njobs ; k++){
				IloExpr u(env);
				u += -x[i][k] + x[j][k] + x[i][j];
				model.add(u <= 1);
			}
		}
	}
	/********** Constraint 4 ***********/
	for(IloInt i = 0 ; i < Njobs ; i++){
    	for(IloInt j = 0 ; j < Njobs ; j++){
    		IloExpr u(env);
    		int M = getsumP();
    		u += S[j] - S[i] + M - (P[i] + M)*x[i][j];
    		model.add(u >= 0);
    	}
    }
    /********** Constraint 5 ***********/
    for(IloInt i = 0 ; i < Njobs-1 ; i++){
    	for(IloInt k = 0 ; k < Nresources ; k++){
    		IloExpr u(env);
    		for(IloInt j = 1 ; j < Njobs ; j++){
    			u += f[i][j][k];
    		}
    		model.add(u==b[i][k]);
    	}
    }
    /********** Constraint 6 ***********/
    for(IloInt j=1 ; j < Njobs ; j++){
    	for(IloInt k = 0 ; k < Nresources ; k++){
    		IloExpr u(env);
    		for(IloInt i = 0 ; i < Njobs-1 ; i++){
    			u += f[i][j][k];
    		}
    		model.add(u==b[j][k]);
    	}
    }
    /********** Constraint 7 ***********/
    for(IloInt i = 1 ; i < Njobs-1 ; i++){
    	for(IloInt j = 1 ; j < Njobs-1 ; j++){
    		for(IloInt k = 0 ; k < Nresources ; k++){
    			IloExpr u(env);
    			u += f[i][j][k] - min(b[i][k],b[j][k])*x[i][j];
    			model.add(u <= 0);
    		}
    	}
    }
    /***********************************/
	model.add(IloMinimize(env, v));
}
double Resolution :: callModelFlot(BoolVarMatrix& x, IloFloatVarArray& S, IntVarMatrix3& f){
	IloEnv env;
	stringstream inf;
	try {
				IloModel model(env);
				IloTimer Time(env);

				ModelFlot(model,x,S,f);

				IloCplex cplex(env);

				int **a;
				a = new int*[critique.size()];
				for(unsigned int i = 0 ; i < critique.size(); i++)
					a[i] = new int[3];

				IloExprArray cuts(env);
				makearg(critique,incompatible,a);
				makecuts(x,a,cuts);
				cplex.use(CtCallback(env,cuts));
				cuts.end();

				cplex.extract(model);
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

string Resolution :: resolvingFlot(){

    BoolVarMatrix x;
    IntVarMatrix3 f;
    IloFloatVarArray S;
    for(IloInt k = 0 ; k < Nresources ; k++){
        b[0][k] = B[k];
        b[Njobs-1][k] = B[k];
    }
    clock_t time = clock();
    double time_init = (double) time/CLOCKS_PER_SEC;

    double res = callModelFlot(x,S,f);

    time = clock();
    double time_resolution = (double) time/CLOCKS_PER_SEC - time_init;

    stringstream chaine;
    if(res == -1)
    	chaine << "0,0";
    else {
		chaine << res << ";" << time_resolution;
    }
    return chaine.str();
}

/*********************************************/
/*                                           */
/*                                           */
/*        Time indexed formulation           */
/*                                           */
/*                                           */
/*********************************************/

void Resolution :: ModelTimeIndexed(IloModel& model,BoolVarMatrix& y, int T){
	IloEnv env = model.getEnv();

	/*********** Allocation ***********/

	y = BoolVarMatrix(env, Njobs);

	for (IloInt j = 0; j < Njobs; j++) {
		y[j] = IloBoolVarArray(env, T);
	}

	/******** Objective function ******/
	IloExpr v(env);
	for(IloInt t = 0 ; t < T ; t++){
		v += t*y[Njobs-1][t];
	}

	/********* Constraint 1 ***********/
	for(IloInt j = 0 ; j < Njobs ; j++){
		IloExpr u(env);
		for(IloInt t = 0 ; t < T ; t++){
			u += y[j][t];
		}
		model.add(u == 1);
	}
	/***** Precedence constraint ******/
	for(IloInt i = 0 ; i < Njobs ; i++){
		for(IloInt j = 0 ; j < Njobs ; j++){
			if( arc(i,j) == 1 ){
				IloExpr u(env);
				for(IloInt t = 0 ; t < T ; t++){
					u += t*(y[j][t]-y[i][t]);
				}
				model.add(u >= P[i]);
			}
		}
	}
	/****** Resource constraint ******/
	for(IloInt k = 0 ; k < Nresources ; k++){
		for(IloInt t = 0 ; t < T ; t++){
			IloExpr u(env);
			for(IloInt j = 0 ; j < Njobs ; j++){
				if(t-P[j]+1 >= 0){
					IloExpr w(env);
					for(IloInt r = t-P[j]+1 ; r <= t ; r++){
						w += y[j][r];
					}
					u += b[j][k] * w;
				}
			}
			model.add(u <= B[k]);
		}
	}

	/**********************************/
	/*        Valid inequality        */
	/**********************************/
	for(IloInt j = 0 ; j < Njobs ; j++){
		for(int t = 0 ; t < ES[j] ; t++){
			IloExpr u(env);
			u = y[j][t];
			model.add(u == 0);
		}
		for(int t = LS[j]+1 ; t < T ; t++){
			IloExpr u(env);
			u = y[j][t];
			model.add(u == 0);

		}
	}

	model.add(IloMinimize(env, v));
}
void Resolution :: ModelTimeIndexedRelaxed(IloModel& model,FloatVarMatrix& y, int T){
	IloEnv env = model.getEnv();

		/*********** Allocation ***********/

		y = FloatVarMatrix(env, Njobs);

		for (IloInt j = 0; j < Njobs; j++) {
			y[j] = IloFloatVarArray(env, T);
			for(IloInt t = 0 ; t < T ; t++){
				y[j][t] = IloFloatVar(env,0,1);
			}
		}

		/******** Objective function ******/
		IloExpr v(env);
		for(IloInt t = 0 ; t < T ; t++){
			v += t*y[Njobs-1][t];
		}

		/********* Constraint 1 ***********/
		for(IloInt j = 0 ; j < Njobs ; j++){
			IloExpr u(env);
			for(IloInt t = 0 ; t < T ; t++){
				u += y[j][t];
			}
			model.add(u == 1);
		}
		/***** Precedence constraint ******/
		for(IloInt i = 0 ; i < Njobs ; i++){
			for(IloInt j = 0 ; j < Njobs ; j++){
				if( arc(i,j) == 1 ){
					IloExpr u(env);
					for(IloInt t = 0 ; t < T ; t++){
						u += t*(y[j][t]-y[i][t]);
					}
					model.add(u >= P[i]);
				}
			}
		}
		/****** Resource constraint ******/
		for(IloInt k = 0 ; k < Nresources ; k++){
			for(IloInt t = 0 ; t < T ; t++){
				IloExpr u(env);
				for(IloInt j = 0 ; j < Njobs ; j++){
					if(t-P[j]+1 >= 0){
						IloExpr w(env);
						for(IloInt r = t-P[j]+1 ; r <= t ; r++){
							w += y[j][r];
						}
						u += b[j][k] * w;
					}
				}
				model.add(u <= B[k]);
			}
		}


		model.add(IloMinimize(env, v));
}
double Resolution :: callModelTimeIndexed(BoolVarMatrix& y){
	IloEnv env;
	try {
			IloModel model(env);
			IloTimer Time(env);
			int T = LS[Njobs-1];
			ModelTimeIndexed(model,y,T);

			IloCplex cplex(env);

			cplex.extract(model);
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
string Resolution :: resolvingTimeIndexed(){

    BoolVarMatrix y;
    IntVarMatrix3 f;
    IloFloatVarArray S;

    clock_t time = clock();
    double time_init = (double) time/CLOCKS_PER_SEC;

    double res = callModelTimeIndexed(y);

    time = clock();
    double time_resolution = (double) time/CLOCKS_PER_SEC - time_init;

    stringstream chaine;
    if(res == -1) // If solution needs more than 30min to be found
    	chaine << "0;0";
    else {
		chaine << res << ";" << time_resolution;
    }
    return chaine.str();
}

/*********************************************/
/*                                           */
/*                                           */
/*            Calculate ES_{i}               */
/*                                           */
/*                                           */
/*********************************************/

void Resolution::earliestStarting(){

	bool *T;
	T = new bool[Njobs];
	for(int i = 0 ; i < Njobs ; i++){
		T[i] = false;
	}
	T[0] = true;
	ES[0] = 0;
	list<int> l;
	l.push_back(0);
	while(!l.empty()){
		list<int>::iterator max = l.begin();
		for(list<int>::iterator iter = l.begin() ; iter != l.end() ; ++iter){
			if(ES[*iter] > ES[*max])
				max = iter;
		}
		for(int i = 0 ; i < Njobs ; i++){
			if(arc(*max,i)==true){
				if(ES[i] < ES[*max] + P[*max]){
					ES[i] = ES[*max] + P[*max];
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
}
/*********************************************/
/*                                           */
/*                                           */
/*            calculate LS_{i}               */
/*         using parallel algorithm          */
/*                                           */
/*********************************************/

void Resolution :: latestStarting(string namefile){
	Heuristique hrs(namefile);
	hrs.creatList();
	hrs.parallele();

	bool *T;
	T = new bool[Njobs];
	for(int i = 0 ; i < Njobs ; i++){
		T[i] = false;
	}
	cout << "Z = " << hrs.getZ() << endl;

	T[Njobs-1] = true;
	LS[Njobs-1] = hrs.getZ();
	list<int> l;
	l.push_back(Njobs-1);
	while(!l.empty()){
		list<int>::iterator min = l.begin();
		for(list<int>::iterator iter = l.begin() ; iter != l.end() ; ++iter){
			if(LS[*iter] < LS[*min])
				min = iter;
		}
		for(int i = 0 ; i < Njobs ; i++){
			if(arc(i,*min)==true){
				if(LS[i] > LS[*min] - P[i]){
					LS[i] = LS[*min] - P[i];
				}
				if(T[i] == false){
					T[i] = true;
					l.push_back(i);
				}
			}
		}
		l.erase(min);
	}
}

/*  return 1 if i must be before j  */

bool Resolution :: arc(int i, int j){ // check if arc (i,j) exists
	for(int k = 0 ; k < nbE[i] ; k++)
		if(E[i][k] == j)
			return true;
	return false;
}
/*  calculate horizon  */
int Resolution :: getsumP() const {
	int sum = 0 ;
	for(int i = 0 ; i < Njobs ; i++)
		sum += P[i];
	return sum;
}
int Resolution :: getmaxP() const {
	int max = 0;
	for(int i = 1 ; i < Njobs ; i++)
		if(P[i] > P[max])
			max = i;
	return P[max];
}
int Resolution :: getminP() const {
	int min = 1;
	for(int i = 2 ; i < Njobs ; i++)
		if(P[min] > P[i])
			min = i;
	return P[min];
}
int Resolution :: calculM(int i, int j) {
	return max(LS[i],LS[j]) - min(P[i],P[j]);
}
bool Resolution :: precedenceFeasibility(bool **x){ // check if precedence constraints are respected by x

	for(int i = 0 ; i < Njobs ; i++){
		for(int j = 0 ; j < Njobs ; j++){
			if(arc(i,j) == true && x[j][i] == true)
				return false;
		}
	}
	return true;
}
string Resolution :: printESandLS(){
	stringstream chaine;
	for(int i = 0 ; i < Njobs ; i++){
		chaine << i << " : " << ES[i] << ";" << LS[i] << endl;
	}
    return chaine.str();
}
