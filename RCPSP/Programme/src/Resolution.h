#ifndef RESOLUTION_H_
#define RESOLUTION_H_

#include<iostream>
#include <ilcplex/ilocplex.h>
#include<ilconcert/iloexpression.h>
#include <fstream>
#include <string>
#include <vector>
#include<list>
#include <sstream>

using namespace std;

typedef IloArray<IloBoolVarArray> BoolVarMatrix;
typedef IloArray<IloFloatVarArray> FloatVarMatrix;

typedef IloArray<IloArray<IloIntVarArray> > IntVarMatrix3;

class Resolution
{
private :
    int Njobs;
	int Nresources;

	int **E;
	int *nbE;
	int *P;
	int *B;
	int **b;

	int *ES;
	int *LS;

	list<list<int> > critique;
	list<list<int> > incompatible;

	bool arc(int i, int j);
	int getmaxP() const ;
	int getminP() const ;
	int getsumP() const ;

	int calculM(int i, int j) ;

public :
	Resolution(string filename);

	/*********** Flood's formulation ************/
	void ModelFlot(IloModel& model,BoolVarMatrix& x, IloFloatVarArray& S, IntVarMatrix3& f);
    void ModelFlotRelaxed(IloModel& model,FloatVarMatrix& x, IloFloatVarArray& S, IntVarMatrix3& f);
	double callModelFlot(BoolVarMatrix& x, IloFloatVarArray& S, IntVarMatrix3& f);
    string resolvingFlot();

    /******* Time indexed's formulation ********/

    void ModelTimeIndexed(IloModel& model,BoolVarMatrix& y, int T);
    void ModelTimeIndexedRelaxed(IloModel& model,FloatVarMatrix& y, int T);
    double callModelTimeIndexed(BoolVarMatrix& y);
    string resolvingTimeIndexed();


    /******** Calculate Time window  ***********/

    void earliestStarting();
    void latestStarting(string);

    /******** Valid inequalities  ***********/

    void incompatibleJobs();
    void critiqueJobs();
    void makecuts(const BoolVarMatrix& x, int **a, IloExprArray cuts);

    bool tupleCanBeTogether(int i, int j, int t, int k);
	bool tupleCanBeTogether(int i, int j, int t);
	bool canBeTogether(int i, int j, int k);
	bool canBeTogether(int i, int j);

    /*******************************************/
    bool precedenceFeasibility(bool **x);

    //void printSolution();
    string printESandLS();

    /*******************************************/
};

#endif // RESOLUTION_H_
