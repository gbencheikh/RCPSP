#ifndef CUTSGENERATION_H_
#define CUTSGENERATION_H_

#include<iostream>
#include <fstream>
#include <string>
#include <vector>
#include<stack>
#include <sstream>
#include<list>
#include"Resolution.h"
#include <ilcplex/ilocplex.h>
#include<ilconcert/iloexpression.h>

using namespace std;

typedef IloArray<IloFloatVarArray> FloatVarMatrix;

class CutsGeneration{

private:
    int Njobs;
	int Nresources;

	int **E;
	int *nbE;
	int *P;
	int *B;
	int **b;

	list<list<int> > critique;
	list<list<int> > incompatible;

	bool arc(int i, int j);
	int getmaxP() const ;
	int getminP() const ;
	int getsumP() const ;

public:

	CutsGeneration(string);
	void incompatibleJobs();
	void critiqueJobs();

	void model(IloModel& model,BoolVarMatrix& x, IloFloatVarArray& S);
    double callModel(BoolVarMatrix& x, IloFloatVarArray& S);
    string resolvingModel();
    int M();
    void makecuts(const BoolVarMatrix& x, int **a, IloExprArray cuts);

    bool tupleCanBeTogether(int i, int j, int t, int k);
	bool tupleCanBeTogether(int i, int j, int t);

	bool canBeTogether(int i, int j, int k);
	bool canBeTogether(int i, int j);

    void printCritiqueJobs();
    void printIncompatibleJobs();
};

#endif // CUTSGENERATION_H_
