#ifndef HEURISTIQUE_H_
#define HEURISTIQUE_H_

#include<iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

using namespace std;

class Heuristique{

private:
    int Njobs;
	int Nresources;

	int **E;
	int *nbE;
	int *P;
	int *B;
	int **b;

	int *liste;
	int *S;
	int Z;

	bool arc(int i, int j);
	int getmaxP() const ;
	int getminP() const ;
	int getsumP() const ;

public:
	Heuristique(string);
	void creatList();
	void serie();
	void parallele();
	bool canBeThere(int i, int t, int **r);

	int getZ(){return Z;}
	int getS(int i){return S[i];}

	void printList();
	void printResource();
	void printSolution();

	bool resource();

};

#endif // HEURISTIQUE_H_
