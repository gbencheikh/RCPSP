//============================================================================
// Name        : RCPSP.cpp
// Author      : GBENCHEIKH and MSADIK
// Version     : 18/01/2016
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include"Data.h"
#include"Resolution.h"
#include"Heuristique.h"
#include"CutsGeneration.h"
#include<string>
#include<sstream>
#include <limits>

using namespace std;

int main() {


	ofstream f("resultats.dat",ios::out);
	//f << "nb;Param;Instance;solution;time"<<endl;

	//for(int nbtaches = 30; nbtaches <= 120 ; nbtaches += 30 ){
	int nbtaches = 30;
		for(int param = 1 ; param <= 48 ; param++){
			for(int instance = 1 ; instance <= 10 ; instance++){
				f << nbtaches << ";" << param << ";" << instance << ";";
				stringstream information;
				string filename;

				information << "Instances/";
				information << nbtaches;
				information << "/";
				information << "j";
				information << nbtaches;
				information << ".sm/j";
				information << nbtaches;
				information << param;
				information <<  "_";
				information << instance;
				information << ".sm";

				filename = information.str();

				/*
				CutsGeneration cuts(filename);
				cuts.critiqueJobs();
				cuts.incompatibleJobs();
				f << cuts.resolvingModel();
				*/

				/*
				Resolution solOpt(filename);
				f << solOpt.resolvingTimeIndexed() << endl;
				*/

				/*
				Heuristique sol(filename);
				sol.creatList();
				sol.parallele();
				*/
			}
		}

	cout << "fin" << endl;

	//string filename = "Instances/30/j30.sm/j301_1.sm";
	//Data data(filename);

	//Heuristique sol(filename);
	//sol.parallele();
	//sol.printSolution();

	//Resolution sol(filename);
	//cout << sol.resolvingTimeIndexed() << endl;
	//cout << sol.resolvingFlot();

	//CutsGeneration cuts(filename);
	//cuts.critiqueJobs();
	//cout << cuts.resolvingModel();

    return 0;
}

