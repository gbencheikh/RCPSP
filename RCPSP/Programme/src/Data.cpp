#include"Data.h"
#include<iostream>
#include<fstream>
#include<string>
#include<sstream>

using namespace std;

Data :: Data(string namefile)
{
	Njobs = 0;
	Nresources = 0;
	horizon = 0;
	Nresources = 0;
	duedate = 0;
	tardcost = 0;
	MPMtime = 0;

	ifstream file(namefile.c_str(),ios::in);
    if(file)
    {
        string c;
        file >> c;
        while (c != "):")
            file >> c ;
        int N;
        file >> N; // read number of jobs
        Njobs = N;
        file >> c >> c ; // read character
        file >> horizon ; // read horizon
        file >> c;
        while ( c != ":")
            file >> c;
        file >> Nresources; // read number of resources

        while (c != "MPM-Time")
            file >> c;
        file >> c >> c >> c;
        file >> duedate >> tardcost >> MPMtime; // Read due date and MPM-Time

        /***** allocation of matrix E, nbE , P , B and b *****/

        E = new int *[Njobs];
        P = new int [Njobs];
        nbE = new int [Njobs];
        B = new int [Nresources];
        b = new int *[Njobs];
        for(unsigned int i = 0 ; i < Njobs ; i++)
            b[i] = new int [Nresources];

        /**********     Continue reading  ***********/
        // Read precedence relations
        while(c != "successors")
             file >> c;
        for(unsigned int i = 0 ; i < Njobs ; i++)
        {
            file >> c >> c ;
            file >> nbE[i]; // read number of successors of job i
            E[i] = new int[nbE[i] ];
            for(int j = 0 ; j < nbE[i] ; j++){
            	int temp;
            	file >> temp;
				E[i][j] = temp-1;
            }

        }
        while (c != "REQUESTS/DURATIONS:")
            file >> c;

        getline(file, c);
        getline(file, c);
        getline(file, c);
        // read duration and consumption
        for(unsigned int i = 0 ; i < Njobs ; i++)
        {
            file >> c >> c;
            file >> P[i];
            for(unsigned int j = 0 ; j < Nresources ; j++)
                file >> b[i][j];
        }
        while( c != "RESOURCEAVAILABILITIES:")
        	file >> c;

        getline(file, c);
        getline(file, c);

        for(unsigned int i = 0 ; i < Nresources ; i++)
            file >> B[i];
    }
    else
        cout << "file cannot be open" << endl;
}
int Data :: getNjobs() const {return Njobs; }
int Data :: getNresources() const {return Nresources; }
int Data :: gethorizon() const {return horizon; }
int Data :: getduedate() const {return duedate ;}
int Data :: gettardcost() const {return tardcost; }
int Data :: getMPMtime() const {return MPMtime; }

int** Data :: getE() const {return E;}
int* Data :: getLineE(unsigned int i) const {return E[i]; }
int Data :: getValueE(unsigned int i,unsigned  int j) const {return E[i][j]; }
int* Data :: getnbE() const {return nbE; }
int Data :: getnbEvelue(unsigned int i) const {return nbE[i]; }

int* Data :: getP() const {return P; }
int Data :: getValueP(unsigned int i) const {return P[i];}
int* Data :: getB() const {return B; }
int Data :: getValueB(unsigned int i) const {return B[i];}
int** Data :: getb() const {return b; }
int* Data :: getLineb(unsigned int i) const {return b[i]; }
int* Data :: getColumnb(unsigned int j) const {
	int *temp;
	temp = new int[Njobs];
	for(unsigned int i = 0 ; i < Njobs ; i++)
		temp[i] = b[i][j];
	return temp;
}
int Data :: getValueb(unsigned int i,unsigned  int j) const {return b[i][j]; }

bool Data :: arc(int i, int j){
	for(int k = 0 ; k < nbE[i] ; k++){
		if(j == E[i][k])
			return true;
	}
	return false;
}
int Data :: getsumP() const{
	int sum = 0;
	for(unsigned int i = 0 ; i < this->Njobs ; i++)
		sum += this->P[i];
	return sum;
}
int Data :: getmaxP() const{
	int max = 0;
	for(unsigned int i = 0 ; i < this->Njobs ; i++)
		if(this->P[i] > this->P[max])
			max = i;
	return this->P[max];
}
int Data :: getminP() const{
	int min = 1;
	for(unsigned int i = 2 ; i < this->Njobs-1 ; i++)
		if(this->P[i] < this->P[min])
			min = i;
	return this->P[min];
}
string Data :: toString()
{
    stringstream information;
    information << "Njobs = " << this->Njobs << "\n";
    information << "Nresources = " << this->Nresources << "\n";
    information << "horizon = " << this->horizon << "\n";
    information << "duedate = " << this->duedate << "\n";
    information << "tard cost = " << this->tardcost << "\n";
    information << "MPM-Time = " << this->MPMtime << "\n";

    information << "Precedence relations : \n";
    information << "job : number of precedences - name of successor jobs \n";
    for(unsigned int i = 0 ; i < this->Njobs ; i++)
    {
        information << i << " : " << nbE[i] << " - ";
        for(int j = 0 ; j < this->nbE[i] ; j++)
            information << E[i][j] << " ";
        information << "\n";
    }
    information << "job : duration - consumption on each resource : \n";
    for(unsigned int i = 0 ; i < Njobs ; i++)
    {
        information << i << " : " << this->P[i] << " - ";
        for(unsigned int j = 0 ; j < this->Nresources ; j++)
            information << b[i][j] << " ";
        information << "\n";
    }
    information << "Resource availabilities :\n";
    for(unsigned int i = 0 ; i < Nresources ; i++)
            information << B[i] << " ";
    information << "\n";
    return information.str();
}
