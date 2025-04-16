#ifndef DATA_H_
#define DATA_H_

#include<string>

using namespace std;
class Data
{
private :
    unsigned int Njobs ;
    unsigned int Nresources ;
    unsigned int horizon ;
    unsigned int duedate ;
    unsigned int tardcost ;
    unsigned int MPMtime ;

    int ** E ;
    int *nbE ;

    int *P ;
    int *B ;
    int **b ;

public :
    Data(const string namefile);

    int getNjobs() const ;
    int getNresources() const ;
    int gethorizon() const ;
    int getduedate() const ;
    int gettardcost() const ;
    int getMPMtime() const ;
    int** getE() const ;
    int* getLineE(unsigned int i) const ;
    int getValueE(unsigned int i,unsigned  int j) const ;
    int* getnbE() const ;
    int getnbEvelue(unsigned int i) const ;
    int* getP() const ;
    int getValueP(unsigned int i) const ;
    int* getB() const ;
    int getValueB(unsigned int i) const ;
    int** getb() const ;
    int* getLineb(unsigned int i) const ;
    int* getColumnb(unsigned int i) const ;
    int getValueb(unsigned int i,unsigned  int j) const ;

    bool arc(int i, int j);
    int getsumP() const;
    int getmaxP() const;
    int getminP() const;

    string toString();

    ~Data()
    {
        for(unsigned int i=0 ; i<Njobs ; ++i){
			delete[] E[i];
			delete[] b[i];
		}
		delete[] E;
		delete[] b;

		delete[] nbE;
		delete[] P;
		delete[] B;
    }

};

#endif // DATA_H_INCLUDED
