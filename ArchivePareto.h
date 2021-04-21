/*
 * ArchivePareto.h
 *
 *  Created on: 06/03/2015
 *      Author: hydra
 */

#ifndef ARCHIVEPARETO_H_
#define ARCHIVEPARETO_H_

#include"Global.h"
#include"Subprob.h"


using namespace std;

class ArchivePareto{
public:
	bool add;
    int EP[ParetoSolutionsMax][objMax]; //structure to maintain the f(x) of the solutions that are non dominated so far
    int auxEP[ParetoSolutionsMax][objMax]; //the auxiliary structure. NOTE: this way is faster than to use another structure as vector or List
    int sizeEP,sizeAux; //parameters for update the EP with the non dominated so far
    void UpdateExternalPareto(Subprob solution); //update the EP according to if the problem is maximization or minimization
    void ShowExternalPareto(int gen); //print in a file the EP solutions of the run

    ArchivePareto();
    ~ArchivePareto();
};

#endif /* ARCHIVEPARETO_H_ */
