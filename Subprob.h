/*
 * Subprob.h
 *
 *  Created on: 27/02/2015
 *      Author: hydra
 */

#ifndef INDIVIDUAL_H_
#define INDIVIDUAL_H_
#include <algorithm>
#include "Global.h"
#include "NTools.h"
#include "Lambda.h"

using namespace std;

class Subprob{


private:
	int index; //index of the subproblem that corresponds to a weigh vector	
    int vectorSolution[varMax]; //this is the vector of the solution
	int function[objMax]; //the image of the vector solution, according to the problem (fitness evaluation, number of objetives)
	double scalarValue; //the scalar value of the solution, i.e., using a scalarizing function (as Tchebyccheff) to estimate a single scalar value


	//this attribute is only for MTSP
	bool noVisited[varMax]; //is used to know what index (cities) was already visit, and can not be used again during the construction of a solution (permutation)
	int initialization; //parameter to know is the initialization is by problem specific or in a random way

	//this is only for MOKP
	int weigthAcc[objMax];
	float heuristicValue[varMax];
public:

	void InitializeSolution(int problem, Lambda lambda); //initialize the solutions according to the problem specified and if is in random way or problem specific
//	void FitnessEvaluation(); //calculate the f(x) of the solutions according to the problem specified
	double ScalarFitnessValueNormalized(Lambda lambda, int scalar); //calculate the scalarizing function according value according to the aggregation function (weighted sum or Tchebycheff)
	double ScalarFitnessValue(Subprob &solution,Lambda lambda, int scalar);

	//gets and sets of index, vectorSolution[], function, scalarValue
	int getIndex(void);
	void setIndex(int ind);
	int getFunction(int i);
	void setFunction(int fc, int i);
	int getVectorSolution(int i);
	void setVectorSolution(int sol, int i);
	double getScalarValue (void);
	void setScalarValue(double sc);

	bool getNoVisited(int i);
	void setNoVisited(bool nv, int i);

	int getWeigthAcc(int i);
	void setWeightAcc(int ww, int i);
	void setHeuristic(float h, int i);
	float getHeuristic(int i);

	int* getPointerToSolution(); 
	int* getPointerToFitness();
	void PrintSolution();
	void SetSoluion(int * vector);
	void SetFitness(int * fitness);
	void SetScalar(double sc);


public:
       
	Subprob();
	~Subprob();

};


#endif /* INDIVIDUAL_H_ */
