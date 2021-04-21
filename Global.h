/*
 * Global.h
 *
 *  Created on: 27/02/2015
 *      Author: hydra
 */



#ifndef GLOBAL_H_
#define GLOBAL_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <ctime>
#include <string.h>
#include <sys/time.h>
#include <vector>


using namespace std;

extern	int numbObj;
extern	int numbVar;
extern  int scalarFunction;
extern	string idxInstance;
extern	int numbGenerations;
extern	int problem;
extern  long int numbEval;
extern int cycle;
extern int machine;

extern float probGenetic;
extern int n_sampling;
extern int n_size_update;
extern int initial_pop;
extern int benchmark;
extern int shake;
extern float probX;

extern float LR;

extern string pho, test;

extern	int run_id;

	const long int intMin = -999999999;
	const long int intMax = 999999999;
	const int objMax = 10;
	const int varMax = 3000;
	const int subProbMax=320;
	const int neighborMax = subProbMax;
	const int numSamplMax = subProbMax;
	const int ParetoSolutionsMax = 20000;
	extern int * ReferencePoint;
	extern int * worstReferencePoint;

	extern 	vector< vector <vector <int> > > matrixOfItems; //matrix for the mUBQP
	extern vector< vector <vector <float> > > matrixOfEdges; //matrix of edges for MTSP.  NOTE: They are global because to pass these large structures as parameter for the functions increase a lot the execution time, thus, by now each problem that has a different structure of data has its declaration
	extern 	vector<int> vectorOfNeighbors;

	extern long int moead_maxEval;
	extern long int MAX_EVAL;

#define MAX(A,B) ( (A > B) ? A : B)
#define MIN(A,B) ( (A < B) ? A : B)

#endif /* GLOBAL_H_ */
