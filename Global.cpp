#include "Global.h"

	vector< vector <vector <int> > > matrixOfItems; //vector for the mUBQP or vector of edges for the MTSP
	vector< vector <vector <float> > > matrixOfEdges; //vector for the edges for MTSP
	vector<int> vectorOfNeighbors;

	int numbObj;
	int numbVar;
	string idxInstance;
	int numbGenerations;
	int problem;
	int scalarFunction;
	long int numbEval;
	int machine;
	int run_id;
	int cycle;
	float probGenetic;
	int n_sampling;
	int n_size_update;
	int initial_pop;
	float LR;
	int benchmark;
	int shake;
	float probX;

	long int moead_maxEval;
	long int MAX_EVAL;

	string pho, test;
	int * ReferencePoint;
	int * worstReferencePoint;





