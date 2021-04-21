/*
 * MODEL.h
 *
 *  Created on: 27/02/2015
 *      Author: hydra
 */

#ifndef MODEL_H_
#define MODEL_H_

#include "Global.h"
#include "Subprob.h"
#include "Distance.h"
#include "Lambda.h"
#include "NTools.h"


using namespace std;

class AbstractModel{
	public: //the parameters that are needed for learning and sampling, that are passed at the initialization of the model
	int numbSubProb; //number of subproblems
	int FreqMatrix[varMax][varMax];
	double* probVec;
	public:
		AbstractModel(); //initialization of the parameters
		virtual ~AbstractModel();
		virtual void ComputeFrequencies(){};
		virtual void Learning(Lambda lambda, int numSamp)=0; //abstract method learning
				virtual void Sampling(Subprob* solution, int s, Subprob &newSol )=0; //abstract method sampling that returns a new subproblem
		virtual void Sampling(Subprob &newSol )=0; //abstract method sampling that returns a new subproblem
		virtual void InitializeModel(){}; //initialization of the model
		virtual void ShowModel(){}; //show some learned information of the model
		virtual void ShowModel(int n){};
		virtual void PrintFreqMatrix(int i){};
		void PrintSelectPop(Lambda lambda, Subprob* solution);
		//void ProbVectorScFun(Lambda lambda, Subprob* solution, double* probVec);
};

//class PBILModel:public AbstractModel{
//public:
//	double prob[varMax]; //do PBIL
//	float alpha;
//	int mutation;
//	int sigma;
//	PBILModel(int nsp, int mut, float lr);
//	~PBILModel();
//	 void ComputeFrequencies(){};
//	 void Learning(Lambda lambda, Subprob* solution, int numSamp);
//	 void Learning(Lambda lambda, Subprob* solution, int numSamp, float theta){};
//	 void Sampling(Subprob* solution, int s, Subprob &newSol);
//	 void Sampling(Subprob &newSol ){};
//
//	 void InitializeModel();
//	 void ShowModel();
//
//	 //specific methods from PBIL
//	 void InitializeProb();
//	 void PrintFreqMatrix(int i){};
//};
//
//class PBILOrderInvers:public AbstractModel{
//public:
//	double prob[varMax]; //do PBIL
//	double alpha;
//	int mutation;
//	int sigma;
//	PBILOrderInvers(int nsp, int mut, float lr);
//	~PBILOrderInvers();
//	void ComputeFrequencies(){};
//	 void Learning(Lambda lambda, Subprob* solution, int numSamp);
//	 void Learning(Lambda lambda, Subprob* solution, int numSamp, float theta){};
//	 void Sampling(Subprob* solution, int s, Subprob &newSol);
//	 void Sampling(Subprob &newSol ){};
//	 void InitializeModel();
//	 void ShowModel();
//
//	 //specific methods from PBIL
//	 void InitializeProb();
//	 void PrintFreqMatrix(int i){};
//};
//
//class PBILOrder:public AbstractModel{
//public:
//	double prob[varMax]; //do PBIL
//	double alpha;
//	int mutation;
//	int sigma;
//	PBILOrder(int nsp, int mut, float lr);
//	~PBILOrder();
//	void ComputeFrequencies(){};
//	 void Learning(Lambda lambda, Subprob* solution, int numSamp);
//	 void Learning(Lambda lambda, Subprob* solution, int numSamp, float theta){};
//	 void Sampling(Subprob* solution, int s, Subprob &newSol);
//	 void Sampling(Subprob &newSol ){};
//	 void InitializeModel();
//	 void ShowModel();
//
//	 //specific methods from PBIL
//	 void InitializeProb();
//	 void PrintFreqMatrix(int i){};
//};
//
//class IUMDA:public AbstractModel{
//public:
//	double prob[varMax]; //do PBIL
//	float factoR;
//	int mutation;
//	bool bolztman;
//	IUMDA(int nsp, int mut, float factor, bool boltz);
//	~IUMDA();
//	void ComputeFrequencies(){};
//	 void Learning(Lambda lambda, Subprob* solution, int numSamp);
//	 void Learning(Lambda lambda, Subprob* solution, int numSamp, float theta){};
//	 void Sampling(Subprob* solution, int s, Subprob &newSol);
//	 void Sampling(Subprob &newSol ){};
//	 void InitializeModel();
//	 void ShowModel();
//
//	 //specific methods from PBIL
//	 void InitializeProb();
//	 void PrintFreqMatrix(int i){};
//};


class GA:public AbstractModel{
public:
	int p1;
	int p2;
	double sigma;
	float mutation;
	int index;
	GA(int nsp, float mut);
	~GA();
	void ComputeFrequencies(){};
	void Learning(Lambda lambda, int numSamp);
	void Sampling(Subprob* solution, int s, Subprob &newSol);
	void Sampling(Subprob &newSol ){};
	 void InitializeModel();
	 void ShowModel(int n);
	 void PrintFreqMatrix(int i){};
};

//////////////////////LOCAL SEARCH/////////////////////////////////

class LS:public AbstractModel{
public:
	int index;

	float probMut;

	LS(int ins);
	~LS();
	void ComputeFrequencies(){};
	void Learning(Lambda lambda, int numSamp);

	void Sampling(Subprob* solution, int s, Subprob &newSol);
	void Sampling(Subprob &newSol){};
	 void InitializeModel(){};
	 void ShowModel(int n);

	 void PrintFreqMatrix(int i){};
};

//////////////////////LOCAL SEARCH/////////////////////////////////

class GLS:public AbstractModel{
public:
	int p1;
	int p2;
	int index;
	float sigma;
	float mutation;
	float probX;
	GLS(float probC, float mut, float pX);
	~GLS();
	void ComputeFrequencies(){};
	void Learning(Lambda lambda, int numSamp);
	void Sampling(Subprob* solution, int s, Subprob &newSol);
	void Sampling(Subprob &newSol ){};
	 void InitializeModel();
	 void ShowModel(int n);
	 void PrintFreqMatrix(int i){};
};

////////////////////// TREE  ////////////////////////////////////

//class TREE: public AbstractModel, public IntTreeModel{
//public:
//		bool bolztman;
//		int mutation;
//		double Prior;
//		int sigma;
//        int typetree;         // Three possible types of models
//                              // typetree = 0, Univariate model (as used by UMDA)
//                              // typetree = 1, a chain-shaped model in which each variable
//                              // depends on the previous in the order given by the representation
//                              //  x1 depends on x0, x2 depends on x1, and xi depends on xi-1
//                              // typetree = 2, tree with disconnected components (forest)
//     TREE(int nsp, int ttree, unsigned int* cards, int mut ,double prior, bool boltz);
//	~TREE();
//	void ComputeFrequencies();
//    void Learning(Lambda lambda, Subprob* solution, int numSamp);
//    void Learning(Lambda lambda, Subprob* solution, int numSamp, float theta){};
//	void Sampling(Subprob* solution, int s, Subprob &newSol);
//	void Sampling(Subprob &newSol ){};
//	void InitializeModel();
//	void ShowModel();
//	void PrintFreqMatrix(int i);
//};





////////////////////// MALLOWS MODEL ////////////////////////////////////

//class MALLOWS: public AbstractModel, public CMallowsModel{
//public:
//
//
//
//
//        float* float_probVec;            // The Mallows code uses float precision
//                                         // so we will have to transform from double to float during learning
//
//
//	 MALLOWS(int nsp,  char* metric_type);
//	~MALLOWS();
//
//    void Learning(Lambda lambda, Subprob* solution, int numSamp){};
//    void Learning(Lambda lambda, Subprob* solution, int numSamp, float theta);
//	void Sampling(Subprob* solution, int s, Subprob &newSol);
//	void Sampling(Subprob &newSol);
//	void InitializeModel();
//	void ShowModel();
//	void PrintFreqMatrix(int i);
//};




#endif /* MODEL_H_ */
