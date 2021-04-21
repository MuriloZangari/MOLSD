/*
 * Problem.h
 *
 *  Created on: 27/02/2015
 *      Author: hydra
 */

#ifndef PROBLEM_H_
#define PROBLEM_H_

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
#include "Subprob.h"
#include "Global.h"
#include "NTools.h"
#include "Lambda.h"
#include "Distance.h"

using namespace std;

class AbstractProblem{
public:
	int numGreedRepair;
	int optimum;
	virtual void Initialize()=0; //initialization of the problem (read the files)
	virtual void LoadInstance()=0; //initialize the problem matrix with the data of the initialization
	virtual void heuristicValue(Lambda* lambda, int numsubProb)=0; //the heuristic value for greedy repair and other methods
	virtual void greedyRepair(Subprob &solution, int indexSub, int &numGR)=0; //greedy repair for constrained binary problems
	virtual void Show()=0; //show the data structure of the problem matrix

	virtual void solutionFitness(Subprob &solution)=0;
	virtual int * solutionFitness(int * gen)=0;
//	virtual void InvSolutionFitness(Subprob &solution)=0;
//	virtual int * InvSolutionFitness(int * gen)=0;

    /*
     * Partial evaluation method for LR constructive heuristic method. The call to this method is expected in the index function.
     */
//    virtual int PartialEvaluation(int * genes, int size, int * new_job_times)=0;
//    virtual int PartialEvaluation(int * genes, int size, int new_job)=0;
//    virtual void IdleTimesAtPosition(int * sequence, int size, int job_i, int * idleTimes)=0;

    /*
     * Fitness for heuristics single objective
     */
    virtual int EvalCmax(int * genes, int size)=0;
    virtual int EvalTFT(int * genes, int size)=0;
    virtual int EvalTWT(int * genes, int size)=0;

	AbstractProblem();
	virtual ~AbstractProblem();

	/*
	 * The processing times matrix.
	 */
	int **m_processingtimes;

	int * m_DueDates;

	int * m_weights;
};

//class mUBQP: public AbstractProblem{
//
//public:
//	int seed;
//	int GainByObjective[objMax];
//	mUBQP(int s);
//	~mUBQP();
//
//public:
//	void Initialize();
//	void LoadInstance ();
//	void heuristicValue(Lambda* lambda, int numsubProb);
//	void greedyRepair(Subprob &solution, int indexSub, int &numGR);
//	void Show();
//	void solutionFitness(Subprob &solution);
//	void FitnessGain(Subprob &solution, Lambda &lambda, int pos);
//	void evalGreed(Subprob &solution, Lambda &lambda);
//    /*
//     * Partial evaluation method for LR constructive heuristic method. The call to this method is expected in the index function.
//     */
//    int PartialEvaluation(int * genes, int size, int * new_job_times){return 0;};
//    int PartialEvaluation(int * genes, int size, int new_job){return 0;};
//    void IdleTimesAtPosition(int * sequence, int size, int job_i, int * idleTimes){};
//
//    void InvSolutionFitness(Subprob &solution){};
//    int * solutionFitness(int * gen){return 0;};
//    int * InvSolutionFitness(int * gen){return 0;};
//};
//
//class Trap5: public AbstractProblem{
//
//public:
//
//	Trap5();
//	~Trap5();
//
//public:
//	void Initialize();
//	void LoadInstance ();
//	void heuristicValue(Lambda* lambda, int numsubProb);
//	void greedyRepair(Subprob &solution, int indexSub, int &numGR);
//	void Show();
//	void solutionFitness(Subprob &solution);
//	void FitnessGain(Subprob &solution, Lambda &lambda, int pos){};
//	void evalGreed(Subprob &solution, Lambda &lambda){};
//
//    /*
//     * Partial evaluation method for LR constructive heuristic method. The call to this method is expected in the index function.
//     */
//    int PartialEvaluation(int * genes, int size, int * new_job_times){return 0;};
//    int PartialEvaluation(int * genes, int size, int new_job){return 0;};
//    void IdleTimesAtPosition(int * sequence, int size, int job_i, int * idleTimes){};
//    void InvSolutionFitness(Subprob &solution){};
//    int * solutionFitness(int * gen){return 0;};
//    int * InvSolutionFitness(int * gen){return 0;};
//};
//
//class MOKP: public AbstractProblem{
//
//private:
//	float capacity[objMax];
//	int profit[objMax][varMax];
//	int weight[objMax][varMax];
//	double heuristic[varMax];
//	int sort[subProbMax][varMax];
//
//public:
//
//	MOKP();
//	~MOKP();
//
//	void setCapacity(float c, int i);
//	float getCapacity(int i);
//	void setProfit(int p, int i, int j);
//	int getProfit(int i,int j);
//	void setWeight(int w, int i, int j);
//	int getWeight(int i,int j);
//	void setHeuristic(double h, int i);
//	double getHeuristic(int i);
//	void setSort(int s, int sub,int v);
//	int getSort(int sub, int v);
//
//
//	void heuristicValue(Lambda* lambda, int numsubProb);
//	bool IsFeasible(Subprob solution);
//	void greedyRepair(Subprob &solution, int indexSub,int &numGR);
//	void solutionFitness(Subprob &solution);
//	void Initialize();
//	void LoadInstance();
//	void Show();
//	void FitnessGain(Subprob &solution, Lambda &lambda, int pos){};
//	void evalGreed(Subprob &solution, Lambda &lambda){};
//
//    /*
//     * Partial evaluation method for LR constructive heuristic method. The call to this method is expected in the index function.
//     */
//    int PartialEvaluation(int * genes, int size, int * new_job_times){return 0;};
//    int PartialEvaluation(int * genes, int size, int new_job){return 0;};
//    void IdleTimesAtPosition(int * sequence, int size, int job_i, int * idleTimes){};
//    void InvSolutionFitness(Subprob &solution){};
//    int * solutionFitness(int * gen){return 0;};
//    int * InvSolutionFitness(int * gen){return 0;};
//};
//
//class City: public AbstractProblem{
//
//private:
//	int x[varMax][objMax];
//	int y[varMax][objMax];
//
//public:
//	City();   // This is the constructor declaration
//	~City();  // This is the destructor: declaration
//
//	void setX( int xis, int i, int j );
//	int getX( int i, int j );
//	void setY( int yis, int i, int j );
//	int getY( int i,int j );
//
//	void solutionFitness(Subprob &solution);
//	void Initialize();
//	void LoadInstance();
//	void heuristicValue(Lambda* lambda, int numsubProb);
//	void greedyRepair(Subprob &solution, int indexSub,int &numGR){}; //greedy repair for constrained binary problems
//	void Show();
//	void FitnessGain(Subprob &solution, Lambda &lambda, int pos){};
//	void evalGreed(Subprob &solution, Lambda &lambda){};
//    /*
//     * Partial evaluation method for LR constructive heuristic method. The call to this method is expected in the index function.
//     */
//    int PartialEvaluation(int * genes, int size, int * new_job_times){return 0;};
//    int PartialEvaluation(int * genes, int size, int new_job){return 0;};
//    void IdleTimesAtPosition(int * sequence, int size, int job_i, int * idleTimes){};
//    void InvSolutionFitness(Subprob &solution){};
//    int * solutionFitness(int * gen){return 0;};
//    int * InvSolutionFitness(int * gen){return 0;};
//};
//
//class LOP: public AbstractProblem{
//
//public:
//	/*
//	 * The size of the problem.
//	 */
//	int m_problemsize;
//
//	LOP();
//	~LOP();
//
//public:
//	void Initialize();
//	void LoadInstance ();
//	void heuristicValue(Lambda* lambda, int numsubProb){};
//	void greedyRepair(Subprob &solution, int indexSub, int &numGR){};
//	void Show();
//	void solutionFitness(Subprob &solution);
//	void FitnessGain(Subprob &solution, Lambda &lambda, int pos){};
//	void evalGreed(Subprob &solution, Lambda &lambda){};
//    /*
//     * Partial evaluation method for LR constructive heuristic method. The call to this method is expected in the index function.
//     */
//    int PartialEvaluation(int * genes, int size, int * new_job_times){return 0;};
//    int PartialEvaluation(int * genes, int size, int new_job){return 0;};
//    void IdleTimesAtPosition(int * sequence, int size, int job_i, int * idleTimes){};
//    void InvSolutionFitness(Subprob &solution){};
//    int * solutionFitness(int * gen){return 0;};
//    int * InvSolutionFitness(int * gen){return 0;};
//};

class PFSP: public AbstractProblem{
public:

    /*
     * Auxiliary vector for inversion.
     */
    int * m_aux;
	/*
	 * The number of jobs of the problem.
	 */
	int m_jobs;

	/*
	 * The number of machines of the problem.
	 */
	int m_machines;

    /*
     * The time table for the processing times.
     */
    int * m_timeTable;

	// The constructor. It initializes a flowshop scheduling problem from a file.
	PFSP();

    // The destructor.
    virtual ~PFSP();

	void Initialize();
	void LoadInstance ();
	void heuristicValue(Lambda* lambda, int numsubProb){};
	void greedyRepair(Subprob &solution, int indexSub, int &numGR){};
	void Show();
	void solutionFitness(Subprob &solution);


	int GetProblemSize();
    /*
     * Partial evaluation method for LR constructive heuristic method. The call to this method is expected in the index function.
     */
//    int PartialEvaluation(int * genes, int size, int * new_job_times);
//    int PartialEvaluation(int * genes, int size, int new_job);
//    void IdleTimesAtPosition(int * sequence, int size, int job_i, int * idleTimes);

//    void InvSolutionFitness(Subprob &solution);
    int * solutionFitness(int * gen);
//    int * InvSolutionFitness(int * gen);

    int EvalCmax(int * genes, int size){return 0;};
    int EvalTFT(int * genes, int size){return 0;};
    int EvalTWT(int * genes, int size){return 0;};

private:

};

class PFSP_tard: public AbstractProblem{
public:

    /*
     * Auxiliary vector for inversion.
     */
    int * m_aux;
	/*
	 * The number of jobs of the problem.
	 */
	int m_jobs;

	/*
	 * The number of machines of the problem.
	 */
	int m_machines;

    /*
     * The time table for the processing times.
     */
    int * m_timeTable;


	// The constructor. It initializes a flowshop scheduling problem from a file.
	PFSP_tard();

    // The destructor.
    virtual ~PFSP_tard();

	void Initialize();
	void LoadInstance ();
	void heuristicValue(Lambda* lambda, int numsubProb){};
	void greedyRepair(Subprob &solution, int indexSub, int &numGR){};
	void Show();
	void solutionFitness(Subprob &solution);


	int GetProblemSize();
    /*
     * Partial evaluation method for LR constructive heuristic method. The call to this method is expected in the index function.
     */
//    int PartialEvaluation(int * genes, int size, int * new_job_times);
//    int PartialEvaluation(int * genes, int size, int new_job);
//    void IdleTimesAtPosition(int * sequence, int size, int job_i, int * idleTimes);

//    void InvSolutionFitness(Subprob &solution);
    int * solutionFitness(int * gen);
//    int * InvSolutionFitness(int * gen);

    int EvalCmax(int * genes, int size){return 0;};
    int EvalTFT(int * genes, int size){return 0;};
    int EvalTWT(int * genes, int size){return 0;};
private:
};

class PFSP_3ob: public AbstractProblem{
public:

    /*
     * Auxiliary vector for inversion.
     */
    int * m_aux;
	/*
	 * The number of jobs of the problem.
	 */
	int m_jobs;

	/*
	 * The number of machines of the problem.
	 */
	int m_machines;

    /*
     * The time table for the processing times.
     */
    int * m_timeTable;


	// The constructor. It initializes a flowshop scheduling problem from a file.
    PFSP_3ob();

    // The destructor.
    virtual ~PFSP_3ob();

	void Initialize();
	void LoadInstance ();
	void heuristicValue(Lambda* lambda, int numsubProb){};
	void greedyRepair(Subprob &solution, int indexSub, int &numGR){};
	void Show();
	void solutionFitness(Subprob &solution);


	int GetProblemSize();
    /*
     * Partial evaluation method for LR constructive heuristic method. The call to this method is expected in the index function.
     */
//    int PartialEvaluation(int * genes, int size, int * new_job_times);
//    int PartialEvaluation(int * genes, int size, int new_job);
//    void IdleTimesAtPosition(int * sequence, int size, int job_i, int * idleTimes);

//    void InvSolutionFitness(Subprob &solution);
    int * solutionFitness(int * gen);
//    int * InvSolutionFitness(int * gen);

    int EvalCmax(int * genes, int size){return 0;};
    int EvalTFT(int * genes, int size){return 0;};
    int EvalTWT(int * genes, int size){return 0;};
private:

};

///////////////////////////////////sdst Cmax and TWT //////////////////////////////////////////////////////////////////////////////

class PFSP_SDST: public AbstractProblem{
public:

    /*
     * Auxiliary vector for inversion.
     */
    int * m_aux;
	/*
	 * The number of jobs of the problem.
	 */
	int m_jobs;

	/*
	 * The number of machines of the problem.
	 */
	int m_machines;


	/*
     * The time table for the processing times.
     */
    int * m_timeTable;

    /*
     * The sequence dependent times for each machine and jobs combination
     */
    int ***m_SDST;


	// The constructor. It initializes a flowshop scheduling problem from a file.
	PFSP_SDST();

    // The destructor.
    virtual ~PFSP_SDST();

	void Initialize();
	void LoadInstance ();
	void heuristicValue(Lambda* lambda, int numsubProb){};
	void greedyRepair(Subprob &solution, int indexSub, int &numGR){};
	void Show();
	void solutionFitness(Subprob &solution);

	int GetProblemSize();
    /*
     * Partial evaluation method for LR constructive heuristic method. The call to this method is expected in the index function.
     */
//    void IdleTimesAtPosition(int * sequence, int size, int job_i, int * idleTimes){};
//    int PartialEvaluation(int * genes, int size, int * new_job_times){return 0;};
//    int PartialEvaluation(int * genes, int size, int new_job){return 0;};


//    void InvSolutionFitness(Subprob &solution);
    int * solutionFitness(int * gen);
//    int * InvSolutionFitness(int * gen);

    int EvalCmax(int * genes, int size);
    int EvalTFT(int * genes, int size){return 0;};
    int EvalTWT(int * genes, int size);
private:
};

//////////////////////////////////////////////////////SDST Cmax and TFT ///////////////////////////////////////////////////////////

class SDST_Cmax_TFT: public AbstractProblem{
public:

    /*
     * Auxiliary vector for inversion.
     */
    int * m_aux;
	/*
	 * The number of jobs of the problem.
	 */
	int m_jobs;

	/*
	 * The number of machines of the problem.
	 */
	int m_machines;


	/*
     * The time table for the processing times.
     */
    int * m_timeTable;

    /*
     * The sequence dependent times for each machine and jobs combination
     */
    int ***m_SDST;


	// The constructor. It initializes a flowshop scheduling problem from a file.
    SDST_Cmax_TFT();

    // The destructor.
    virtual ~SDST_Cmax_TFT();

	void Initialize();
	void LoadInstance ();
	void heuristicValue(Lambda* lambda, int numsubProb){};
	void greedyRepair(Subprob &solution, int indexSub, int &numGR){};
	void Show();
	void solutionFitness(Subprob &solution);

	int GetProblemSize();
    /*
     * Partial evaluation method for LR constructive heuristic method. The call to this method is expected in the index function.
     */
//    void IdleTimesAtPosition(int * sequence, int size, int job_i, int * idleTimes){};
//    int PartialEvaluation(int * genes, int size, int * new_job_times){return 0;};
//    int PartialEvaluation(int * genes, int size, int new_job){return 0;};


//    void InvSolutionFitness(Subprob &solution);
    int * solutionFitness(int * gen);
//    int * InvSolutionFitness(int * gen);

    int EvalCmax(int * genes, int size);
    int EvalTFT(int * genes, int size);
    int EvalTWT(int * genes, int size){return 0;};
private:
};
#endif /* PROBLEM_H_ */
