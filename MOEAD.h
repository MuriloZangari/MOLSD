/*
 * MOEAD.h
 *
 *  Created on: 27/02/2015
 *      
 */

#ifndef MOEAD_H_
#define MOEAD_H_
#include <algorithm>
#include <math.h>
#include <limits>
#include "Global.h"
#include "NTools.h"
#include "Problem.h"
#include "Lambda.h"
#include "Distance.h"
#include "Subprob.h"
#include "NMODEL.h"
#include "Group.h"
#include "ArchivePareto.h"
#include "NEH.h"


using namespace std;


class CMOEAD{


private:


public:
	int numbSubProb; //number of subproblems
    int numbGroup; //number of probabilistic model
    int neighborhoodSizeUpdate; //number of neighbors for update
    int size_selection;
    int maxNumbOfUpdate; //maximum number that one solution can update a set of neighbors
    int numbSamplingSol; //flexible number, each group can have a specific number of sampling solution
    int update; //type of the mechanism of update
    int ScalarFunction; // 0 means weighted sum, and 1 means theebicheffi
    int model;
    float mutation; // parameter to set if use mutation or not for the GA and simple UMDA

    AbstractProblem* problemTest;
    Subprob * solution;		//declaration of a set of subproblems
    Subprob * auxSolution;
    Lambda * lambdaVector;
    ArchivePareto * externalPareto;
    Group * group;
	int * invertedBB;
	int * sol1;
	int * sol2;
	int * compositionn;

	int ** previousFitnessVale;

	int  m_restart_countdown;
	int * m_count_down_for_restart;
	int * m_max_restars;

	int m_MOEAD_restar_shake_power;
	int m_VNS_restart_shake_power;

	bool boltz;
    int t_neighborhood;
    int typeSampling;
    char type[2];

    int theta;

public:
    void InitializeReferencePoint(); //iniatilize the reference point according to if the problem is maximization or minization
    void InitializeWeightVectors(int obj, int numSubProb );//read the file of weight vector according to the number of subproblems and number of objetives
    double EuclideanDist(Lambda lambdaAtual, Lambda lambda2, int obj); //function to calculate the distance between the weight vectors
    void shellSort(Distance* dist, int solutionNumber); // to order the neighboors of each weight vector
    void InitializeNeighorhood(int obj, int numSubProb);  //set the neighborhood of each subproblem, according to the order of the most close neighbors
    void InitializeGroups(Group* group, int model); //initialize the groups, also the model (reproduction mechanism) is set
    void InitializeSolutions(); //initialize the solutions according to the problem specific

    void ShowSolutions(); //show the current set of solutions
    void PrintObjSpaceSol(int i); //show the image of the solutions
    void UpdateReferencePoint(Subprob solution); //update the reference point, according to each new solution sampled (generated) and if the problem is maximization or minimization
    void UpdateNeighboringSolution(Subprob newSol, int index); //conventional update
    void UpdateNeighboringSolutionGR1(Subprob newSol, int index); //Global Replacement (GR), random choose
    void UpdateSelfish(Subprob newsol, int id);
    void UpdateNeighborhoodSize();
//    void UpdateNeighboringSolutionGR2(Lambda* lambda, Subprob newSol, Subprob* solution); //Global Replacement (GR), most appropriate subproblem, than update its neighbors
//    void UpdateNeighboringSolutionDD(int centroid, Lambda* lambda, Subprob newSol, Subprob* solution); //DD, the new solution is associated to the most appropriate subregion, than update the neighborhood of subproblem associate to this subregion
    void ShowVar(int i);
    void Run(); //run the MOEAD with the configuration of parameters
    void ShowExecutionTime(double tmili, int i);
    void ShakeSolution(int index);
    void ShowExtraInformation(int * sh);

       CMOEAD(int n_sub, int neighborSizeSel, int modell, AbstractProblem* problemTestt, int typeNeighborhood, float mut, int ts);
       virtual ~CMOEAD();


};


#endif /* MOEAD_H_ */

