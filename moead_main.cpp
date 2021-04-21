#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <algorithm>
#include <ctime>
#include <string.h>
#include <sys/time.h>
#include <vector>
#include "Global.h"
#include "MOEAD.h"
#include "Problem.h"
#include "Subprob.h"


using namespace std;

int main(int argc, char* argv[]){

	//imput parameters, type of model, neighborhood size, type of neighborhood, type of sampling, etc

	int model, neighborhoodSizeSelection,typeNeighbor, seed_id, typeSampling; //type of the model of learning/reproduction //0 means PBIL, 1 means GA, 2 = UMDA, 3= CHAINS model, 4= TREE model
	float mutation;
	if(argc!=21){
		cout << "Specify the Scalar " << endl;
		cout << "Specify the problem " << endl;
		cout << "Insert the number of variables (only for problems different of MoFSSP" << endl;
		cout << "Insert the number of objectives " << endl;
		cout << "Insert specfic oeprator (insert or exchange)" << endl;
		cout << "Insert the number of index of the Flow Shop instance from 001 to 110" << endl;
		cout << "Insert the benchmark for SDST 50 or 125" << endl;
		cout << "Insert the model " << endl;
		cout << "Insert the neighbhorhood size for mating selection" << endl;
		cout << "Insert the type of neigbhorhood, 0=constant, 1=adaptive" << endl;
		cout << "Insert a description" << endl;
		cout << "Insert Run id for the randon seed" << endl;
		cout << "Insert prob mutation" << endl;
		cout << "Insert the adaptive update scheme 0=linear, 1=exponential" << endl;
		cout << "n size update " << endl;
		cout << "maximum n update" << endl;
		cout << "Initial POP 0=without heuristic, 1=heuristic" << endl;
		cout << "Probability crossover for genetic operator" << endl;
		cout << "Shaking procedure 0=no, 1=yes" << endl;
		cout << "Type of LS opereator. PX is prob of interchange" << endl;
		exit (EXIT_FAILURE);

	} else{
		scalarFunction = atoi(argv[1]);
		problem = atoi(argv[2]);
		numbVar = atoi(argv[3]);
		numbObj = atoi(argv[4]);
		pho = argv[5];
		idxInstance = argv[6];
		benchmark = atoi(argv[7]);
		model = atoi(argv[8]); //model
		neighborhoodSizeSelection = atoi(argv[9]); //set the number of neighbor for the mating selection
		typeNeighbor = atoi(argv[10]);
		test = argv[11];
		seed_id = atoi(argv[12]);
		mutation = atof(argv[13]);
		typeSampling = atoi(argv[14]);
		n_size_update=atoi(argv[15]);
		n_sampling=atoi(argv[16]);
		initial_pop=atoi(argv[17]);
		probGenetic=atof(argv[18]);
		shake=atoi(argv[19]);
		probX=atof(argv[20]);

	}

	//----------------------------------------------------------------
	srand((unsigned)time(0)*(seed_id*1000+1)); //seed to random values
	cout << "valor seed " << (unsigned)time(0)*(seed_id*1000+1) << endl;

	//Initialization of the problem
	AbstractProblem * problemTest;
//	if(problem==0) problemTest = new mUBQP(0); //seed 0
//	 else if(problem==1)problemTest = new MOKP();
//	 else if(problem==2)problemTest = new City();
//	 else if(problem==3)problemTest = new Trap5();
//	 else if(problem==4)problemTest = new LOP();
//	 else
	if(problem==5){
		 problemTest = new PFSP();
		 numbObj=2;
	 }
	 else if(problem==6){
		 problemTest = new PFSP_tard();
		 numbObj=2;
	 }
	 else if(problem==7){
		 problemTest = new PFSP_3ob();
		 numbObj=3;
	 }
	 else if(problem==8){
		 problemTest = new PFSP_SDST();
		 numbObj=2;
	 }
	 else if(problem==9){
		 problemTest = new SDST_Cmax_TFT();
		 numbObj=2;
	 }
	else{
		cout << "error in to set the number of the problem" << endl;
		exit (EXIT_FAILURE);
	}
	problemTest->Initialize(); //initialize the vector structure according to the test instance, i.e., according to number of objectives and number of variables
	problemTest->LoadInstance(); //load the data of the test instance
	problemTest->Show(); //print the data


	cout << "shaking " << mutation << " model " << model << " variables " << numbVar << " initial " << initial_pop;

	//getchar();
	//set the number of suproblems according to the number of objectives

	int numbSubProb;
	if(numbObj==2)numbSubProb=101;
	else if (numbObj==3)numbSubProb=231;
	else if (numbObj==4)numbSubProb=220;
	else if (numbObj==6)numbSubProb=126;
	else if (numbObj==8)numbSubProb=120;
	else if (numbObj==10)numbSubProb=220;
	else{
		cout << "Error setting the number of subproblem" << endl;;
		exit (EXIT_FAILURE);
	}
	cout << numbSubProb << endl;

	//The main class is invoked
		CMOEAD * MOEAD;
		MOEAD = new CMOEAD(numbSubProb, neighborhoodSizeSelection, model, problemTest, typeNeighbor, mutation, typeSampling);

	// The main method from the class MOEAD is invoked
		MOEAD->Run();


}
