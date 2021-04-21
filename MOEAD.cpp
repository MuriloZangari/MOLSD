/////////////////////////////////////////////
//MOEAD Framework for discrete multiobjective optimization problems
////////////////////////////////////////////
#include <algorithm>
#include "MOEAD.h"
#include "Global.h"
#include "NTools.h"



using namespace std;


CMOEAD::CMOEAD(int n_sub, int neighorSizeSel, int modell, AbstractProblem* problemTestt, int typeNeighborhood, float mut, int ts) {
	numbSubProb = n_sub;
	numbGroup = numbSubProb;

	maxNumbOfUpdate = n_sampling;
	size_selection = neighorSizeSel;
	//neighborhoodSizeSelection = neighorSizeSel;
	neighborhoodSizeUpdate = n_size_update; //to use the different types of update mechanism, the neighborhood to update is set to all solution (subproblems)
	numbSamplingSol = 1; //(int)(numbSubProb/numbGroup); //to a fair comparison, each generation (cycle) generates (N = numSubProb) new solutions as conventional MOEAD algorithms
	update =0; //0 means conventional update according to the neighbors, 1 means means global replacement in a random way, 2 means global replacement most appropriate subproblem (according to the scalar function) and its neighbors, 3 means subregion approapriate
	ScalarFunction = scalarFunction; // 0 means weighted sum, and 1 means theebicheffi
	model = modell;
	problemTest = problemTestt;

	solution= new Subprob[numbSubProb];
	auxSolution = new Subprob[1000];
	group = new Group[numbGroup];
	lambdaVector =new Lambda[numbSubProb];
	externalPareto = new ArchivePareto;
	invertedBB = new int[numbVar];
	sol1 = new int[numbVar];
	sol2 = new int[numbVar];
	compositionn = new int[numbVar];

	previousFitnessVale = new int *[numbSubProb];


	m_restart_countdown = numbVar;

	m_max_restars = new int[numbSubProb];
	m_count_down_for_restart = new int[numbSubProb];
	for(int i=0;i<numbSubProb;i++){
		previousFitnessVale[i] = new int[numbObj];
		m_count_down_for_restart[i]=m_restart_countdown;
		m_max_restars[i]=10000;
	}

	//ONLY MOEA/D
	numbGenerations = numbVar*500;
	moead_maxEval = (int)(numbSubProb*numbGenerations);
	m_MOEAD_restar_shake_power = 14;
	//m_VNS_restart_shake_power = 0;
	cout << "only MOEA/D" << endl;

	MAX_EVAL= moead_maxEval;
	cout << "max evaluations: " << MAX_EVAL << endl;
	cout << "the problem index is " << problem << endl;

	t_neighborhood=typeNeighborhood;
	boltz = 0;
	mutation = mut;
	typeSampling = ts;
	strcpy(type,pho.c_str());
	//getchar();

}

CMOEAD::~CMOEAD() {
	delete [] solution;
	delete problemTest;

	delete [] invertedBB;
	delete [] compositionn;
	delete [] sol1;
	delete [] sol2;

}

void CMOEAD::InitializeReferencePoint() { //Initialization of the reference point depends on the problems be maximization or minimization
	ReferencePoint = new int[numbObj];
	worstReferencePoint = new int[numbObj];

	if(problem==0 or problem==1 or problem==3 or problem==4){ //maximization - problems: mUBQP, MOKP, LOP
		for(int k=0;k<numbObj;k++){
			ReferencePoint[k]= intMin;
			worstReferencePoint[k] = intMax;
			cout << ReferencePoint[k] << endl;
		}
	}else if(problem==2 or problem==5 or problem==6 or problem==7 or problem==8 or problem==9){ //minimization - problems: MTSP
		for(int k=0;k<numbObj;k++){
			ReferencePoint[k]=intMax; //max, because its is a minimization problem
			worstReferencePoint[k]=intMin;
			cout << ReferencePoint[k] << endl;
		}
	}else{
		cout << "error function InitializeReferencePoint to set the problem to the reference point " << endl;
	}
}

void CMOEAD::InitializeWeightVectors(int obj, int numSubProb) { //Load the weight vectors
	char stringPeso[10];
		ifstream ler;
		string cu = static_cast<ostringstream*>( &(ostringstream() << obj) )->str();
		string sub = static_cast<ostringstream*>( &(ostringstream() << numSubProb) )->str();
		ler.open(("../weightsIshi/wks" + sub + cu + ".dat").c_str());
		if(ler==NULL){
			cout << "Error opening file" << endl;
			exit (EXIT_FAILURE);
		}
		ler >> stringPeso;
		if(atof(stringPeso)!= numSubProb){
//			cout << "Number of weights is different of number of subproblems" << endl;;
			exit (EXIT_FAILURE);
		}
		for(int i=0;i<numSubProb;i++){
			for(int j=0;j<obj;j++){
				ler >> stringPeso;
				lambdaVector[i].setWeight(atof(stringPeso),j);
//				cout << lambda[i].getWeight(j) << " ";
			}
//			cout << endl;
		}

}

double CMOEAD::EuclideanDist(Lambda lambdaAtual, Lambda lambda2, int obj){ //Calculate the distance between the weight vectors
	double dist =0;
	for(int k=0; k<obj; k++){
		double diff = lambdaAtual.getWeight(k) - lambda2.getWeight(k);
		dist+= diff*diff;
	}

	return sqrt(dist);
}

void CMOEAD::shellSort(Distance* dist, int solutionNumber){
	int i, flag = 1;
	Distance temp;
	int d = solutionNumber;
	while( flag || (d > 1)){  // boolean flag (true when not equal to 0)
	flag = 0; // reset flag to 0 to check for future swaps
	d = (d+1) / 2;
		for (i = 0; i < (solutionNumber - d); i++){
			if (dist[i+d].getDist() < dist[i].getDist()){
			  temp = dist[i+d]; // swap positions i+d and i
			  dist[i+d]=dist[i];
			  dist[i]=temp;
			  flag = 1; // tells swap has occurred
			}
		}
	}

}

void CMOEAD::InitializeNeighorhood(int obj, int numSubProb) { //Set the index of neighbors for each subproblem according to the neighborhood size and the euclidean distance
	//float uniformDist = constantDistanceParameter;
	for(int pop=0;pop<numSubProb; pop++){
		Distance dist[numSubProb];
		Lambda lambdaAtual = lambdaVector[pop];
		Lambda lambda2;
		for(int pop2=0;pop2<numSubProb;pop2++){ 	//calcule the distance between all population (subproblems)
			lambda2 = lambdaVector[pop2];
			dist[pop2].setDist(this->EuclideanDist(lambdaAtual,lambda2, obj));
			dist[pop2].setIndiceViz(pop2);
		}
		this->shellSort(dist,numSubProb); 			//ordern the array dist and set the neighbors according to the closest neighbors
//		cout << "list for " << pop << " ";
//		for(int j=0;j<constantSizeSelection;j++){
//			cout << dist[j].getDist() << " and index " << dist[j].getIndiceViz() << "; ";
//
//		}
//		cout << endl;
		lambdaVector[pop].NeigbhborSizeSimga=size_selection;


		for(int i=0;i<numSubProb;i++){
			lambdaVector[pop].setIndexOfNeighbors(dist[i].getIndiceViz(),i);
//				cout << lambdaVector[pop].getIndexOfNeighbors(i) << " ";
			}
//			cout << endl;

		if(t_neighborhood==0) lambdaVector[pop].NeighborSizeUpdate= neighborhoodSizeUpdate;
		else if(t_neighborhood==1) lambdaVector[pop].NeighborSizeUpdate=1;
		else{
			cout << "Error setting the type of neighborhood for updating" << endl;;
			exit (EXIT_FAILURE);
		}

		lambdaVector[pop].numbSampling=numbSamplingSol;
//		cout << pop << " " << lambdaVector[pop].NeigbhborSizeSimga << endl;
		}
}

void CMOEAD::InitializeGroups(Group* group, int model) { //Here we initialize a model for each SUBPROBLEM (could be for a group of subproblems but we have not used this scenario)
	int distCentroides = (int) (numbSubProb/numbGroup);
	for(int g=0;g<numbGroup;g++){
		if(g==0){
			group[0].setCentroid((int) distCentroides /2);
		}else{
			group[g].setCentroid(group[g-1].getCentroid() + distCentroides);
		}
//		group[g].setNumbOfSampling(numbSamplingSol);
//		group[g].setNumbOfSolutionGroup(1);
		char type[2];
		strcpy(type,pho.c_str());
//		//initialize the model of each group
//		cout << "model " << model << " problem " << problem << endl;
			 if (model==0 and (problem==5 or problem==6 or problem==7 or problem==8 or problem==9))group[g].modelGroup = new GA(numbSubProb, mutation);
		//else if (model==9 and (problem==5 or problem==6 or problem==7 or problem==8 or problem==9))group[g].modelGroup = new MALLOWS(numbSubProb,type);
		else if (model==10 and (problem==5 or problem==6 or problem==7 or problem==8 or problem==9))group[g].modelGroup = new LS(numbSubProb);
		else if (model==11 and (problem==5 or problem==6 or problem==7 or problem==8 or problem==9))group[g].modelGroup = new GLS(probGenetic,mutation, probX);

		else{
			cout << "Error setting the model for the problem" << endl;;
			exit (EXIT_FAILURE);
		}

		group[g].modelGroup->InitializeModel();
//		group[g].modelGroup->ShowModel();
//		cout << endl;

	}
}

void CMOEAD::InitializeSolutions() {

	cout << "initi " << endl;
//	int gen[20]={3, 15, 17, 9, 8, 13, 1, 16, 10, 19, 14, 2, 4, 7, 11, 6, 20, 12, 5, 18};

	if(problem==5 or problem==6 or problem==7 or problem==8 or problem==9){

		if(initial_pop==0){
			cout << "initialization " << endl;
			for(int j=0;j<numbSubProb;j++){
				vector<int> initialPermutation;
				for(int i=0;i<numbVar;i++){
					initialPermutation.push_back(i);
				}
				random_shuffle(initialPermutation.begin(),initialPermutation.end());
				for(int i=0;i<numbVar;i++){
					solution[j].setVectorSolution(initialPermutation.at(i),i);
					cout << initialPermutation.at(i) << " ";
				}
				cout << endl;
				problemTest->solutionFitness(solution[j]);
				//solution[j].PrintSolution();
				this->UpdateReferencePoint(solution[j]);
			}

		}else{
					cout << "NEH nao apta no momento" << endl; exit (EXIT_FAILURE);
					//first, create a solution using heuristics
					int * genesNEH = new int[numbVar];
					int * genesNEHedd = new int[numbVar];
					int * fit = new int[numbObj];
					//int singleFit;

					NEH m_nehHeuristic;

					genesNEH = m_nehHeuristic.nehHeuristicCmax(problemTest);
					PrintArray(genesNEH,numbVar,"NEH: ");
					fit = problemTest->solutionFitness(genesNEH); PrintArray(fit,numbObj,"F(x): "); //getchar();

					genesNEH = m_nehHeuristic.Improvements(genesNEH, fit[0],problemTest,0); //getchar();

					cout << "end heuristic Cmax" << endl;

					genesNEHedd = m_nehHeuristic.nehHeuristicTWT(problemTest);
					PrintArray(genesNEHedd,numbVar,"NEHedd ");
					fit = problemTest->solutionFitness(genesNEHedd); PrintArray(fit,numbObj,"F(x): "); //getchar();

					genesNEHedd = m_nehHeuristic.Improvements(genesNEHedd,fit[1],problemTest,1); //getchar();

					solution[0].SetSoluion(genesNEHedd);
					solution[0].setIndex(0);
					problemTest->solutionFitness(solution[0]);
					this->UpdateReferencePoint(solution[0]); //update the reference point for each solution initialized

					solution[(numbSubProb-1)].SetSoluion(genesNEH);
					solution[(numbSubProb-1)].setIndex((numbSubProb-1));
					problemTest->solutionFitness(solution[(numbSubProb-1)]);
					this->UpdateReferencePoint(solution[(numbSubProb-1)]); //update the reference point for each solution initialized

					//getchar();

					vector<int> initialPermutation;
					for(int i=0;i<numbVar;i++){
						initialPermutation.push_back(i);
					}

					for(int j=1;j<(numbSubProb-1);j++){
						int * cpygene = new int[numbVar];
						if(j%2==0){
							memcpy(cpygene,genesNEH,sizeof(int)*numbVar);

						}else{
							memcpy(cpygene,genesNEHedd,sizeof(int)*numbVar);
						}

						Shake_Insert(cpygene,m_MOEAD_restar_shake_power);

						int index = j;
						solution[index].setIndex(index);
						solution[index].SetSoluion(cpygene);
						problemTest->solutionFitness(solution[index]);
						this->UpdateReferencePoint(solution[index]); //update the reference point for each solution initialized
						cout << "index shake H " << index << ": ";solution[index].PrintSolution();
						delete [] cpygene;

//						}else{
//
//							random_shuffle(initialPermutation.begin(),initialPermutation.end());
//							int indexR=j;
//							for(int i=0;i<numbVar;i++){
//								solution[indexR].setIndex(indexR);
//								solution[indexR].setVectorSolution(initialPermutation.at(i),i);
//							}
//
//							problemTest->solutionFitness(solution[indexR]);
//							this->UpdateReferencePoint(solution[indexR]);
//							cout << "index random " << indexR << ": "; solution[indexR].PrintSolution();
					}
					delete [] genesNEHedd;
					delete [] genesNEH;
					delete [] fit;
		}
	}
}

void CMOEAD::UpdateReferencePoint(Subprob solution) { //Update the reference point according to the type of problem (minimization or maximization
//	if(problem==0 or problem==1 or problem==3 or problem==4){ //maximization - problems: 0 - mUBQP, MOKP
//		for(int obj=0;obj<numbObj;obj++){
//			if(solution.getFunction(obj)>ReferencePoint[obj]) ReferencePoint[obj] = solution.getFunction(obj);
//			if(solution.getFunction(obj)<worstReferencePoint[obj]) worstReferencePoint[obj] = solution.getFunction(obj);
//
//		}
//	}else
//	if(problem ==2 or problem==5 or problem==6 or problem==7 or problem==8 or problem==9){ //minimization - problems: 1 - MOTSP, 5 - MFSSP
		for(int obj=0;obj<numbObj;obj++){
			if(solution.getFunction(obj)<ReferencePoint[obj]) ReferencePoint[obj] = solution.getFunction(obj);
			if(solution.getFunction(obj)>worstReferencePoint[obj]) worstReferencePoint[obj] = solution.getFunction(obj);

		}
//	}else {
//		cout << "Error to set the problem to update reference point" << endl;
//		exit (EXIT_FAILURE);
//	}
}

void CMOEAD::UpdateNeighboringSolution(Subprob newSol, int index) { //update 0 - Conventional update of the subproblems
	int count=0; 	//iterator of max number of update, if achieve the max number of update, then stop the update
//	cout << "centroid for update " << centroid << endl;

	for(int n=0;n<lambdaVector[index].NeighborSizeUpdate;n++){
		double f1, f2;

		int id = lambdaVector[index].getIndexOfNeighbors(n);
		f1 = newSol.ScalarFitnessValueNormalized(lambdaVector[id],ScalarFunction);
		f2 = solution[id].ScalarFitnessValueNormalized(lambdaVector[id],ScalarFunction);
//		cout << f1 << " - " << f2 << endl;
//		cout << newSol.getScalarValue() << " - " << solution[id].getScalarValue() << endl;
//		solution[id].PrintSolution();
		if(f1<=f2){ //minimization, verify this operator
			solution[id]=newSol;
			solution[id].setIndex(id);
			count++;
//			solution[id].PrintSolution();
		}
		if(count>=maxNumbOfUpdate){

			break;
		}
	}
}

void CMOEAD::UpdateNeighboringSolutionGR1(Subprob newSol, int index) { //update 1 (alternative version)
	//iterator of max number of update, if achieve the max number of update, then stop the update
//	cout << "centroide para atualizacao " << centroid << endl;

	random_shuffle(vectorOfNeighbors.begin(),vectorOfNeighbors.end());
	int count=0;

	for(int n=0;n<numbSubProb;n++){ //update all the solution random way
		double f1, f2;
		int id = vectorOfNeighbors.at(n);
		f1 = newSol.ScalarFitnessValueNormalized(lambdaVector[id],ScalarFunction);
		f2 = solution[id].ScalarFitnessValueNormalized(lambdaVector[id],ScalarFunction);
		if(f1<=f2){
//			cout << "solution generated " << newSol.getIndex() << " update " << id << endl;
			solution[id]=newSol;
			solution[id].setIndex(id);
			count++;
		}
		if(count>=maxNumbOfUpdate){
			break;
		}
	}
}
void CMOEAD::UpdateSelfish(Subprob newSol, int id){ // Alternative version of update

		double f1, f2;

		f1 = newSol.ScalarFitnessValueNormalized(lambdaVector[id],ScalarFunction);
		f2 = solution[id].ScalarFitnessValueNormalized(lambdaVector[id],ScalarFunction);
//		cout << f1 << " - " << f2 << endl;

		if(f1<=f2){ //minimization, verify this operator
			solution[id]=newSol;
			solution[id].setIndex(id);
		}
}

void CMOEAD::UpdateNeighborhoodSize(){

	int adaptiveSizeUpdate;

	//linear or sigmoid

	if(typeSampling==0) adaptiveSizeUpdate = ceil((double)cycle * (double)numbSubProb / (double)numbGenerations);
	else if(typeSampling==1) adaptiveSizeUpdate = ceil((double)((double)(exp((5.0 * (double)cycle)/(double)numbGenerations)-1.0)*(double)numbSubProb)/(double)(exp(5.0)-1.0));
	else{
		cout << "Error setting the adaptive update" << endl;;
		exit (EXIT_FAILURE);
	}
	cout << " adptive " << adaptiveSizeUpdate;
	for(int n=0;n<numbSubProb;n++){
		lambdaVector[n].NeighborSizeUpdate=adaptiveSizeUpdate;
	}

}



void CMOEAD::Run() { // Main method

	struct timeval inicio, final; //execution time
	gettimeofday(&inicio, NULL);

	this->InitializeReferencePoint();							//initialize the reference point according to if the problem is maximization or minimization

	this->InitializeWeightVectors(numbObj,numbSubProb); //initialize the weight vectors
	this->InitializeNeighorhood(numbObj,numbSubProb);	//initialize neighborhood

	problemTest->heuristicValue(lambdaVector, numbSubProb); //an information of the problem if saved according to the weighted vectors (subproblem) if it is necessary for the problem, such as MOKP that needs a heuristic value for greedy repair, this is called once
	//int numGRinvoked=0; int maxGreedRepair=	numbVar*numbSubProb; //initialize the counter of number of greed repair invoked
	this->InitializeSolutions();			//Initialize the solutions of each subpproblem in a random way
	//this->ShowSolutions();

	//declaration of the archive of pareto
	for(int i=0;i<numbSubProb;i++) externalPareto->UpdateExternalPareto(solution[i]); 	cout << "current size EP " << externalPareto->sizeEP  << endl; //update the external pareto with the initialized solutions

	this->InitializeGroups(group, model); //centroid is the reference subproblem to a group. In this case it is the same index for the subproblem, because we are not using more than onw suprobem for each group

	numbEval=0;
	cycle=0;
	int * totalShaking = new int[numbSubProb]; for(int s=0;s<numbSubProb;s++) totalShaking[s]=0;
	int plot_gen=0;
//	int starts_explore= (int)(numbGenerations*0.0);

//	ofstream fitGen;
//	fitGen.open(("FitGen/"+ test+".dat").c_str());

	for(int i=0;i<numbSubProb;i++){

		vectorOfNeighbors.push_back(i);
	}

	//double * thetaLocal = new double[numbSubProb]; for(int g=0;g<numbSubProb;g++) thetaLocal[g]=theta;

	while (cycle<numbGenerations){	//while a max numb of generation
		cycle++;

		for(int g=0;g<numbGroup;g++){ 		//for each group, construct a model (learning) (or update the model) and sampling a specific number of new solutions
			int centroid = g; 				//index of the centroid of the group, which is used to set the population of the group according to the neighbors
			//////////////LEARNING////////////////////
		//	cout << "subprob " << g << endl;
			group[centroid].modelGroup->Learning(lambdaVector[centroid],lambdaVector[centroid].numbSampling);

			/////////////SAMPLING////////////////////

			int numbSamp = lambdaVector[centroid].numbSampling;
			Subprob newSol[numbSamp];
			for(int s=0;s<numbSamp;s++)
			{	//SAMPLING/generate new solutions according to the parameter numbOfSampling (numbSubProb = numbGroup + numbSampling)

				//bool sampleNew=false;
				//int count = 0;
				group[centroid].modelGroup->Sampling(solution, s, newSol[s]);

				//cout << "the new solution AFTER sampling " << endl;

				newSol[s].setIndex(centroid);
//	        	if (((string)pho)=="k"){
//	        		problemTest->InvSolutionFitness(newSol[s]);
//	            }else {
	            	problemTest->solutionFitness(newSol[s]);
//	            }
	        	numbEval++;
				//cout << "samp " << s << " ";  newSol[s].PrintSolution();

				this->UpdateReferencePoint(newSol[s]); 	//update the reference point for each new solution generated

				/////////DIFFERENT UPDATE MECHANISM FOR EACH NEW SOLUTION SAMPLED/////////////////////
				this->UpdateNeighboringSolution(newSol[s],centroid); 	//update the neighbors according to the centroid of the group, each new solution can update until achieve the parameter maxNumberOfUpdate
//				else if(update==1) {
//					this->UpdateSelfish(newSol[s], centroid);
//					this->UpdateNeighboringSolutionGR1(newSol[s],centroid); 		//Global Replacement (GR), random choose
//
//				}
//				else if (update==2) this->UpdateNeighboringSolutionGR2(lambdaVector,newSol,solution); 		//Global Replacement (GR), most appropriate subproblem, than update its neighbors
//				else if(update==3) this->UpdateNeighboringSolutionDD(centroid,lambdaVector,newSol,solution); //DD, the new solution is associated to the most appropriate subregion, than update the neighborhood of subproblem associate to this subregion
//				else{
//					cout << "Error setting the update mechanism" << endl;
//					exit (EXIT_FAILURE);
//				}
				//////////////////////////////////////////////////////////////////////////////////

			}
			if(shake==1){ //check the evolution to execute the shaking procedure

				int * prev = previousFitnessVale[centroid];
				int * current = solution[centroid].getPointerToFitness();

				if(memcmp(current,prev,sizeof(int)*numbObj)==0 && m_max_restars[centroid]>0){
	//				cout << "index " << centroid << " ";
	//				for(int k=0;k<numbObj;k++) cout << current[k] << " " << prev[k] << " - ";
	//				cout << "count " << m_count_down_for_restart[centroid] << endl;
					if(m_count_down_for_restart[centroid] == 0){
						m_count_down_for_restart[centroid] = m_restart_countdown;
						m_max_restars[centroid]--;
	//					if(m_max_restars[centroid]==0){
	//						cout << "max restarts " << centroid << endl;
	//						getchar();
	//					}
						//cout << "shakeeee !!!" << endl;
						//cout << "vou dar um shake" << endl;

						ShakeSolution(centroid);
						totalShaking[centroid]++;

					}else{
//						if(memcmp(current,prev,sizeof(int)*numbObj)==0){
//							thetaLocal[centroid] = thetaLocal[centroid] - (1/m_restart_countdown);
//						}
						m_count_down_for_restart[centroid]--;
					}

				}else{
					m_count_down_for_restart[centroid]=m_restart_countdown;
//					thetaLocal[centroid] = theta;
				}
				for(int k=0;k<numbObj;k++) previousFitnessVale[centroid][k]=solution[centroid].getFunction(k);


				//solution[centroid].PrintSolution();
			}

		}
		for(int i=0;i<numbSubProb;i++) externalPareto->UpdateExternalPareto(solution[i]); 	cout << "EP " << externalPareto->sizeEP;		//update the external pareto with the updated set of solutions

		if(t_neighborhood==1) UpdateNeighborhoodSize();

//		PRINTING SOME EXTRA INFORMATIONS AFTER EACH GENERATION if used FLOW SHOP

//		if((model==10) and (cycle  % (int)(numbVar*100) ==0) and (cycle != numbGenerations)){
//			plot_gen++;
//			externalPareto->ShowExternalPareto(plot_gen);
//			//this->ShowExtraInformation();
//		}


		cout << " cycle " << cycle << " " << "eval " << numbEval << endl;
		//END PRINTING SOME EXTRA INFORMATIONS

	} ///END WHILE

	cout << "cycle " << cycle << " " << "num evaluations  MOEAD " << numbEval << endl;

	this->ShowExtraInformation(totalShaking);


	this->ShowSolutions();
	cout << "end run " << run_id << endl;
	gettimeofday(&final, NULL);
	double tmili = (1000* (final.tv_sec - inicio.tv_sec) + (final.tv_usec - inicio.tv_usec)/1000);
	this->ShowExecutionTime(tmili,1);								//print the execution time in seconds
	//this->PrintObjSpaceSol(1); 							//print the archive with the final solutions
	externalPareto->ShowExternalPareto(10);				//print the archive with the external pareto solutions
//	this->ShowVar(1);
	cout << "model " << model <<endl;

}

void CMOEAD::ShakeSolution(int index){
    int * inverted= new int[numbVar];
//    Invert(solution[index].getPointerToSolution(), numbVar, inverted);
    memcpy(inverted,solution[index].getPointerToSolution(),sizeof(int)*numbVar);
    int i,j;

    int min_range,max_range,val;
    int half_orbit=5;
    for (int iter=0;iter<m_MOEAD_restar_shake_power;iter++)
    {
        //permute randomly genes in position i and j to scape the stackeness in both neighborhood.
        i = rand() % numbVar;
        //j = rand() % m_size;

        min_range=MAX(i-half_orbit,0);
        max_range=MIN(i+half_orbit,numbVar-1);
        val = rand() % (max_range-min_range);
        j= min_range+val;

        InsertAt(inverted,i,j,numbVar);
    }
//    int * aux = new int[numbVar];
//    Invert(inverted, numbVar, aux);
    solution[index].SetSoluion(inverted);
    problemTest->solutionFitness(solution[index]);

    delete[] inverted;
//    delete[] aux;
}

void CMOEAD::ShowSolutions() {
	for(int i=0;i<numbSubProb;i++){
		cout << "sub " << i << ": ";
		for(int j=0;j<numbObj;j++){
			cout << solution[i].getFunction(j) << " ";
		}
		cout  << solution[i].getScalarValue() << endl;
	}
}
void CMOEAD::ShowVar(int i) {
	ofstream final;
	final.open(("VAR/"+ test+".dat").c_str());

	for(int i=0;i<numbSubProb;i++){

		for(int j=0;j<numbVar;j++){
			final << solution[i].getVectorSolution(j) << " ";
		}

		final << endl;
	}
	final.close();
}

void CMOEAD::PrintObjSpaceSol(int i) {

	ofstream FUN;
	FUN.open(("FUN/"+ test+".dat").c_str());
	for(int i=0;i<numbSubProb;i++){
		for(int obj=0;obj<numbObj;obj++){
			FUN << solution[i].getFunction(obj) << " ";
		}
		FUN << endl;
	}
	FUN.close();
}

void CMOEAD::ShowExecutionTime(double tmili, int i) {

	ofstream TIME;
	TIME.open(("TIME/"+ test+".dat").c_str());
	TIME << tmili/1000 << endl;
	cout << "execution time: " << tmili/1000 << " seconds"<< endl; //print the time
	cout << "end run " << run_id << endl;

}

void CMOEAD::ShowExtraInformation(int * sh) {

	if(problem==5 or problem==6 or problem==7 or problem==8 or problem==9){

	    ofstream distanceVar;
	    distanceVar.open(("DistanceVar/"+test+".dat").c_str());

		float dist,sumDist=0.0;
		for(int i=0;i<numbSubProb-1;i++){
			sol1 = solution[i].getPointerToSolution();
			sol2 = solution[i+1].getPointerToSolution();
			Invert(sol2,numbVar,invertedBB);
		    Compose(sol1, invertedBB, compositionn, numbVar);
		    dist=0.0;
		    if(((string)pho)=="k"){
			    int * vetor = new int[1000];
		        vVector_Fast2(vetor,compositionn,numbVar);

		        for (int j = 0; j < numbVar-1; j++){
		        	dist += (float)vetor[j];
		        }
		        delete [] vetor;
		    }else{

			    dist=CalculateDistanceAndX2(compositionn, NULL);
		    }

		    //sumDist=sumDist+dist;
		    distanceVar << dist << endl;

		}

	    //distanceVar << sumDist/(numbSubProb-1) << endl;
	    distanceVar.close();

	    ofstream totalSH;
	    totalSH.open(("totalSH/"+test+".dat").c_str());
	    for(int i=0;i<numbSubProb;i++){
	    	totalSH << sh[i] << endl;
	    }
	    totalSH.close();

	}
}
