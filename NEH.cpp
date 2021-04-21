
#include "NEH.h"
#include "NTools.h"
#include "Problem.h"
#include "Global.h"

NEH::NEH()
{
//	cout << "entrou no construtor" << endl; getchar();
	processingTime = new double[numbVar];
	gen = new int[numbVar];


}

NEH::~NEH()
{
}

int * NEH::Improvements(int * permutation, int current_fit, AbstractProblem *m_fsp, int indexFnc){

	int * permut = new int[numbVar];
	int * fit = new int[numbObj];
	memcpy(permut,permutation,sizeof(int)*numbVar);

	int c_fit = current_fit;
	permut = this->JIBIS(permut,c_fit,m_fsp,indexFnc);
	PrintArray(permut,numbVar,"NEH + JIBIS ");
	fit = m_fsp->solutionFitness(permut); PrintArray(fit,numbObj,"F(x): ");

	c_fit = fit[indexFnc];
	permut = this->OSSBIS(permut,c_fit,m_fsp,indexFnc);
	PrintArray(permut,numbVar,"NEH + JIBIS + OSSBIS: ");
	fit = m_fsp->solutionFitness(permut); PrintArray(fit,numbObj,"F(x): ");

	c_fit = fit[indexFnc];
	permut = this->JIBSS(permut,c_fit,m_fsp,indexFnc);
	PrintArray(permut,numbVar,"NEH + JIBIS + OSSBIS + JIBSS: ");
	fit = m_fsp->solutionFitness(permut); PrintArray(fit,numbObj,"F(x): ");

	c_fit = fit[indexFnc];
	permut = this->JIBIS(permut,c_fit,m_fsp,indexFnc);
	PrintArray(permut,numbVar,"NEH + JIBIS + OSSBIS + JIBSS + JIBIS: ");
	fit = m_fsp->solutionFitness(permut); PrintArray(fit,numbObj,"F(x): ");

	return permut;

//	getchar();

}
int * NEH::OSSBIS(int * permutation, int current_fit, AbstractProblem * m_fsp, int indexFnc){

	int * cpy_permutation = new int[numbVar];;
	int * aux_permutation = new int[numbVar];
	int * best_permutation = new int[numbVar];

	memcpy(best_permutation,permutation,sizeof(int)*numbVar);
	memcpy(aux_permutation,permutation,sizeof(int)*numbVar);

	int fitness;
	int best_fitness = current_fit;

	for(int i=0;i<numbVar;i++){
		memcpy(aux_permutation,best_permutation,sizeof(int)*numbVar);

		for(int k=0;k<numbVar;k++){
			if(k!=i){
				memcpy(cpy_permutation,aux_permutation,sizeof(int)*numbVar);

				InsertAt(cpy_permutation,i,k,numbVar);

				if(indexFnc==0)fitness = m_fsp->EvalCmax(cpy_permutation,numbVar);
				else fitness = m_fsp->EvalTWT(cpy_permutation,numbVar);

				if(fitness<best_fitness){
					memcpy(best_permutation,cpy_permutation,sizeof(int)*numbVar);
					best_fitness = fitness;
					//PrintArray(best_permutation,numbVar,"new: ");
				}
			}
		}
	}

	delete [] cpy_permutation;
	delete [] aux_permutation;


	return best_permutation;

}

int * NEH::JIBIS(int * permutation, int current_fit, AbstractProblem * m_fsp, int indexFnc){

	int * aux_permutation = new int[numbVar];
	int * best_permutation = new int[numbVar];

	memcpy(best_permutation,permutation,sizeof(int)*numbVar);
	memcpy(aux_permutation,permutation,sizeof(int)*numbVar);


	int fitness;
	int best_fitness = current_fit;

	for(int i=0;i<numbVar;i++){

		for(int k=0;k<numbVar;k++){
			if(k!=i){

				InsertAt(aux_permutation,i,k,numbVar);

				if(indexFnc==0)fitness = m_fsp->EvalCmax(aux_permutation,numbVar);
				else fitness = m_fsp->EvalTWT(aux_permutation,numbVar);

				if(fitness<best_fitness){
					memcpy(best_permutation,aux_permutation,sizeof(int)*numbVar);

					best_fitness = fitness;
					//PrintArray(best_permutation,numbVar,"new: ");
				}else{
					memcpy(aux_permutation,best_permutation,sizeof(int)*numbVar);

				}
			}
		}
	}

	delete [] aux_permutation;


	return best_permutation;
}

int * NEH::JIBSS(int * permutation, int current_fit, AbstractProblem * m_fsp, int indexFnc){

	int * aux_permutation = new int[numbVar];
	int * best_permutation = new int[numbVar];

	memcpy(best_permutation,permutation,sizeof(int)*numbVar);
	memcpy(aux_permutation,permutation,sizeof(int)*numbVar);


	int fitness;
	int best_fitness = current_fit;

	for(int i=0;i<numbVar;i++){

		for(int k=0;k<numbVar;k++){
			if(k!=i){

				Swap(aux_permutation,i,k);

				if(indexFnc==0)fitness = m_fsp->EvalCmax(aux_permutation,numbVar);
				else fitness = m_fsp->EvalTWT(aux_permutation,numbVar);

				if(fitness<best_fitness){
					memcpy(best_permutation,aux_permutation,sizeof(int)*numbVar);

					best_fitness = fitness;
					//PrintArray(best_permutation,numbVar,"new: ");
				}else{
					memcpy(aux_permutation,best_permutation,sizeof(int)*numbVar);

				}
			}
		}
	}

	delete [] aux_permutation;


	return best_permutation;
}


void NEH::sortProcessingTimesCmax(AbstractProblem * m_fsp)
{
	for(int j=0;j<numbVar;j++){
		processingTime[j]=0;
		gen[j]=j; // gen with the all the jobs
	}

	for(int i=0;i<numbVar;i++){
		for(int k=0;k<machine;k++){
			processingTime[i]+=m_fsp->m_processingtimes[k][i]; //processing time on the machines for each job from 1 to n
		}
	}
	sortDsc();


}

void NEH::sortProcessingTimesCmax_SKE(AbstractProblem * m_fsp)
{
	double * AVG = new double[numbVar];
	double * STD = new double[numbVar];
	double * SKE = new double[numbVar];
	int * ptimes = new int[numbVar];
	double * aux1 = new double[numbVar];
	double * aux2 = new double[numbVar];
	for(int j=0;j<numbVar;j++){
		processingTime[j]=0.0;
		AVG[j]=0.0; STD[j]=0.0; SKE[j]=0.0; ptimes[j]=0; aux1[j]=0.0; aux2[j]=0.0;
		gen[j]=j; // gen with the all the jobs
	}

	cout << "chegou antes aqui " << endl;
	for(int i=0;i<numbVar;i++){
		for(int k=0;k<machine;k++){
			ptimes[i]+=m_fsp->m_processingtimes[k][i]; //processing time on the machines for each job from 1 to n
		}

		AVG[i] = (double)((double)ptimes[i]/(double)machine);

		for(int k=0;k<machine;k++){
			aux1[i] += (double) pow(m_fsp->m_processingtimes[k][i]-AVG[i],2);
			aux2[i] += (double) pow(m_fsp->m_processingtimes[k][i]-AVG[i],3);
		}
		STD[i] = sqrt((double)(double(1/(double)(machine-1))*aux1[i]));
		SKE[i] = (double)(aux2[i]/(double)machine) / (double)(pow((double)(sqrt(aux1[i]/(double)machine)),3));
		cout << SKE[i] << endl;
		processingTime[i] = AVG[i] + STD[i] + SKE[i];
	}

	cout << "chegou aqui " << endl;
	sortDsc();

	delete [] AVG;
	delete [] STD;
	delete [] SKE;
	delete [] ptimes;


}

int * NEH::nehHeuristicCmax(AbstractProblem * m_fsp)
{
	//step 1  order the jobs by non-increasing sums of processing times on the machines
	sortProcessingTimesCmax_SKE(m_fsp);
	PrintArray(processingTime,numbVar,"processing times sortedddddd:     ");
	PrintArray(gen,numbVar,"gen1:     ");

	int bestFitness,fitness;

	//Step 2: evaluate the first two jobs scheduled
	bestFitness = m_fsp->EvalCmax(gen,2);

	int * temp_gen = new int[numbVar];
	memcpy(temp_gen,gen,sizeof(int)*numbVar);
	//swap the first two jobs and evaluate it
	Swap(temp_gen,0,1);
	fitness = m_fsp->EvalCmax(temp_gen,2);
	if(fitness<bestFitness) memcpy(gen,temp_gen,sizeof(int)*numbVar);

	//Step 3 for size=3 to numbVar do Step 4
	for(int size=3;size<numbVar;size++){
		bestFitness = m_fsp->EvalCmax(gen,size);
		for(int pos=0;pos<size;pos++){ //Step 4 insert the kth job a the place, which minimizes the partial evaluation fitness among the k possible ones
			memcpy(temp_gen,gen,sizeof(int)*numbVar);
			InsertAt(temp_gen,size,pos,numbVar);
			fitness = m_fsp->EvalCmax(temp_gen,size);
			if(fitness<bestFitness){
				bestFitness = fitness;
				memcpy(gen,temp_gen,sizeof(int)*numbVar);
			}
		}
	}
	PrintArray(gen,numbVar,"gen4:     ");

	return gen;
}

int * NEH::sortProcessingTimesTWT_EDD(AbstractProblem * m_fsp)
{
	for(int j=0;j<numbVar;j++){
		processingTime[j]=0.0;
		gen[j]=j; // gen with the all the jobs
	}

	for(int i=0;i<numbVar;i++) processingTime[i]=(double)(m_fsp->m_weights[i]*m_fsp->m_DueDates[i]); //EDD with processing time

	sortAsc();

	return gen;


}


int * NEH::sortProcessingTimesTWT_EDDP(AbstractProblem * m_fsp)
{
	for(int j=0;j<numbVar;j++){
		processingTime[j]=0.0;
		gen[j]=j; // gen with the all the jobs
	}

	for(int i=0;i<numbVar;i++){
		for(int k=0;k<machine;k++){
			processingTime[i]+=m_fsp->m_processingtimes[k][i]; //processing time on the machines for each job from 1 to n
		}
	}

	for(int i=0;i<numbVar;i++) processingTime[i]=(double)(m_fsp->m_DueDates[i]/processingTime[i]); //EDD with processing time

	sortAsc();

	return gen;


}

int * NEH::sortProcessingTimesTWT_LPT(AbstractProblem * m_fsp)
{
	for(int j=0;j<numbVar;j++){
		processingTime[j]=0.0;
		gen[j]=j; // gen with the all the jobs
	}

	for(int i=0;i<numbVar;i++){
		for(int k=0;k<machine;k++){
			processingTime[i]+=m_fsp->m_processingtimes[k][i]; //processing time on the machines for each job from 1 to n
		}
	}

	sortDsc();

	return gen;
}

int * NEH::nehHeuristicTWT(AbstractProblem * m_fsp)
{

	int bestFitness,fitness;
	int * sortPermut = new int[numbVar];
	int * bestSortPermut = new int[numbVar];

	sortPermut = sortProcessingTimesTWT_EDDP(m_fsp); //first sort<<<<<<<<<<<<<<<<<<<<<<
	PrintArray(processingTime,numbVar,"EDDP: ");
	PrintArray(sortPermut,numbVar,"gen_EDDP:     ");

	bestFitness = m_fsp->EvalTWT(sortPermut,numbVar);

	memcpy(bestSortPermut,sortPermut,sizeof(int)*numbVar);

	cout << bestFitness << endl;

	sortPermut = sortProcessingTimesTWT_LPT(m_fsp); //second sort<<<<<<<<<<<<<<<<<<<<<<<
//	PrintArray(processingTime,numbVar,"LPT:     ");
	PrintArray(sortPermut,numbVar,"gen_LPT:     ");
	fitness = m_fsp->EvalTWT(sortPermut,numbVar);
	cout << bestFitness << " " << fitness << endl;
	if(fitness < bestFitness){
		memcpy(bestSortPermut,sortPermut,sizeof(int)*numbVar);
		bestFitness = fitness;
	}
	sortPermut = sortProcessingTimesTWT_EDD(m_fsp); //third sort<<<<<<<<<<<<<<<<<<<<<<<<<
//	PrintArray(processingTime,numbVar,"EDD     ");
	PrintArray(sortPermut,numbVar,"gen_EDD:     ");
	fitness = m_fsp->EvalTWT(sortPermut,numbVar);
	cout << bestFitness << " " << fitness << endl;
	if(fitness < bestFitness){
		memcpy(bestSortPermut,sortPermut,sizeof(int)*numbVar);
		bestFitness = fitness;
	}

	PrintArray(bestSortPermut,numbVar,"Best permut:     "); //best sort chosen

	//Step 2: after select the best sort, according to the TWT, evaluate the first two jobs scheduled
	bestFitness = m_fsp->EvalTWT(bestSortPermut,2);
	int * temp_gen = new int[numbVar];
	memcpy(temp_gen,bestSortPermut,sizeof(int)*numbVar);
	//swap the first two jobs and evaluate it
	Swap(temp_gen,0,1);
	fitness = m_fsp->EvalTWT(temp_gen,2);
	if(fitness<bestFitness) memcpy(bestSortPermut,temp_gen,sizeof(int)*numbVar);

	//Step 3 for size=3 to numbVar do Step 4
	for(int size=3;size<numbVar;size++){
		bestFitness = m_fsp->EvalTWT(bestSortPermut,size);
		for(int pos=0;pos<size;pos++){ //Step 4 insert the kth job a the place, which minimizes the partial evaluation fitness among the k possible ones
			memcpy(temp_gen,bestSortPermut,sizeof(int)*numbVar);
			InsertAt(temp_gen,size,pos,numbVar);
			fitness = m_fsp->EvalTWT(temp_gen,size);
			if(fitness<bestFitness){
				bestFitness = fitness;
				memcpy(bestSortPermut,temp_gen,sizeof(int)*numbVar);
			}
		}
	}
	PrintArray(bestSortPermut,numbVar,"gen_Final:     ");

	delete [] sortPermut;

	return bestSortPermut;
}

void NEH::sortAsc(){

	//begin Sort the processingTimes in a decreasing order
	int i, flag = 1;
	double tempProcess;
	double tempIndex;
	int d = numbVar;
	while( flag || (d > 1)){  // boolean flag (true when not equal to 0)
	flag = 0; // reset flag to 0 to check for future swaps
	d = (d+1) / 2;
		for (i = 0; i < (numbVar - d); i++){ //Earliest Due Dates (EDD)
			if (processingTime[i+d] < processingTime[i]){
			  tempProcess = processingTime[i+d]; // swap positions i+d and i
			  tempIndex = gen[i+d];

			  processingTime[i+d]=processingTime[i];
			  gen[i+d]=gen[i];

			  processingTime[i]=tempProcess;
			  gen[i]=tempIndex;
			  flag = 1; // tells swap has occurred
			}
		}
	} //end shellsort


}

void NEH::sortDsc(){

	//begin Sort the processingTimes in a decreasing order
	int i, flag = 1;
	double tempProcess;
	double tempIndex;
	int d = numbVar;
	while( flag || (d > 1)){  // boolean flag (true when not equal to 0)
	flag = 0; // reset flag to 0 to check for future swaps
	d = (d+1) / 2;
		for (i = 0; i < (numbVar - d); i++){ //Earliest Due Dates (EDD)
			if (processingTime[i+d] > processingTime[i]){
			  tempProcess = processingTime[i+d]; // swap positions i+d and i
			  tempIndex = gen[i+d];

			  processingTime[i+d]=processingTime[i];
			  gen[i+d]=gen[i];

			  processingTime[i]=tempProcess;
			  gen[i]=tempIndex;
			  flag = 1; // tells swap has occurred
			}
		}
	} //end shellsort

}

