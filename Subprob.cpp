#include "Subprob.h"

int Subprob::getIndex(void) {
	return index;
}

void Subprob::setIndex(int ind){
	index = ind;
}

int Subprob::getFunction(int i) {
	return function[i];
}

void Subprob::setFunction(int fc, int i) {
	function[i]=fc;
}

int Subprob::getVectorSolution(int i) {
	return vectorSolution[i];
}

void Subprob::setVectorSolution(int sol, int i) {
	vectorSolution[i]=sol;
}

double Subprob::getScalarValue(void) {
	return scalarValue;
}

void Subprob::setScalarValue(double sc) {
	scalarValue = sc;
}

bool Subprob::getNoVisited(int i) {
	return noVisited[i];
}

void Subprob::setNoVisited(bool nv, int i) {
	noVisited[i] = nv;
}
int Subprob::getWeigthAcc(int i) {
	return weigthAcc[i];
}

void Subprob::setWeightAcc(int ww, int i) {
	weigthAcc[i]=ww;
}


double Subprob::ScalarFitnessValue(Subprob &solution, Lambda lambda, int scalar) {


//	if(scalar==0){ //weighted sum normalizded ------------------------------------
//		this->setScalarValue(0.0);
//		for(int k=0;k<numbObj;k++){
//			double s = lambda.getWeight(k)*((double)this->getFunction(k));
//
//			this->setScalarValue(this->getScalarValue() + s);
//		}
//
//	}
////	else if(scalar==01){ //weight sum alternative same as MOLSD
////		this->setScalarValue(0.0);
////		for(int k=0;k<numbObj;k++){
////			this->setScalarValue(this->getScalarValue()+ lambda.getWeight(k)*(this->getFunction(k)/ReferencePoint[k]));
////		}
////	}
////	else if(scalar==02){ //weighted sum normalizded and reference point * constant------------------------------------
////		this->setScalarValue(0.0);
////		for(int k=0;k<numbObj;k++){
////			double s = lambda.getWeight(k)*(fabs((double)this->getFunction(k)-(0.6*ReferencePoint[k]))/scale[k]);
////
////			this->setScalarValue(this->getScalarValue() + s);
////		}
////
////	}
//	else if(scalar==1){ //thebycheff
//		this->setScalarValue(intMin); //techycheff
//		for(int k=0;k<numbObj;k++){
//			double dif = (this->getFunction(k) - 1.1*ReferencePoint[k]);
//			dif = fabs(dif);
//			double s = lambda.getWeight(k)*dif;
//			if(s > this->getScalarValue()) this->setScalarValue(s);
//		}
//	}
//	else if(scalar==2){ //normalized tchebbycheff same as MOEADD
//
//		this->setScalarValue(intMin);
//		for(int k=0;k<numbObj;k++){
//			double dif = fabs((this->getFunction(k)- ReferencePoint[k])); //1.1*ReferencePoint[k] - this->getFunction(k);
//			double s;
//
//			 s = dif * lambda.getWeight(k);
//			if(s > this->getScalarValue()) this->setScalarValue(s);
//		}
//
//	}
//	else if(scalar==3){ //PBI
//		double theta=0.1;
//		vector <double> nambda(numbObj);
//		for(int k=0;k<numbObj;k++){
//			nambda[k]=lambda.getWeight(k);
//		}
//		//normalize the weight vector (line segment)
//		double nd = Norma(nambda);
//		for(int k=0;k<numbObj;k++){
//			nambda[k] = lambda.getWeight(k)/nd;
//		}
//		vector <double> realA(numbObj);
//		vector <double> realB(numbObj);
//		//difference between reference point and current point (maximization)
//		for(int k=0;k<numbObj;k++){
//			realA[k] = ((double)ReferencePoint[k] - this->getFunction(k));
//		}
//		//distance along the line segment
//		double d1 = fabs(innerProduct(realA, nambda));
//		//distance to the line segment
//		for(int k=0;k<numbObj;k++){
//			realB[k]=(this->getFunction(k) - ((double)ReferencePoint[k] - d1*nambda[k]));
//		}
//		double d2 = Norma(realB);
//		double pbi = (double)(d1 + theta*d2);
//		this->setScalarValue(pbi);
//
//	}
//	else{
		cout << "Error setting the agregation function approach" << endl;;
		exit (EXIT_FAILURE);
//	}
//	return this->getScalarValue();
}

double Subprob::ScalarFitnessValueNormalized(Lambda lambda, int scalar) {

	double scale[numbObj];
//	double minn[numbObj];
//	for(int k=0;k<numbObj;k++){ //get max and min from the current population
//		int min = intMax, max = intMin;
//		for(int j=0;j<subprob;j++){
//
//			int tp =  solution[j].getFunction(k);
//
//			if(tp > max) max = tp;
//			if(tp < min) min = tp;
//		}
//		int tp = this->getFunction(k);
//		if(tp > max) max = tp;
//		if(tp < min) min = tp;
//		scale[k]=(double)max - min;
//		minn[k]=(double)min;
//	}
	for(int k=0;k<numbObj;k++) scale[k] = (double) worstReferencePoint[k]-ReferencePoint[k];

	if(scalar==0){ //weighted sum normalizded ------------------------------------
		this->setScalarValue(0.0);
		for(int k=0;k<numbObj;k++){
			double s = lambda.getWeight(k)*(((double)this->getFunction(k)-ReferencePoint[k])/scale[k]);

			this->setScalarValue(this->getScalarValue() + s);
		}

	}
	else if(scalar==2){ //weight sum alternative same as MOLSD
		this->setScalarValue(0.0);
		for(int k=0;k<numbObj;k++){
//			int factor; if(k==0) factor=9; else factor=1;
			double s = lambda.getWeight(k)*(double)((double)(this->getFunction(k)/(double)ReferencePoint[k]));
			this->setScalarValue(this->getScalarValue()+ s);
		}
	}
	else if(scalar==3){ //weighted sum normalizded and reference point * constant------------------------------------
		this->setScalarValue(0.0);
		for(int k=0;k<numbObj;k++){
			double s = lambda.getWeight(k)*(fabs((double)this->getFunction(k)-(0.6*(double)ReferencePoint[k]))/scale[k]);

			this->setScalarValue(this->getScalarValue() + s);
		}

	}
	else if(scalar==4){ //thebycheff
		this->setScalarValue(intMin); //techycheff
		for(int k=0;k<numbObj;k++){
			double dif = (double)((double)this->getFunction(k) - 0.6*(double)ReferencePoint[k]);
			dif = fabs(dif);
			double s = (double)lambda.getWeight(k)*dif;
			if(s > this->getScalarValue()) this->setScalarValue(s);
		}
	}
	else if(scalar==5){ //normalized tchebbycheff same as MOEADD

		this->setScalarValue((double)intMin);
		for(int k=0;k<numbObj;k++){
			double dif = fabs((double)((double)this->getFunction(k)-(0.8*(double)ReferencePoint[k]))) / (double)scale[k]; //1.1*ReferencePoint[k] - this->getFunction(k);
			double s;

			 s = dif * lambda.getWeight(k);
			if(s > (double)this->getScalarValue()) this->setScalarValue(s);
		}

	}
	else if(scalar==6){ //PBI
		double theta=5;
		vector <double> nambda(numbObj);
		for(int k=0;k<numbObj;k++){
			nambda[k]=lambda.getWeight(k);
		}
		//normalize the weight vector (line segment)
		double nd = Norma(nambda);
		for(int k=0;k<numbObj;k++){
			nambda[k] = lambda.getWeight(k)/nd;
		}
		vector <double> realA(numbObj);
		vector <double> realB(numbObj);
		//difference between reference point and current point (maximization)
		for(int k=0;k<numbObj;k++){
			realA[k] = (0.6*(double)ReferencePoint[k] - this->getFunction(k))/ (double)scale[k];
		}
		//distance along the line segment
		double d1 = fabs(innerProduct(realA, nambda));
		//distance to the line segment
		for(int k=0;k<numbObj;k++){
			realB[k]=(this->getFunction(k) - (0.6*(double)ReferencePoint[k] - d1*nambda[k])) / (double)scale[k];
		}
		double d2 = Norma(realB);
		double pbi = (double)(d1 + theta*d2);
		this->setScalarValue(pbi);

	}
	else{
		cout << "Error setting the agregation function approach" << endl;;
		exit (EXIT_FAILURE);
	}
	return this->getScalarValue();
}

void Subprob::InitializeSolution(int problem, Lambda lambda) {
	if(problem==0 or problem==1 or problem==3){ //mUBQP MOKP Trap5
		for(int i=0;i<numbVar;i++){
			double r = (double)rand()/(double) RAND_MAX;
			if(r < 0.5){
				this->setVectorSolution(1,i);
			}else{
				this->setVectorSolution(0,i);
			}
		}
	}else if(problem ==2){ //MTSP
		if(initialization==0){
			for(int i=0;i<numbVar;i++){
				this->setNoVisited(true,i);
			}
			int currentCity = rand()%numbVar; //randomly choose the first city to start the tour
			int nextCity;
			this->setNoVisited(false,currentCity);
			for(int i=1;i<=numbVar;i++){ //while there are cities available to be add
				double minDist =0.0;
				double dist[varMax];
				for(int j=0;j<numbVar;j++){
					dist[j]=0.0;
					if(this->getNoVisited(j)==true){ //verify if the city j is not visited yet, if not, so, calcule the probability to choose it
						for(int k=0;k<numbObj;k++){
							dist[j]=dist[j]+(lambda.getWeight(k) * matrixOfEdges[k][currentCity][j]); //calculation of the heuristic value
						}
						dist[j] = (float)(1/dist[j]);
						if(dist[j]>minDist){
							minDist=dist[j];
							nextCity = j;
						}
					}
				}
				this->setVectorSolution(nextCity,currentCity); //add the next city in the index of the current city
				this->setNoVisited(false,nextCity); //set the index of the next city as visited
				currentCity = nextCity;
			}
//			this->FitnessEvaluation(problem); // function f(x) of each new solution initialized
		}else if(initialization==1){
			vector<int> vectorOfVar;
			for(int i=0;i<numbVar;i++){
				vectorOfVar.push_back(i);
			}
			random_shuffle(vectorOfVar.begin(),vectorOfVar.end());
			for(int i=0;i<numbVar;i++){
				int currentCity = vectorOfVar.at(i);
				this->setVectorSolution(currentCity,i);
			}
		}else{
			cout << "Error setting the the type of initialization for permutation" << endl;;
			exit (EXIT_FAILURE);
		}
		//fitness evaluation
//		this->FitnessEvaluation(problem);
	}else if(problem==4 or problem==5){
		vector<int> initialPermutation;
		for(int i=0;i<numbVar;i++){
			initialPermutation.push_back(i);
		}
		random_shuffle(initialPermutation.begin(),initialPermutation.end());
		for(int i=0;i<numbVar;i++){
			this->setVectorSolution(initialPermutation.at(i),i);
			cout << initialPermutation.at(i) << " ";
		}
		cout << endl;
	}
	else{
		cout << "Error setting the the initialization" << endl;;
		exit (EXIT_FAILURE);
	}

}

int* Subprob::getPointerToSolution() {
  return vectorSolution;           // Pointer to the solution
}

int* Subprob::getPointerToFitness() {
  return function;           // Pointer to the solution
}

//void Subprob::setPointerToSolution(int* vector) {
//	vectorSolution=vector;
//}


Subprob::Subprob() {
	initialization = 0; //0 means problem specific and 1 means random
}

Subprob::~Subprob() {
}

void Subprob::PrintSolution() {
	for(int i=0;i<numbVar;i++){
		cout << this->getVectorSolution(i) << " ";
	}
	cout << (double)(this->getFunction(0))<< " " << (double)(this->getFunction(1)) << endl;
}

void Subprob::SetSoluion(int * vector){
	for(int i=0;i<numbVar;i++){
		this->setVectorSolution(vector[i],i);
	}
}

void Subprob::SetFitness(int * fitness){
	for(int i=0;i<numbObj;i++){
		this->setFunction(fitness[i],i);
	}
}

void Subprob::SetScalar(double sc){
	this->scalarValue=sc;
}
