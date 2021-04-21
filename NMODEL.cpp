
#include "NMODEL.h"




AbstractModel::AbstractModel() {


//	probVec = new double[subProbMax];
//	numbSubProb = nsp;
//    for(int i=0;i<numbVar;i++){
//    	for(int j=0;j<numbVar;j++){
//    		this->FreqMatrix[i][j]=0;
//    	}
//    }
}

AbstractModel::~AbstractModel() {
	delete[] probVec;
}

//void AbstractModel::ProbVectorScFun(Lambda lambda, Subprob* solution, double* probVec) {
//
//	double min=(double)intMax, max=(double)intMin;
//	for(int i=0;i<lambda.NeigbhborSizeSimga;i++){ //step1: compute the scalar function for each neigbhor and put on probVec
//		int p1 = lambda.getIndexOfNeighbors(i);
//		probVec[i]= (double) solution[p1].getScalarValue();
//
//		if(probVec[i]<min){
//			min = probVec[i];
//		}
//		if(probVec[i]>max){
//			max = probVec[i];
//		}
//	}
//	double sumExp=0.0;
//	for(int i=0;i<lambda.NeigbhborSizeSimga;i++){//step 2: put the value between 0 and 1 and exponential
//		probVec[i]= (probVec[i]-min) / (max-min);
//		probVec[i]=exp(probVec[i]);
//		sumExp = sumExp + probVec[i];
//	}
//	for(int i=0;i<lambda.NeigbhborSizeSimga;i++){
//		probVec[i]=(double)(probVec[i]/sumExp);
//	}
//
//}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
GA::GA(int nsp, float mut):AbstractModel() {
	sigma = 1.0;
	mutation = mut;
}

GA::~GA(){
}
void GA::Learning(Lambda lambda, int numSamp){

		int id1 = (int) rand() % lambda.NeigbhborSizeSimga;
		int id2 = (int) rand() % lambda.NeigbhborSizeSimga;
//		cout << "list of ngbr " ; for(int j=0;j<lambda.NeigbhborSizeSimga;j++) cout << lambda.getIndexOfNeighbors(j) << " ";
//		cout << endl;
//		cout << "choosen " << id1 << " " << id2 << endl;
		p1 = lambda.getIndexOfNeighbors(id1);
		p2 = lambda.getIndexOfNeighbors(id2);
		index = lambda.getIndexOfNeighbors(0);
		//cout << "sub " << lambda.getIndexOfNeighbors(0) << " pais: " << p1[s] << " " << p2[s] << endl;

//	this->ShowModel(numSamp);
//	this->PrintSelectPop(lambda, solution);
}

void GA::Sampling(Subprob* solution, int s, Subprob &newSol) {

	double rd = (double)rand()/(double)RAND_MAX;
	int half_orbit=5;
	int min_range,max_range,val;
	int son[numbVar];
	if(rd <= sigma){

	//		int daughter[numbVar];
		int mother[numbVar];
		int father[numbVar];
	//		cout << "f "; solution[this->p1[s]].PrintSolution();
	//		cout << "m "; solution[this->p2[s]].PrintSolution();
		for(int i=0;i<numbVar;i++) mother[i] = solution[p1].getVectorSolution(i);
		for(int i=0;i<numbVar;i++) father[i] = solution[p2].getVectorSolution(i);
		int aux;
		//PrintArray(father,IND_SIZE,"in father TTP: ");
		//PrintArray(mother,IND_SIZE,"in mother TTP: ");

		//Auxiliary structures.
		int mother_aux[numbVar];
	//		int father_aux[numbVar];
		for (int i=0;i<numbVar;i++)
		{
			mother_aux[i]=1;
	//			father_aux[i]=1;
		}

		//first crossover
		int pa1 = rand() % (numbVar);
		int pa2 = rand() % (numbVar);
		if (pa1>pa2)
		{
			aux=pa1;
			pa1=pa2;
			pa2=aux;
		}

		//cout<<"Crossover points: "<<p1<<" "<<p2<<endl;
		//block from father
		for (int i=pa1;i<=pa2;i++){
			son[i]=father[i];
			mother_aux[father[i]]=0;
		}
		//PrintArray(mother_aux , IND_SIZE, "Aux data: ");
		//remaining from mother
		//first section
		//cout<<"first section"<<endl;
		int z=0;
		for (int j=0;j<pa1;j++)
		{
			if (mother_aux[mother[z]]==1)
			{
				son[j]=mother[z];
				mother_aux[mother[z]]=0;
			}
			else
				j--;
			z++;
		}

		//PrintArray(son,IND_SIZE,"middle: ");
		//PrintArray(mother_aux , IND_SIZE, "Aux data in the middle: ");
		//cout<<"second section"<<endl;
		//second section
		for (int j=pa2+1;j<numbVar;j++)
		{
		//	cout<<"next mother value "<<mother[z];
			if (mother_aux[mother[z]]==1)
			{
		//		cout<<" is free"<<endl;
				son[j]=mother[z];
				mother_aux[mother[z]]=0;
		//		PrintArray(son,IND_SIZE,"ins middle: ");
		//		PrintArray(mother_aux , IND_SIZE, "ins Aux data in the middle: ");
			}
			else{
		//		cout<<" is NOT free"<<endl;
				j--;}
			z++;
		}
		//PrintArray(son,IND_SIZE,"out son TTP: ");
		//exit(1);

	//		//second crossover
	//		int p3 = rand() % (numbVar);
	//	    int p4 = rand() % (numbVar);
	//		if (p3>p4)
	//		{
	//			aux=p3;
	//			p3=p4;
	//			p4=aux;
	//		}
	//		//cout<<"Crossover points: "<<p3<<" "<<p4<<endl;
	//		//block from mother
	//		for (int i=p3;i<=p4;i++){
	//			daughter[i]=mother[i];
	//			father_aux[mother[i]]=0;
	//		}
	//		//remaining from father
	//		//first section
	//		z=0;
	//		for (int j=0;j<p3;j++)
	//		{
	//			if (father_aux[father[z]]==1)
	//			{
	//				daughter[j]=father[z];
	//				father_aux[father[z]]=0;
	//			}
	//			else
	//				j--;
	//			z++;
	//		}
	//		//second section
	//		for (int j=p4+1;j<numbVar;j++)
	//		{
	//			if (father_aux[father[z]]==1)
	//			{
	//				daughter[j]=father[z];
	//				father_aux[father[z]]=0;
	//			}
	//			else
	//				j--;
	//			z++;
	//		}

	}else{
		for(int i=0;i<numbVar;i++) son[i] = solution[index].getVectorSolution(i);
	}

	//PrintArray(daughter,IND_SIZE,"out daughter TTP: ");
//		cout << "s ";
//		for(int i=0;i<numbVar;i++) cout << son[i] << " ";
//		cout << endl << "d ";
//		for(int i=0;i<numbVar;i++) cout << daughter[i] << " ";
//		cout << endl;


	// INSERT MUTATION
	rd = (double)rand()/(double)RAND_MAX;
	if(rd<=0.5){
		int ii = rand() % (numbVar);
		//int jj = rand() % (numbVar);
        min_range=MAX(ii-half_orbit,0);
        max_range=MIN(ii+half_orbit,numbVar-1);
        val = rand() % (max_range-min_range);
        int jj= min_range+val;
		InsertAt(son,ii,jj,numbVar);
	}
	for(int i=0;i<numbVar;i++) newSol.setVectorSolution(son[i],i);
//			cout << "choosed son" << endl;

//		cout << "son: "; newSol.PrintSolution();


}
void GA::InitializeModel(){

}

void GA::ShowModel(int n) {

	cout << "parents from " << index << ": " << p1 << " " << p2 << endl;

}

//////////////////////////////GENETIC LOCAL SEARCH ////////////////////////////////////////////////

GLS::GLS(float probC, float mut, float pX):AbstractModel() {
	sigma = probC;
	mutation = mut;
	probX = pX;

//	cout << sigma << " " << mutation << " " << probX << endl;
}

GLS::~GLS(){
}
void GLS::Learning(Lambda lambda, int numSamp){

		int id1 = (int) rand() % lambda.NeigbhborSizeSimga;
		int id2;
		do{
			 id2 = (int) rand() % lambda.NeigbhborSizeSimga;
		}while(id2!=id1);

		p1 = lambda.getIndexOfNeighbors(id1);
		p2 = lambda.getIndexOfNeighbors(id2);
		index = lambda.getIndexOfNeighbors(0);
		//cout << "sub " << lambda.getIndexOfNeighbors(0) << " pais: " << p1[s] << " " << p2[s] << endl;

//	this->ShowModel(numSamp);
//	this->PrintSelectPop(lambda, solution);
}

void GLS::Sampling(Subprob* solution, int s, Subprob &newSol) {

	double rd = (double)rand()/(double)RAND_MAX;
	//int half_orbit=5;
	//int min_range,max_range,val;
	int son[numbVar];

	if(rd <= (double)sigma){

	//		int daughter[numbVar];
		int mother[numbVar];
		int father[numbVar];
	//		cout << "f "; solution[this->p1[s]].PrintSolution();
	//		cout << "m "; solution[this->p2[s]].PrintSolution();
		for(int i=0;i<numbVar;i++) mother[i] = solution[p1].getVectorSolution(i);
		for(int i=0;i<numbVar;i++) father[i] = solution[p2].getVectorSolution(i);
		int aux;
		//PrintArray(father,IND_SIZE,"in father TTP: ");
		//PrintArray(mother,IND_SIZE,"in mother TTP: ");

		//Auxiliary structures.
		int mother_aux[numbVar];
	//		int father_aux[numbVar];
		for (int i=0;i<numbVar;i++)
		{
			mother_aux[i]=1;
	//			father_aux[i]=1;
		}

		//first crossover
		int pa1 = rand() % (numbVar);
		int pa2 = rand() % (numbVar);
		if (pa1>pa2)
		{
			aux=pa1;
			pa1=pa2;
			pa2=aux;
		}

		//cout<<"Crossover points: "<<p1<<" "<<p2<<endl;
		//block from father
		for (int i=pa1;i<=pa2;i++){
			son[i]=father[i];
			mother_aux[father[i]]=0;
		}
		//PrintArray(mother_aux , IND_SIZE, "Aux data: ");
		//remaining from mother
		//first section
		//cout<<"first section"<<endl;
		int z=0;
		for (int j=0;j<pa1;j++)
		{
			if (mother_aux[mother[z]]==1)
			{
				son[j]=mother[z];
				mother_aux[mother[z]]=0;
			}
			else
				j--;
			z++;
		}

		//PrintArray(son,IND_SIZE,"middle: ");
		//PrintArray(mother_aux , IND_SIZE, "Aux data in the middle: ");
		//cout<<"second section"<<endl;
		//second section
		for (int j=pa2+1;j<numbVar;j++)
		{
		//	cout<<"next mother value "<<mother[z];
			if (mother_aux[mother[z]]==1)
			{
		//		cout<<" is free"<<endl;
				son[j]=mother[z];
				mother_aux[mother[z]]=0;
		//		PrintArray(son,IND_SIZE,"ins middle: ");
		//		PrintArray(mother_aux , IND_SIZE, "ins Aux data in the middle: ");
			}
			else{
		//		cout<<" is NOT free"<<endl;
				j--;}
			z++;
		}

	}else{
			for(int i=0;i<numbVar;i++) son[i] = solution[index].getVectorSolution(i);
	}

	//// mutation

	rd = (double)rand()/(double)RAND_MAX;
	if(rd<=(double)mutation){ //if mutation
		int i,j;
		double rd2 = (double)rand()/(double)RAND_MAX;

		if(rd2<=(double)probX){

			i = rand() % numbVar;
			j = rand() % numbVar;
			Swap(son,i,j);
		}else{

			i = rand() % numbVar;
			j = rand() % numbVar;
			InsertAt(son,i,j,numbVar);
		}
	}

	newSol.SetSoluion(son);

}
void GLS::InitializeModel(){

}

void GLS::ShowModel(int n) {

	cout << "GLS index " << index << endl;

}


////////////////////////////LOCAL SEARCH /////////////////////////////////////////////////////////////

LS::LS(int ins):AbstractModel() {
probMut = 0.0;
}

LS::~LS(){
}

void LS::Learning(Lambda lambda, int numSamp){
		index = lambda.getIndexOfNeighbors(0);
		//ShowModel(index);

}


void LS::Sampling(Subprob* solution, int s, Subprob &newSol) {
    int * inverted= new int[numbVar];

    memcpy(inverted,solution[index].getPointerToSolution(),sizeof(int)*numbVar);
    int i,j;

   //int numSwap;
   if((string)pho=="c"){
		i = rand() % numbVar;
		j = rand() % numbVar;
		Swap(inverted,i,j);

    }else if((string)pho=="u"){
		i = rand() % numbVar;
		j = rand() % numbVar;

		InsertAt(inverted,i,j,numbVar);

	}else {cout << "error seting distance metric " << endl; exit (EXIT_FAILURE);}

    newSol.SetSoluion(inverted);

    delete[] inverted;


}


void LS::ShowModel(int n) {

//	cout << "LS index " << n << endl;

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//PBILModel::PBILModel(int nsp, int mut, float lr):AbstractModel(nsp) {
//	alpha=lr;
//	mutation=mut;
//	sigma = 1;
//
//}
//
//PBILModel::~PBILModel(){
//}
//
//void PBILModel::Learning(Lambda lambda, Subprob* solution, int numSamp) {
////	this->ProbVectorScFun(lambda, solution, this->probVec);
//	double rd = (double)rand()/(double)RAND_MAX;
//
//	if(rd<= sigma){
//		for(int k=0;k<lambda.NeigbhborSizeSimga;k++){
//			int p1 = lambda.getIndexOfNeighbors(k);
//			for(int i=0;i<numbVar;i++){
//				this->prob[i] = (double) (this->prob[i]* (1.0 - alpha)) + (double) (alpha*solution[p1].getVectorSolution(i));
//			}
//		}
//	}else{
//		for(int k=0;k<lambda.NeigbhborSizeSimga;k++){
//			int p1 = (int) rand() % numbSubProb;
//			for(int i=0;i<numbVar;i++){
//				this->prob[i] = (double) (this->prob[i] * (1.0 - alpha)) + (double) (alpha*solution[p1].getVectorSolution(i));
//			}
//		}
//	}
////this->ShowModel();
//
//}
//
//void PBILModel::Sampling(Subprob* solution, int s, Subprob &newSol) {
//	for(int i=0;i<numbVar;i++){
//		double r = (double)rand()/(double)RAND_MAX;
//		if(r < prob[i]){
//			newSol.setVectorSolution(1,i);
//		}else{
//			newSol.setVectorSolution(0,i);
//		}
//		if(mutation==1){
//			r = (double)rand()/(double)RAND_MAX;
//			if(r < (double)(1.0/(double)numbVar)){ //mutation
//				//cout << (double)(1.0/(double)numbVar) << endl;
//				newSol.setVectorSolution((int)1-newSol.getVectorSolution(i),i);
//			}
//		}
//	}
////	cout << "new sampled solution: ";
////	for(int j=0;j<numbVar;j++){
////		cout << newSol.getVectorSolution(j) << " ";
////	}
////	cout << endl;
//
//}
//
//void PBILModel::InitializeModel(){
//	this->InitializeProb();
//}
//
//void PBILModel::InitializeProb(){
//	for(int i=0;i<numbVar;i++){
//		this->prob[i]=0.5;
//	}
//}
//
//void PBILModel::ShowModel(){
//	cout << "PBIL " << endl;
//	for(int i=0;i<numbVar;i++){
//		cout << this->prob[i] << " ";
//	}
//}
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//PBILOrderInvers::PBILOrderInvers(int nsp, int mut, float lr):AbstractModel(nsp) {
//	alpha=lr;
//	mutation=mut;
//	sigma = 1;
//}
//
//PBILOrderInvers::~PBILOrderInvers(){
//}
//
//void PBILOrderInvers::Learning(Lambda lambda, Subprob* solution, int numSamp) {
//	double rd = (double)rand()/(double)RAND_MAX;
//	if(rd<= sigma){
//		double sumScalarValue = 0.0;
//		Distance aggregationFunctionOrder[lambda.NeigbhborSizeSimga];
//		for(int i=0;i<lambda.NeigbhborSizeSimga;i++){
//			int p1 = lambda.getIndexOfNeighbors(i);
//			aggregationFunctionOrder[i].setIndiceViz(p1);
//			double f = solution[p1].ScalarFitnessValue(solution[p1],lambda,0);
//			aggregationFunctionOrder[i].setDist(f);
//			sumScalarValue = sumScalarValue + f;
//		}
//		//shell sort
//		int i, flag = 1;
//		Distance temp;
//		int d = lambda.NeigbhborSizeSimga;
//		while( flag || (d > 1)){  // boolean flag (true when not equal to 0)
//		flag = 0; // reset flag to 0 to check for future swaps
//		d = (d+1) / 2;
//			for (i = 0; i < (lambda.NeigbhborSizeSimga - d); i++){
//				if (aggregationFunctionOrder[i+d].getDist() > aggregationFunctionOrder[i].getDist()){ //change the sinal
//				  temp = aggregationFunctionOrder[i+d]; // swap positions i+d and i
//				  aggregationFunctionOrder[i+d]=aggregationFunctionOrder[i];
//				  aggregationFunctionOrder[i]=temp;
//				  flag = 1; // tells swap has occurred
//				}
//			}
//		}
//		for(int k=0;k<lambda.NeigbhborSizeSimga;k++){
//			int p1 = aggregationFunctionOrder[k].getIndiceViz();
//				double w = solution[p1].getScalarValue()/sumScalarValue;
//				double lr = alpha + (double)(w/(double)lambda.NeigbhborSizeSimga);
//			//	cout << lambda.getIndexOfNeighbors(0) << " " << k << " " <<  w << " " << lr << endl;
//
//				for(int i=0;i<numbVar;i++){
//					this->prob[i] = (double) (this->prob[i] * (1.0 - lr)) + (double) (lr*solution[p1].getVectorSolution(i));
//				}
//		}
//	}else{
//		for(int k=0;k<lambda.NeigbhborSizeSimga;k++){
//			int p1 = (int) rand() % numbSubProb;
//			for(int i=0;i<numbVar;i++){
//				this->prob[i] = (double) (this->prob[i] * (1.0 - alpha)) + (double) (alpha*solution[p1].getVectorSolution(i));
//			}
//		}
//	}
//}
//
//void PBILOrderInvers::Sampling(Subprob* solution, int s, Subprob &newSol) {
//	for(int i=0;i<numbVar;i++){
//		double r = (double)rand()/(double)RAND_MAX;
//		if(r < prob[i]){
//			newSol.setVectorSolution(1,i);
//		}else{
//			newSol.setVectorSolution(0,i);
//		}
//		if(mutation==1){
//			r = (double)rand()/(double)RAND_MAX;
//			if(r < (double)(1.0/(double)numbVar)){ //mutation
//				//cout << "entrou mutation" << endl;
//				newSol.setVectorSolution((int)1-newSol.getVectorSolution(i),i);
//			}
//		}
//	}
//}
//
//void PBILOrderInvers::InitializeModel(){
//	this->InitializeProb();
//}
//
//void PBILOrderInvers::InitializeProb(){
//	for(int i=0;i<numbVar;i++){
//		this->prob[i]=0.5;
//	}
//}
//
//void PBILOrderInvers::ShowModel(){
//	cout << endl;
//	for(int i=0;i<numbVar;i++){
//		cout << this->prob[i] << " ";
//	}
//}
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//PBILOrder::PBILOrder(int nsp, int mut, float lr):AbstractModel(nsp) {
//	alpha=lr;
//	mutation=mut;
//	sigma = 1;
//}
//
//PBILOrder::~PBILOrder(){
//}
//
//void PBILOrder::Learning(Lambda lambda, Subprob* solution, int numSamp) {
//	double rd = (double)rand()/(double)RAND_MAX;
//	if(rd<= sigma){
//		double sumScalarValue = 0.0;
//		Distance aggregationFunctionOrder[lambda.NeigbhborSizeSimga];
//		for(int i=0;i<lambda.NeigbhborSizeSimga;i++){
//			int p1 = lambda.getIndexOfNeighbors(i);
//			aggregationFunctionOrder[i].setIndiceViz(p1);
//			double f = solution[p1].ScalarFitnessValue(solution[p1],lambda,0);
//			aggregationFunctionOrder[i].setDist(f);
//			sumScalarValue = sumScalarValue + f;
//		}
//		//shell sort
//		int i, flag = 1;
//		Distance temp;
//		int d = lambda.NeigbhborSizeSimga;
//		while( flag || (d > 1)){  // boolean flag (true when not equal to 0)
//		flag = 0; // reset flag to 0 to check for future swaps
//		d = (d+1) / 2;
//			for (i = 0; i < (lambda.NeigbhborSizeSimga - d); i++){
//				if (aggregationFunctionOrder[i+d].getDist() < aggregationFunctionOrder[i].getDist()){ //change the sinal
//				  temp = aggregationFunctionOrder[i+d]; // swap positions i+d and i
//				  aggregationFunctionOrder[i+d]=aggregationFunctionOrder[i];
//				  aggregationFunctionOrder[i]=temp;
//				  flag = 1; // tells swap has occurred
//				}
//			}
//		}
//		for(int k=0;k<lambda.NeigbhborSizeSimga;k++){
//			int p1 = aggregationFunctionOrder[k].getIndiceViz();
//				double w = solution[p1].getScalarValue()/sumScalarValue;
//				double lr = alpha + (double)(w/(double)lambda.NeigbhborSizeSimga);
//			//	cout << lambda.getIndexOfNeighbors(0) << " " << k << " " <<  w << " " << lr << endl;
//
//				for(int i=0;i<numbVar;i++){
//					this->prob[i] = (double) (this->prob[i] * (1.0 - lr)) + (double) (lr*solution[p1].getVectorSolution(i));
//				}
//		}
//	}else{
//		for(int k=0;k<lambda.NeigbhborSizeSimga;k++){
//			int p1 = (int) rand() % numbSubProb;
//			for(int i=0;i<numbVar;i++){
//				this->prob[i] = (double) (this->prob[i] * (1.0 - alpha)) + (double) (alpha*solution[p1].getVectorSolution(i));
//			}
//		}
//	}
//}
//
//void PBILOrder::Sampling(Subprob* solution, int s, Subprob &newSol) {
//	for(int i=0;i<numbVar;i++){
//		double r = (double)rand()/(double)RAND_MAX;
//		if(r < prob[i]){
//			newSol.setVectorSolution(1,i);
//		}else{
//			newSol.setVectorSolution(0,i);
//		}
//		if(mutation==1){
//			r = (double)rand()/(double)RAND_MAX;
//			if(r < (double)(1.0/(double)numbVar)){ //mutation
//				//cout << "entrou mutation" << endl;
//				newSol.setVectorSolution((int)1-newSol.getVectorSolution(i),i);
//			}
//		}
//	}
//}
//
//void PBILOrder::InitializeModel(){
//	this->InitializeProb();
//}
//
//void PBILOrder::InitializeProb(){
//	for(int i=0;i<numbVar;i++){
//		this->prob[i]=0.5;
//	}
//}
//
//void PBILOrder::ShowModel(){
//	cout << endl;
//	for(int i=0;i<numbVar;i++){
//		cout << this->prob[i] << " ";
//	}
//}
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//IUMDA::IUMDA(int nsp, int mut, float factor, bool boltz ):AbstractModel(nsp) {
//	factoR=factor;
//	mutation = mut;
//	bolztman=boltz;
//
//}
//
//IUMDA::~IUMDA(){
//}
//
//void IUMDA::Learning(Lambda lambda, Subprob* solution, int numSamp) {
//	if(bolztman){
//		this->ProbVectorScFun(lambda, solution, this->probVec);
//
//	}else{
//		for(int i=0;i<numbSubProb;i++){
//			this->probVec[i]=(double)((double)1/(double)lambda.NeigbhborSizeSimga);
////			cout << probVec[i] << endl;
//		}
//	}
//
//	for(int i=0;i<numbVar;i++){
//		double frequency = 0.0;
//		for(int k=0;k<lambda.NeigbhborSizeSimga;k++){
//			int p1 = lambda.getIndexOfNeighbors(k);
////			cout << probVec[k] << endl;
//			frequency = frequency + (double)(solution[p1].getVectorSolution(i)*probVec[k]);
//
//		}
//
////		frequency = (double)frequency/lambda.NeigbhborSizeSimga;
////		cout << frequency << endl;
//		this->prob[i] = (double) this->prob[i]*(1.0-factoR) + (double)factoR*frequency;
//	}
////	ShowModel();
//
//}
//
//void IUMDA::Sampling(Subprob* solution, int s, Subprob &newSol) {
//	for(int i=0;i<numbVar;i++){
//		double r = (double)rand()/(double)RAND_MAX;
//		if(r < prob[i]){
//			newSol.setVectorSolution(1,i);
//		}else{
//			newSol.setVectorSolution(0,i);
//		}
//		if(mutation==1){
//			r = (double)rand()/(double)RAND_MAX;
//			if(r < (double)(1.0/(double)numbVar)){ //mutation
//				//cout << "entrou mutation" << endl;
//				newSol.setVectorSolution((int)1-newSol.getVectorSolution(i),i);
//			}
//		}
//	}
//}
//
//void IUMDA::InitializeModel(){
//	this->InitializeProb();
//}
//
//void IUMDA::InitializeProb(){
//	for(int i=0;i<numbVar;i++){
//		this->prob[i]=0.5;
//	}
//}
//
//void IUMDA::ShowModel(){
//	cout << endl;
//	for(int i=0;i<numbVar;i++){
//		cout << this->prob[i] << " ";
//	}
//}
//
//
///////////////////////////////////////////////////////////////////////////
///*
//   Tree Class. It is implemented using multiple-inheritance.
//   The class TREE,  will inherit from AbstractModel and from IntTreeModel
//
//*/
//
//
//TREE::TREE(int nsp, int ttree, unsigned int* cards, int mut, double prior, bool boltz):AbstractModel(nsp),IntTreeModel(numbVar,0.75,50,cards){
//	Prior = prior;
//	bolztman=boltz;
//	sigma = 1;
//    typetree = ttree;
//    mutation = mut;
//
//    //HERE INITIALIZE THE AllProb
//
//}
//
//
//TREE::~TREE(){
//}
//void TREE::ComputeFrequencies(){
//	for(int j=0; j<length;j++){
//		if(Tree[Queue[j]]!=-1){
//			this->FreqMatrix[Queue[j]][Tree[Queue[j]]] = this->FreqMatrix[Queue[j]][Tree[Queue[j]]] + 1;
//			this->FreqMatrix[Tree[Queue[j]]][Queue[j]] = this->FreqMatrix[Tree[Queue[j]]][Queue[j]] + 1;
//		}
//	}
//}
//void TREE::Learning(Lambda lambda, Subprob* solution, int numSamp){
//
//     rootnode = RandomRootNode();                                   // A node in the tree is randomly selected
//
// 	if(bolztman){
//
// 		this->ProbVectorScFun(lambda, solution, this->probVec);
// 	}else{
// 		for(int i=0;i<numbSubProb;i++){
// 			this->probVec[i]=(double)((double)1/(double)lambda.NeigbhborSizeSimga);
// //			cout << probVec[i] << endl;
// 		}
// 	}
//
//     CalProbFvect(lambda, solution,lambda.NeigbhborSizeSimga, probVec);                 // The univariate and bivariate probabilities are learned
//                                                                  // from the solutions in the neighborhood
//     if(typetree==2)                                               // Disconnected tree (forest)
//       {
//        CalMutInf();                                              // Mutual information is computed
//        MakeTree(rootnode);                                       // The tree is constructed
//       }
//     else if(typetree==1)                                          // Chain-shaped tree
//       {
//        MakeTreeMarkovChain();
//       }
//     else if(typetree==0)                                          // Univariate model
//       MakeTreeIndependentVars();
//
////     ShowModel();
//     ComputeFrequencies();
//}
//
//void TREE::Sampling(Subprob* solution, int s, Subprob &newSol) {
//	GenIndividual(newSol, Prior, mutation);
//
//
//}
//void TREE::InitializeModel(){
////	cout << "neighborhood size " <<  neighborSizeSelection << endl;
//}
//
//void TREE::ShowModel() {
//	for (int i=0;i<length; i++ ) cout<<" "<<AllProb[IndexUnivEntries[i]]<<" "<<AllProb[IndexUnivEntries[i]+1];
//	cout<<endl;
//
//}
//void TREE::PrintFreqMatrix(int g){
//	string group = static_cast<ostringstream*>( &(ostringstream() << g) )->str();
//	ofstream freq;
//	freq.open(("FREQ/"+ test+"_"+group+".dat.txt").c_str());
//
//    for(int i=0;i<numbVar;i++){
//    	for(int j=0;j<numbVar;j++){
//    		freq << this->FreqMatrix[i][j] << " ";
//    	}
//    	freq << endl;
//    }
//    freq.close();
//}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/* 
   Mallows Class. It is implemented using multiple-inheritance. 
   The class MALLOWS,  will inherit from AbstractModel and from NMallowsModel

*/

// MALLOWS::MALLOWS(int nsp, char* metric_type):AbstractModel(),CMallowsModel(numbVar,metric_type){
//
//
//
//    float_probVec = new float[subProbMax];  // The Mallows code uses float precision
//                                            // so we will have to transform from double to float during learning
//
//}
//
//
//MALLOWS::~MALLOWS(){
//  delete[] float_probVec;
//}
//
//void MALLOWS::Learning(Lambda lambda, Subprob* solution, int numSamp, float theta){
//
////   	if(bolztman){
//// 		this->ProbVectorScFun(lambda, this->probVec);
//// 	}else{
//// 		for(int i=0;i<lambda.NeigbhborSizeSimga;i++)this->probVec[i]=(double)((double)1/(double)lambda.NeigbhborSizeSimga);
//// 	}
////
////    for(int i=0;i<lambda.NeigbhborSizeSimga;i++) float_probVec[i] = (float) probVec[i];
//
//	//Here we can pass different neigbhor size for the learning, instead of what is set on lambda.
//	LearnForMOEAD(lambda, solution, lambda.NeigbhborSizeSimga, theta); // This is the learning step of Mallows
//
//
//}
//
//void MALLOWS::Sampling(Subprob* solution, int s, Subprob &newSol) {
//  Sample(newSol.getPointerToSolution());
//
//}
//
//void MALLOWS::Sampling(Subprob &newSol) {
//  Sample(newSol.getPointerToSolution());
//
//}
//void MALLOWS::InitializeModel(){
//
//}
//
//void MALLOWS::ShowModel() {
//  // To be defined for Mallows
//
//}
//void MALLOWS::PrintFreqMatrix(int g){
//  // Probably there are not frequency matrices here
//}

void AbstractModel::PrintSelectPop(Lambda lambda, Subprob* solution) {
    for(int i=0;i<lambda.NeigbhborSizeSimga;i++){
    	int p1 = lambda.getIndexOfNeighbors(i);
    	cout << "i: " << p1 << " vec: ";
    	solution[p1].PrintSolution();
    }
}
