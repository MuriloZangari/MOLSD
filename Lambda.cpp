#include"Lambda.h"

Lambda::Lambda(){
	NeigbhborSizeSimga=0;

}
Lambda::~Lambda(void)
{
	//cout << "Object is being deleted" << endl;
}

void Lambda::setWeight(double w, int i){
	weight[i] = w;
}

float Lambda::getWeight(int i){
	return weight[i];
}

void Lambda::setIndexOfNeighbors(int n, int i){
	indexOfNeighbors[i]=n;
}
int Lambda::getIndexOfNeighbors(int i){
	return indexOfNeighbors[i];
}

void Lambda::setDistOfNeighbors(double d, int i) {
	distOfNeighbors[i] = d;
}

double Lambda::getDistOfNeighbors(int i) {
	return distOfNeighbors[i];
}
//params - the set of solutions to be sorted
//		 - the index of the objective in wich the solutions should be sorted
//		 - the number of solutions to be sorted


