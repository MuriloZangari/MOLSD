#include"Distance.h"
#include <iostream>



Distance::Distance()
{
	//cout << "Object is being created" << endl;
}
Distance::~Distance(void)
{
	//cout << "Object is being deleted" << endl;
}

void Distance::setDist(double distt){
	dist = distt;
}
double Distance::getDist(void){
	return dist;
}
void Distance::setIndiceViz(int ind){
	indiceViz = ind;
}
int Distance::getIndiceViz(void){
	return indiceViz;
}
