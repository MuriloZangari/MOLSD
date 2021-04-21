/*
 * Lambda.h
 *
 *  Created on: 03/03/2015
 *      Author: hydra
 */

#ifndef LAMBDA_H_
#define LAMBDA_H_

#include <iostream>
#include"Global.h"
#include"Distance.h"



using namespace std;

class Lambda{ //each weight vector is associate to a subproblem. NOTE: how to know it? the subproblem and the weight vector have the same INDEX

private:
	double weight[objMax]; //the weights for each objective
	int indexOfNeighbors[neighborMax]; //the index of neighbors for each subproblem. NOTE: there is no index, because the first neighbor is itself (because the shortest distance is zero), i.e., the index of the weight vector is the index of the first neighbor
	double distOfNeighbors[neighborMax];
public:
		int NeigbhborSizeSimga;
		int NeighborSizeUpdate;
		int numbSampling;
		//gets and sets of the weight and indexOfNeighbors
	 	void setWeight(double w, int i);
	 	float getWeight(int i);
		void setIndexOfNeighbors(int n, int i);
		int getIndexOfNeighbors(int i);
	 	void setDistOfNeighbors(double d, int i);
	 	double getDistOfNeighbors(int i);

	Lambda();
	~Lambda();

};


#endif /* LAMBDA_H_ */
