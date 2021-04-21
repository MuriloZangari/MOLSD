/*
 * Distance.h
 *
 *  Created on: 03/03/2015
 *      Author: hydra
 */

#ifndef DISTANCE_H_
#define DISTANCE_H_

using namespace std;

class Distance{

private:
	double dist;
	int indiceViz;

public:
	 	void setDist(double dist);
		double getDist(void);
	 	void setIndiceViz( int ind );
		int getIndiceViz( void );


	Distance();
	~Distance();

};


#endif /* DISTANCE_H_ */
