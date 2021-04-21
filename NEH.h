/*
 *  LR.h
 *  AGA
 *
 *  Created by Josu Ceberio Uribe on 1/9/12.
 *  Copyright 2012 University of the Basque Country. All rights reserved.
 *
 */

#include "Problem.h"
#include "Global.h"
#include "NTools.h"

using namespace std;

class NEH{

private:


public:

	double * processingTime;
	int * gen;

	int * JIBIS(int * permutation, int current_fit, AbstractProblem * m_fsp, int IndexFnc);

	int * OSSBIS(int * permutation, int current_fit, AbstractProblem * m_fsp, int indexFnc);

	int * JIBSS(int * permutation, int current_fit, AbstractProblem * m_fsp, int indexFnc);

	int * Improvements(int * permutation, int current_fit, AbstractProblem * m_fsp, int indexFnc);

	void sortProcessingTimesCmax(AbstractProblem * m_fsp);

	void sortProcessingTimesCmax_SKE(AbstractProblem * m_fsp);

	int * nehHeuristicCmax(AbstractProblem * m_fsp);

	int * sortProcessingTimesTWT_EDDP(AbstractProblem * m_fsp);

	int * sortProcessingTimesTWT_EDD(AbstractProblem * m_fsp);

	int * sortProcessingTimesTWT_LPT(AbstractProblem * m_fsp);

	int * nehHeuristicTWT(AbstractProblem * m_fsp);

	void sortAsc();

	void sortDsc();



	NEH();
	~NEH();

};




