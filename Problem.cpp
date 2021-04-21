#include "Problem.h"

AbstractProblem::AbstractProblem() {

}

AbstractProblem::~AbstractProblem() {
}

PFSP::PFSP(){
//	m_machines = machine;
//	m_jobs = numbVar;

}

PFSP::~PFSP(){
	for	(int i=0;i<m_machines;i++)
	delete [] m_processingtimes[i];
	delete [] m_processingtimes;
    delete [] m_timeTable;
    delete [] m_aux;
    delete [] m_DueDates;
}

void PFSP::Initialize() {

	//BUILD JOB PROCESSING MATRIX
	//cout << "--> BUILDING JOB PROCESSING MATRIX" << endl;


}

void PFSP::LoadInstance() {

	fstream dataIn;
	dataIn.open(("../FSSP/DD_Ta"+idxInstance+".txt").c_str(),ios::in);
	if(dataIn==NULL){
		cout << "Error opening file FSSP" << endl;
		exit(EXIT_FAILURE);
	}else{
		cout << "test instance found" << endl;
	}
	string line;
	int value;
	dataIn >> line; istringstream(line) >> value; numbVar = value; m_jobs = value; cout << numbVar << endl;
	dataIn >> line; istringstream(line) >> value; machine = value; m_machines = value; cout << machine << endl;

	m_processingtimes = new int*[m_machines];
	m_timeTable= new int[m_machines];
	m_aux= new int[m_jobs];
	for (int i=0;i<m_machines;i++)
	{
		m_processingtimes[i]= new int[m_jobs];
	}
	m_DueDates = new int[m_jobs];

	for(int i=0;i<m_jobs;i++){
		for(int j=0;j<m_machines;j++){
			dataIn >> line; istringstream(line) >> value; cout << value << " ";
			dataIn >> line; istringstream(line) >> value; cout << value << endl;
			m_processingtimes[j][i]=value;
		}
	}

	dataIn >> line;
	if(line=="Reldue"){
		cout << line << endl;
		for(int i=0;i<m_jobs;i++){
			dataIn >> line;
			dataIn >> line; istringstream(line) >> value; cout << value << endl; m_DueDates[i]=value;
			dataIn >> line; dataIn >> line;
		}
	}

	dataIn.close();

}

void PFSP::Show() {
	//print the flow shop


		for(int i=0;i<m_machines;i++){
			for(int j=0;j<m_jobs;j++){
				cout << m_processingtimes[i][j] << " ";
			}
			cout << endl;
		}
		cout << "Due Dates" << endl;
		for(int j=0;j<m_jobs;j++){
			cout << m_DueDates[j] << " ";
		}
		cout << endl;
}

int * PFSP::solutionFitness(int * gen){
	int * fitness = new int[numbObj];
	//objective 1 = Evaluates the given solution with the total flow time criterion.

		for (int i=0;i<m_machines;i++) m_timeTable[i]=0;
		int j,z, job;
		int machine;
		int prev_machine=0;

		int first_gene=gen[0];
		m_timeTable[0]=m_processingtimes[0][first_gene];
		for (j=1;j<m_machines;j++)
		{
			m_timeTable[j]=m_timeTable[j-1]+m_processingtimes[j][first_gene];
		}

		double fitness0=m_timeTable[m_machines-1];
		for (z=1;z<m_jobs;z++)
		{
			job=gen[z];

			//machine 0 is always incremental, so:
			m_timeTable[0]+=m_processingtimes[0][job];
			prev_machine=m_timeTable[0];
			for (machine=1;machine<m_machines;machine++)
			{
				m_timeTable[machine]= MAX(prev_machine,m_timeTable[machine])+ m_processingtimes[machine][job];
				prev_machine=m_timeTable[machine];
			}

			fitness0+=m_timeTable[m_machines-1];
		}


		fitness[1]=(int)fitness0;

	//objective 2 = Make Span

		double fitness1=0;
		double timeTable[m_machines];

		for (int j=0;j<m_machines;j++)
		{
			timeTable[j]=0;
		}
		for (int z=0;z<m_jobs;z++)
		{
			int job=gen[z];
			//cout<<"Job "<<job<<endl;
			for (int machine=0;machine<m_machines;machine++)
			{
				double processingTime=m_processingtimes[machine][job];
				if (machine==0)
				{
					timeTable[machine]+=processingTime;
				}
				else
				{
					if (timeTable[machine-1]<timeTable[machine])
					{
						timeTable[machine]=timeTable[machine]+processingTime;
					}
					else
					{
						timeTable[machine]= timeTable[machine-1]+processingTime;
					}
				}
				//cout<<"M "<<machine<<" Time "<<timeTable[machine]<<endl;
			}
		}

		fitness1=timeTable[m_machines-1];
		fitness[0]=(int)fitness1;
//		cout << solution.getFunction(0) << " " << solution.getFunction(1) << endl;
		return fitness;

}

void PFSP::solutionFitness(Subprob& solution) {

	//objective 1 = Evaluates the given solution with the total flow time criterion.

		for (int i=0;i<m_machines;i++) m_timeTable[i]=0;
		int j,z, job;
		int machine;
		int prev_machine=0;

		int first_gene=solution.getVectorSolution(0);
		m_timeTable[0]=m_processingtimes[0][first_gene];
		for (j=1;j<m_machines;j++)
		{
			m_timeTable[j]=m_timeTable[j-1]+m_processingtimes[j][first_gene];
		}

		double fitness0=m_timeTable[m_machines-1];
		for (z=1;z<m_jobs;z++)
		{
			job=solution.getVectorSolution(z);

			//machine 0 is always incremental, so:
			m_timeTable[0]+=m_processingtimes[0][job];
			prev_machine=m_timeTable[0];
			for (machine=1;machine<m_machines;machine++)
			{
				m_timeTable[machine]= MAX(prev_machine,m_timeTable[machine])+ m_processingtimes[machine][job];
				prev_machine=m_timeTable[machine];
			}

			fitness0+=m_timeTable[m_machines-1];
		}

		solution.setFunction((int)fitness0,1);

	//objective 2 = Make Span

		double fitness1=0;
		double timeTable[m_machines];

		for (int j=0;j<m_machines;j++)
		{
			timeTable[j]=0;
		}
		for (int z=0;z<m_jobs;z++)
		{
			int job=solution.getVectorSolution(z);
			//cout<<"Job "<<job<<endl;
			for (int machine=0;machine<m_machines;machine++)
			{
				double processingTime=m_processingtimes[machine][job];
				if (machine==0)
				{
					timeTable[machine]+=processingTime;
				}
				else
				{
					if (timeTable[machine-1]<timeTable[machine])
					{
						timeTable[machine]=timeTable[machine]+processingTime;
					}
					else
					{
						timeTable[machine]= timeTable[machine-1]+processingTime;
					}
				}
				//cout<<"M "<<machine<<" Time "<<timeTable[machine]<<endl;
			}
		}

		fitness1=timeTable[m_machines-1];
		solution.setFunction((int)fitness1,0);
//		cout << solution.getFunction(0) << " " << solution.getFunction(1) << endl;
}

int PFSP::GetProblemSize()
{
    return m_jobs;
}

/////////////////////////////////////////////////////////////////////////////////// PFSP tard /////////////////////////////////////////////////////////////

PFSP_tard::PFSP_tard(){
//	m_machines = machine;
//	m_jobs = numbVar;

}

PFSP_tard::~PFSP_tard(){
	for	(int i=0;i<m_machines;i++)
	delete [] m_processingtimes[i];
	delete [] m_processingtimes;
    delete [] m_timeTable;
    delete [] m_aux;
    delete [] m_DueDates;
}

void PFSP_tard::Initialize() {

	//BUILD JOB PROCESSING MATRIX
	//cout << "--> BUILDING JOB PROCESSING MATRIX" << endl;


}

void PFSP_tard::LoadInstance() {

	fstream dataIn;
	dataIn.open(("../FSSP/DD_Ta"+idxInstance+".txt").c_str(),ios::in);
	if(dataIn==NULL){
		cout << "Error opening file FSSP" << endl;
		exit(EXIT_FAILURE);
	}else{
		cout << "test instance found" << endl;
	}
	string line;
	int value;
	dataIn >> line; istringstream(line) >> value; numbVar = value; m_jobs = value; cout << numbVar << endl;
	dataIn >> line; istringstream(line) >> value; machine = value; m_machines = value; cout << machine << endl;

	m_processingtimes = new int*[m_machines];
	m_timeTable= new int[m_machines];
	m_aux= new int[m_jobs];
	for (int i=0;i<m_machines;i++)
	{
		m_processingtimes[i]= new int[m_jobs];
	}
	m_DueDates = new int[m_jobs];

	for(int i=0;i<m_jobs;i++){
		for(int j=0;j<m_machines;j++){
			dataIn >> line; istringstream(line) >> value; cout << value << " ";
			dataIn >> line; istringstream(line) >> value; cout << value << endl;
			m_processingtimes[j][i]=value;
		}
	}

	dataIn >> line;
	if(line=="Reldue"){
		cout << line << endl;
		for(int i=0;i<m_jobs;i++){
			dataIn >> line;
			dataIn >> line; istringstream(line) >> value; cout << value << endl; m_DueDates[i]=value;
			dataIn >> line; dataIn >> line;
		}
	}

	dataIn.close();
}

void PFSP_tard::Show() {
	//print the flow shop


		for(int i=0;i<m_machines;i++){
			for(int j=0;j<m_jobs;j++){
				cout << m_processingtimes[i][j] << " ";
			}
			cout << endl;
		}
		cout << "Due Dates" << endl;
		for(int j=0;j<m_jobs;j++){
			cout << m_DueDates[j] << " ";
		}
		cout << endl;
}

int * PFSP_tard::solutionFitness(int * gen){
	int * fitness = new int[numbObj];
	//objective 1 = Evaluates the given solution with the total flow time criterion.

		for (int i=0;i<m_machines;i++) m_timeTable[i]=0;
		int j,z, job;
		int machine;
		int tardness;
		int prev_machine=0;

		int first_gene=gen[0];
		m_timeTable[0]=m_processingtimes[0][first_gene];
		for (j=1;j<m_machines;j++)
		{
			m_timeTable[j]=m_timeTable[j-1]+m_processingtimes[j][first_gene];
		}
		tardness = MAX(m_timeTable[m_machines-1]-m_DueDates[first_gene],0);
		double fitness0=tardness;
		for (z=1;z<m_jobs;z++)
		{
			job=gen[z];

			//machine 0 is always incremental, so:
			m_timeTable[0]+=m_processingtimes[0][job];
			prev_machine=m_timeTable[0];
			for (machine=1;machine<m_machines;machine++)
			{
				m_timeTable[machine]= MAX(prev_machine,m_timeTable[machine])+ m_processingtimes[machine][job];
				prev_machine=m_timeTable[machine];
			}
			tardness = MAX(m_timeTable[m_machines-1]-m_DueDates[job],0);
			fitness0+=tardness;
		}


		fitness[1]=(int)fitness0;

	//objective 2 = Make Span

		double fitness1=0;
		double timeTable[m_machines];

		for (int j=0;j<m_machines;j++)
		{
			timeTable[j]=0;
		}
		for (int z=0;z<m_jobs;z++)
		{
			int job=gen[z];
			//cout<<"Job "<<job<<endl;
			for (int machine=0;machine<m_machines;machine++)
			{
				double processingTime=m_processingtimes[machine][job];
				if (machine==0)
				{
					timeTable[machine]+=processingTime;
				}
				else
				{
					if (timeTable[machine-1]<timeTable[machine])
					{
						timeTable[machine]=timeTable[machine]+processingTime;
					}
					else
					{
						timeTable[machine]= timeTable[machine-1]+processingTime;
					}
				}
				//cout<<"M "<<machine<<" Time "<<timeTable[machine]<<endl;
			}
		}

		fitness1=timeTable[m_machines-1];
		fitness[0]=(int)fitness1;
//		cout << solution.getFunction(0) << " " << solution.getFunction(1) << endl;
		return fitness;

}

void PFSP_tard::solutionFitness(Subprob& solution) {

	//objective 1 = Evaluates the given solution with the total flow time criterion.

		for (int i=0;i<m_machines;i++) m_timeTable[i]=0;
		int j,z, job;
		int machine;
		int prev_machine=0;
		int tardness;
		int first_gene=solution.getVectorSolution(0);
		m_timeTable[0]=m_processingtimes[0][first_gene];
		for (j=1;j<m_machines;j++)
		{
			m_timeTable[j]=m_timeTable[j-1]+m_processingtimes[j][first_gene];
		}
		tardness = MAX(m_timeTable[m_machines-1]-m_DueDates[first_gene],0);
		double fitness0=tardness;
		for (z=1;z<m_jobs;z++)
		{
			job=solution.getVectorSolution(z);

			//machine 0 is always incremental, so:
			m_timeTable[0]+=m_processingtimes[0][job];
			prev_machine=m_timeTable[0];
			for (machine=1;machine<m_machines;machine++)
			{
				m_timeTable[machine]= MAX(prev_machine,m_timeTable[machine])+ m_processingtimes[machine][job];
				prev_machine=m_timeTable[machine];
			}

			tardness = MAX(m_timeTable[m_machines-1]-m_DueDates[job],0);
			fitness0+=tardness;
		}

		solution.setFunction((int)fitness0,1);

	//objective 2 = Make Span

		double fitness1=0;
		double timeTable[m_machines];

		for (int j=0;j<m_machines;j++)
		{
			timeTable[j]=0;
		}
		for (int z=0;z<m_jobs;z++)
		{
			int job=solution.getVectorSolution(z);
			//cout<<"Job "<<job<<endl;
			for (int machine=0;machine<m_machines;machine++)
			{
				double processingTime=m_processingtimes[machine][job];
				if (machine==0)
				{
					timeTable[machine]+=processingTime;
				}
				else
				{
					if (timeTable[machine-1]<timeTable[machine])
					{
						timeTable[machine]=timeTable[machine]+processingTime;
					}
					else
					{
						timeTable[machine]= timeTable[machine-1]+processingTime;
					}
				}
				//cout<<"M "<<machine<<" Time "<<timeTable[machine]<<endl;
			}
		}

		fitness1=timeTable[m_machines-1];
		solution.setFunction((int)fitness1,0);
//		cout << solution.getFunction(0) << " " << solution.getFunction(1) << endl;
}

int PFSP_tard::GetProblemSize()
{
    return m_jobs;
}

//////////////////////////////////////////////////////////////// PFSP 3 obj /////////////////////////////////////////////////////////////////////////////////////////////////////
PFSP_3ob::PFSP_3ob(){
//	m_machines = machine;
//	m_jobs = numbVar;

}

PFSP_3ob::~PFSP_3ob(){
	for	(int i=0;i<m_machines;i++)
	delete [] m_processingtimes[i];
	delete [] m_processingtimes;
    delete [] m_timeTable;
    delete [] m_aux;
    delete [] m_DueDates;
}

void PFSP_3ob::Initialize() {

	//BUILD JOB PROCESSING MATRIX
	//cout << "--> BUILDING JOB PROCESSING MATRIX" << endl;


}

void PFSP_3ob::LoadInstance() {

	fstream dataIn;
	dataIn.open(("../FSSP/DD_Ta"+idxInstance+".txt").c_str(),ios::in);
	if(dataIn==NULL){
		cout << "Error opening file FSSP" << endl;
		exit(EXIT_FAILURE);
	}else{
		cout << "test instance found" << endl;
	}
	string line;
	int value;
	dataIn >> line; istringstream(line) >> value; numbVar = value; m_jobs = value; cout << numbVar << endl;
	dataIn >> line; istringstream(line) >> value; machine = value; m_machines = value; cout << machine << endl;

	m_processingtimes = new int*[m_machines];
	m_timeTable= new int[m_machines];
	m_aux= new int[m_jobs];
	for (int i=0;i<m_machines;i++)
	{
		m_processingtimes[i]= new int[m_jobs];
	}
	m_DueDates = new int[m_jobs];

	for(int i=0;i<m_jobs;i++){
		for(int j=0;j<m_machines;j++){
			dataIn >> line; istringstream(line) >> value; cout << value << " ";
			dataIn >> line; istringstream(line) >> value; cout << value << endl;
			m_processingtimes[j][i]=value;
		}
	}

	dataIn >> line;
	if(line=="Reldue"){
		cout << line << endl;
		for(int i=0;i<m_jobs;i++){
			dataIn >> line;
			dataIn >> line; istringstream(line) >> value; cout << value << endl; m_DueDates[i]=value;
			dataIn >> line; dataIn >> line;
		}
	}

	dataIn.close();
}

void PFSP_3ob::Show() {
	//print the flow shop


		for(int i=0;i<m_machines;i++){
			for(int j=0;j<m_jobs;j++){
				cout << m_processingtimes[i][j] << " ";
			}
			cout << endl;
		}
		cout << "Due Dates" << endl;
		for(int j=0;j<m_jobs;j++){
			cout << m_DueDates[j] << " ";
		}
		cout << endl;
}

int * PFSP_3ob::solutionFitness(int * gen){
	int * fitness = new int[numbObj];
	//objective 1 = Evaluates the given solution with the total flow time criterion.

		for (int i=0;i<m_machines;i++) m_timeTable[i]=0;
		int j,z, job;
		int machine;
		int tardness=0, tft=0;
		int prev_machine=0;

		int first_gene=gen[0];
		m_timeTable[0]=m_processingtimes[0][first_gene];
		for (j=1;j<m_machines;j++)
		{
			m_timeTable[j]=m_timeTable[j-1]+m_processingtimes[j][first_gene];
		}
		tardness = MAX(m_timeTable[m_machines-1]-m_DueDates[first_gene],0);

		tft = m_timeTable[m_machines-1];
		for (z=1;z<m_jobs;z++)
		{
			job=gen[z];

			//machine 0 is always incremental, so:
			m_timeTable[0]+=m_processingtimes[0][job];
			prev_machine=m_timeTable[0];
			for (machine=1;machine<m_machines;machine++)
			{
				m_timeTable[machine]= MAX(prev_machine,m_timeTable[machine])+ m_processingtimes[machine][job];
				prev_machine=m_timeTable[machine];
			}
			tardness += MAX(m_timeTable[m_machines-1]-m_DueDates[job],0);

			tft += m_timeTable[m_machines-1];
		}

		fitness[1]=(int)tft;
		fitness[2]=(int)tardness;

	//objective 2 = Make Span

		double fitness1=0;
		double timeTable[m_machines];

		for (int j=0;j<m_machines;j++)
		{
			timeTable[j]=0;
		}
		for (int z=0;z<m_jobs;z++)
		{
			int job=gen[z];
			//cout<<"Job "<<job<<endl;
			for (int machine=0;machine<m_machines;machine++)
			{
				double processingTime=m_processingtimes[machine][job];
				if (machine==0)
				{
					timeTable[machine]+=processingTime;
				}
				else
				{
					if (timeTable[machine-1]<timeTable[machine])
					{
						timeTable[machine]=timeTable[machine]+processingTime;
					}
					else
					{
						timeTable[machine]= timeTable[machine-1]+processingTime;
					}
				}
				//cout<<"M "<<machine<<" Time "<<timeTable[machine]<<endl;
			}
		}

		fitness1=timeTable[m_machines-1];
		fitness[0]=(int)fitness1;
//		cout << solution.getFunction(0) << " " << solution.getFunction(1) << endl;
		return fitness;

}

void PFSP_3ob::solutionFitness(Subprob& solution) {

	//objective 1 = Evaluates the given solution with the total flow time criterion.

		for (int i=0;i<m_machines;i++) m_timeTable[i]=0;
		int j,z, job;
		int machine;
		int prev_machine=0;
		int tardness=0, tft=0;
		int first_gene=solution.getVectorSolution(0);
		m_timeTable[0]=m_processingtimes[0][first_gene];
		for (j=1;j<m_machines;j++)
		{
			m_timeTable[j]=m_timeTable[j-1]+m_processingtimes[j][first_gene];
		}
		tardness = MAX(m_timeTable[m_machines-1]-m_DueDates[first_gene],0);

		tft = m_timeTable[m_machines-1];
		for (z=1;z<m_jobs;z++)
		{
			job=solution.getVectorSolution(z);

			//machine 0 is always incremental, so:
			m_timeTable[0]+=m_processingtimes[0][job];
			prev_machine=m_timeTable[0];
			for (machine=1;machine<m_machines;machine++)
			{
				m_timeTable[machine]= MAX(prev_machine,m_timeTable[machine])+ m_processingtimes[machine][job];
				prev_machine=m_timeTable[machine];
			}

			tardness += MAX(m_timeTable[m_machines-1]-m_DueDates[job],0);

			tft += m_timeTable[m_machines-1];
		}

		solution.setFunction((int)tft,1);
		solution.setFunction((int)tardness,2);
	//objective 2 = Make Span

		double fitness1=0;
		double timeTable[m_machines];

		for (int j=0;j<m_machines;j++)
		{
			timeTable[j]=0;
		}
		for (int z=0;z<m_jobs;z++)
		{
			int job=solution.getVectorSolution(z);
			//cout<<"Job "<<job<<endl;
			for (int machine=0;machine<m_machines;machine++)
			{
				double processingTime=m_processingtimes[machine][job];
				if (machine==0)
				{
					timeTable[machine]+=processingTime;
				}
				else
				{
					if (timeTable[machine-1]<timeTable[machine])
					{
						timeTable[machine]=timeTable[machine]+processingTime;
					}
					else
					{
						timeTable[machine]= timeTable[machine-1]+processingTime;
					}
				}
				//cout<<"M "<<machine<<" Time "<<timeTable[machine]<<endl;
			}
		}

		fitness1=timeTable[m_machines-1];
		solution.setFunction((int)fitness1,0);
//		cout << solution.getFunction(0) << " " << solution.getFunction(1) << endl;
}

int PFSP_3ob::GetProblemSize()
{
    return m_jobs;
}

/////////////////////////////////////////////////////////////////////////////////// PFSP-SDST cmax and weight tard/////////////////////////////////////////////////////////////

PFSP_SDST::PFSP_SDST(){
//	m_machines = machine;
//	m_jobs = numbVar;

}

PFSP_SDST::~PFSP_SDST(){
	for	(int i=0;i<m_machines;i++)
	delete [] m_processingtimes[i];
	delete [] m_processingtimes;
    delete [] m_timeTable;
    delete [] m_aux;


}

void PFSP_SDST::Initialize() {

	//BUILD JOB PROCESSING MATRIX
	//cout << "--> BUILDING JOB PROCESSING MATRIX" << endl;


}

void PFSP_SDST::LoadInstance() {
	string st = static_cast<ostringstream*>( &(ostringstream() << benchmark) )->str();
	fstream dataIn;
	dataIn.open(("../PFSP_SDST/ssd"+st+"/DD_SDST"+st+"_ta"+idxInstance+".txt").c_str(),ios::in);
	if(dataIn==NULL){
		cout << "Error opening file FSSP" << endl;
		exit(EXIT_FAILURE);
	}else{
		cout << "test instance found" << endl;
	}
	string line;
	int value;
	dataIn >> line; istringstream(line) >> value; numbVar = value; m_jobs = value; cout << numbVar << endl;
	dataIn >> line; istringstream(line) >> value; machine = value; m_machines = value; cout << machine << endl;

	m_processingtimes = new int*[m_machines];
	m_timeTable= new int[m_machines];
	m_aux= new int[m_jobs];
	for (int i=0;i<m_machines;i++)
	{
		m_processingtimes[i]= new int[m_jobs];
	}

	m_SDST = new int **[m_machines];
	for (int i=0;i<m_machines;i++)
	{
		m_SDST[i] = new int * [m_jobs];
		for (int j=0;j<m_jobs;j++)
		{
			m_SDST[i][j] = new int[m_jobs];
		}
	}

	m_DueDates = new int[m_jobs];
	m_weights = new int[m_jobs];

	for(int i=0;i<m_jobs;i++){
		for(int j=0;j<m_machines;j++){
			dataIn >> line; istringstream(line) >> value;
			dataIn >> line; istringstream(line) >> value;
			m_processingtimes[j][i]=value;
		}
	}

	dataIn >> line;
	if(line=="SSD"){

		for(int k=0;k<m_machines;k++){
			dataIn >> line;
			cout <<"M"<<k<<endl;
			for(int i=0;i<m_jobs;i++){
				for(int j=0;j<m_jobs;j++){
					dataIn >> line; istringstream(line) >> value; m_SDST[k][i][j]=value; //cout << m_SDST[k][i][j] << " ";
				}
				//cout << endl;
			}

		}
	}
	dataIn >> line;
	if(line=="Reldue"){

		for(int i=0;i<m_jobs;i++){
			dataIn >> line;
			dataIn >> line; istringstream(line) >> value;  m_DueDates[i]=value;
			dataIn >> line; dataIn >> line; istringstream(line) >> value; m_weights[i]=value;
		}
	}

	dataIn.close();
}

void PFSP_SDST::Show() {
	//print the flow shop


		for(int i=0;i<m_machines;i++){
			for(int j=0;j<m_jobs;j++){
				cout << m_processingtimes[i][j] << " ";
			}
			cout << endl;
		}
//		cout << "Due Dates" << endl;
//		for(int j=0;j<m_jobs;j++){
//			cout << m_DueDates[j] << " weigth " << m_weights[j] << endl;
//		}
//
//		for(int k=0;k<m_machines;k++){
//			cout <<"M"<<k<<endl;
//			for(int i=0;i<m_jobs;i++){
//				for(int j=0;j<m_jobs;j++){
//					cout << m_SDST[k][i][j] << " ";
//				}
//				cout << endl;
//			}
//		}
}

int * PFSP_SDST::solutionFitness(int * gen){
	int * fitness = new int[numbObj];
	//objective 1 = Evaluates the given solution with the total tardiness criterion
		for (int i=0;i<m_machines;i++) m_timeTable[i]=0;
		int j,z, job, prev_job;
		int machine;
		int tardness;
		int prev_machine=0;

		//first job, which it does not have SDST
		int first_gene=gen[0];
		m_timeTable[0]=m_processingtimes[0][first_gene]; //first machine
		for (j=1;j<m_machines;j++)
		{
			m_timeTable[j]=m_timeTable[j-1]+m_processingtimes[j][first_gene];
		}
		tardness = MAX(m_timeTable[m_machines-1]-m_DueDates[first_gene],0);
		double fitness0=tardness;
		//end first job, which it does not have SDST

		//for the second job until the last one, consider the SDST
		for (z=1;z<m_jobs;z++)
		{
			job=gen[z];
			prev_job=gen[z-1];
			//machine 0 is always incremental, so:
			m_timeTable[0]+=m_processingtimes[0][job]+m_SDST[0][prev_job][job];
			prev_machine=m_timeTable[0];
			for (machine=1;machine<m_machines;machine++)
			{
				m_timeTable[machine]+=m_SDST[machine][prev_job][job];
				m_timeTable[machine]= MAX(prev_machine,m_timeTable[machine])+ m_processingtimes[machine][job];
				prev_machine=m_timeTable[machine];
			}
			tardness = m_weights[job]*MAX((m_timeTable[m_machines-1]-m_DueDates[job]),0);
			fitness0+=tardness;
		}


		fitness[1]=(int)fitness0;

	//objective 2 = Make Span

		double fitness1=0;
		double timeTable[m_machines];

		for (int j=0;j<m_machines;j++)
		{
			timeTable[j]=0;
		}
		for (int z=0;z<m_jobs;z++)
		{
			job=gen[z];
			//cout<<"Job "<<job<<endl;
			for (int machine=0;machine<m_machines;machine++)
			{
				double processingTime=m_processingtimes[machine][job];
				if (machine==0)
				{
					timeTable[machine]+=processingTime;

					if(z!=0)
					{
						prev_job=gen[z-1];
						timeTable[machine]+=m_SDST[machine][prev_job][job];
					}
				}
				else
				{
					if (z!=0)
					{
						prev_job = gen[z-1];
						timeTable[machine]+=m_SDST[machine][prev_job][job];
					}

					if (timeTable[machine-1]<timeTable[machine])
					{
						timeTable[machine]=timeTable[machine]+processingTime;
					}
					else
					{
						timeTable[machine]= timeTable[machine-1]+processingTime;
					}
				}
				//cout<<"M "<<machine<<" Time "<<timeTable[machine]<<endl;
			}
		}

		fitness1=timeTable[m_machines-1];
		fitness[0]=(int)fitness1;
		//cout << fitness[0] << " " << fitness[1] << endl; getchar();
		return fitness;

}

void PFSP_SDST::solutionFitness(Subprob& solution) {
	//objective 1 = Evaluates the given solution with the total tardiness criterion

		for (int i=0;i<m_machines;i++) m_timeTable[i]=0;
		int j,z, job, prev_job;
		int machine;
		int prev_machine=0;
		int tardness;
		int first_gene=solution.getVectorSolution(0);
		m_timeTable[0]=m_processingtimes[0][first_gene];
		for (j=1;j<m_machines;j++)
		{
			m_timeTable[j]=m_timeTable[j-1]+m_processingtimes[j][first_gene];
		}
		tardness = MAX(m_timeTable[m_machines-1]-m_DueDates[first_gene],0);
		double fitness0=tardness;
		//end first job, which it does not have SDST

		//for the second job until the last one, consider the SDST
		for (z=1;z<m_jobs;z++)
		{
			job=solution.getVectorSolution(z);
			prev_job=solution.getVectorSolution(z-1);
			//machine 0 is always incremental, so:
			m_timeTable[0]+=m_processingtimes[0][job]+m_SDST[0][prev_job][job];
			prev_machine=m_timeTable[0];
			for (machine=1;machine<m_machines;machine++)
			{
				m_timeTable[machine]+=m_SDST[machine][prev_job][job];
				m_timeTable[machine]= MAX(prev_machine,m_timeTable[machine])+ m_processingtimes[machine][job];
				prev_machine=m_timeTable[machine];
			}
			tardness = m_weights[job]*MAX((m_timeTable[m_machines-1]-m_DueDates[job]),0);
			fitness0+=tardness;
		}

		solution.setFunction((int)fitness0,1);

	//objective 2 = Make Span

		double fitness1=0;
		double timeTable[m_machines];

		for (int j=0;j<m_machines;j++)
		{
			timeTable[j]=0;
		}
		for (int z=0;z<m_jobs;z++)
		{
			job=solution.getVectorSolution(z);
			//cout<<"Job "<<job<<endl;
			for (int machine=0;machine<m_machines;machine++)
			{
				double processingTime=m_processingtimes[machine][job];
				if (machine==0)
				{
					timeTable[machine]+=processingTime;

					if(z!=0)
					{
						prev_job=solution.getVectorSolution(z-1);
						timeTable[machine]+=m_SDST[machine][prev_job][job];
					}
				}
				else
				{
					if (z!=0)
					{
						prev_job = solution.getVectorSolution(z-1);
						timeTable[machine]+=m_SDST[machine][prev_job][job];
					}

					if (timeTable[machine-1]<timeTable[machine])
					{
						timeTable[machine]=timeTable[machine]+processingTime;
					}
					else
					{
						timeTable[machine]= timeTable[machine-1]+processingTime;
					}
				}
				//cout<<"M "<<machine<<" Time "<<timeTable[machine]<<endl;
			}
		}


		fitness1=timeTable[m_machines-1];
		solution.setFunction((int)fitness1,0);
	//	cout << solution.getFunction(0) << " " << solution.getFunction(1) << endl; getchar();
}

int PFSP_SDST::GetProblemSize()
{
    return m_jobs;
}


int PFSP_SDST::EvalCmax(int * genes, int size)
{
	int fitness;
	int job, prev_job;
//	int prev_machine=0;
	int partition=size;

		//objective Make Span SDST <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		double timeTable[m_machines];

		for (int j=0;j<m_machines;j++)
		{
			timeTable[j]=0;
		}
		for (int z=0;z<partition;z++)
		{
			job=genes[z];
			//cout<<"Job "<<job<<endl;
			for (int machine=0;machine<m_machines;machine++)
			{
				double processingTime=m_processingtimes[machine][job];
				if (machine==0)
				{
					timeTable[machine]+=processingTime;

					if(z!=0)
					{
						prev_job=genes[z-1];
						timeTable[machine]+=m_SDST[machine][prev_job][job];
					}
				}
				else
				{
					if (z!=0)
					{
						prev_job = genes[z-1];
						timeTable[machine]+=m_SDST[machine][prev_job][job];
					}

					if (timeTable[machine-1]<timeTable[machine])
					{
						timeTable[machine]=timeTable[machine]+processingTime;
					}
					else
					{
						timeTable[machine]= timeTable[machine-1]+processingTime;
					}
				}
				//cout<<"M "<<machine<<" Time "<<timeTable[machine]<<endl;
			}
		}

		fitness=timeTable[m_machines-1];

		//cout << fitness << " " << endl; getchar();

    return fitness;
}

int PFSP_SDST::EvalTWT(int * genes, int size)
{
	//objective 1 = Evaluates the given solution with the total tardiness criterion
	int fitness;
	for (int i=0;i<m_machines;i++) m_timeTable[i]=0;
	int j,z, job, prev_job;
	int machine;
	int tardness;
	int prev_machine=0;

	//first job, which it does not have SDST
	int first_gene=genes[0];

	m_timeTable[0]=m_processingtimes[0][first_gene]; //first machine
	for (j=1;j<m_machines;j++)
	{
		m_timeTable[j]=m_timeTable[j-1]+m_processingtimes[j][first_gene];
	}
	tardness = MAX(m_timeTable[m_machines-1]-m_DueDates[first_gene],0);

	fitness=tardness;
	//end first job, which it does not have SDST

	//for the second job until the last one, consider the SDST
	for (z=1;z<size;z++)
	{
		job=genes[z];
		prev_job=genes[z-1];
		//machine 0 is always incremental, so:
		m_timeTable[0]+=m_processingtimes[0][job]+m_SDST[0][prev_job][job];
		prev_machine=m_timeTable[0];
		for (machine=1;machine<m_machines;machine++)
		{
			m_timeTable[machine]+=m_SDST[machine][prev_job][job];
			m_timeTable[machine]= MAX(prev_machine,m_timeTable[machine])+ m_processingtimes[machine][job];
			prev_machine=m_timeTable[machine];
		}
		tardness = m_weights[job]*MAX((m_timeTable[m_machines-1]-m_DueDates[job]),0);
		fitness+=tardness;
	}
	return fitness;
}

/////////////////////////////////////////////////////////////////////////////////// PFSP-SDST cmax and total flow time/////////////////////////////////////////////////////////////

SDST_Cmax_TFT::SDST_Cmax_TFT(){
//	m_machines = machine;
//	m_jobs = numbVar;

}

SDST_Cmax_TFT::~SDST_Cmax_TFT(){
	for	(int i=0;i<m_machines;i++)
	delete [] m_processingtimes[i];
	delete [] m_processingtimes;
    delete [] m_timeTable;
    delete [] m_aux;


}

void SDST_Cmax_TFT::Initialize() {

	//BUILD JOB PROCESSING MATRIX
	//cout << "--> BUILDING JOB PROCESSING MATRIX" << endl;


}

void SDST_Cmax_TFT::LoadInstance() {
	string st = static_cast<ostringstream*>( &(ostringstream() << benchmark) )->str();
	fstream dataIn;
	dataIn.open(("../PFSP_SDST/ssd"+st+"/DD_SDST"+st+"_ta"+idxInstance+".txt").c_str(),ios::in);
	if(dataIn==NULL){
		cout << "Error opening file FSSP" << endl;
		exit(EXIT_FAILURE);
	}else{
		cout << "test instance found" << endl;
	}
	string line;
	int value;
	dataIn >> line; istringstream(line) >> value; numbVar = value; m_jobs = value; cout << numbVar << endl;
	dataIn >> line; istringstream(line) >> value; machine = value; m_machines = value; cout << machine << endl;

	m_processingtimes = new int*[m_machines];
	m_timeTable= new int[m_machines];
	m_aux= new int[m_jobs];
	for (int i=0;i<m_machines;i++)
	{
		m_processingtimes[i]= new int[m_jobs];
	}

	m_SDST = new int **[m_machines];
	for (int i=0;i<m_machines;i++)
	{
		m_SDST[i] = new int * [m_jobs];
		for (int j=0;j<m_jobs;j++)
		{
			m_SDST[i][j] = new int[m_jobs];
		}
	}

	m_DueDates = new int[m_jobs];
	m_weights = new int[m_jobs];

	for(int i=0;i<m_jobs;i++){
		for(int j=0;j<m_machines;j++){
			dataIn >> line; istringstream(line) >> value;
			dataIn >> line; istringstream(line) >> value;
			m_processingtimes[j][i]=value;
		}
	}

	dataIn >> line;
	if(line=="SSD"){

		for(int k=0;k<m_machines;k++){
			dataIn >> line;
			cout <<"M"<<k<<endl;
			for(int i=0;i<m_jobs;i++){
				for(int j=0;j<m_jobs;j++){
					dataIn >> line; istringstream(line) >> value; m_SDST[k][i][j]=value; //cout << m_SDST[k][i][j] << " ";
				}
				//cout << endl;
			}

		}
	}
	dataIn >> line;
	if(line=="Reldue"){

		for(int i=0;i<m_jobs;i++){
			dataIn >> line;
			dataIn >> line; istringstream(line) >> value;  m_DueDates[i]=value;
			dataIn >> line; dataIn >> line; istringstream(line) >> value; m_weights[i]=value;
		}
	}

	dataIn.close();
}

void SDST_Cmax_TFT::Show() {
	//print the flow shop


		for(int i=0;i<m_machines;i++){
			for(int j=0;j<m_jobs;j++){
				cout << m_processingtimes[i][j] << " ";
			}
			cout << endl;
		}
//		cout << "Due Dates" << endl;
//		for(int j=0;j<m_jobs;j++){
//			cout << m_DueDates[j] << " weigth " << m_weights[j] << endl;
//		}
//
//		for(int k=0;k<m_machines;k++){
//			cout <<"M"<<k<<endl;
//			for(int i=0;i<m_jobs;i++){
//				for(int j=0;j<m_jobs;j++){
//					cout << m_SDST[k][i][j] << " ";
//				}
//				cout << endl;
//			}
//		}
}

int * SDST_Cmax_TFT::solutionFitness(int * gen){
	int * fitness = new int[numbObj];
	//objective 1 = Evaluates the given solution with the total flow time

	for (int i=0;i<m_machines;i++) m_timeTable[i]=0;
		int j,z, job, prev_job;
		int machine;
		int prev_machine=0;

		//first job, which it does not have SDST
		int first_gene=gen[0];
		m_timeTable[0]=m_processingtimes[0][first_gene]; //first machine
		for (j=1;j<m_machines;j++)
		{
			m_timeTable[j]=m_timeTable[j-1]+m_processingtimes[j][first_gene];
		}

		double fitness0= m_timeTable[m_machines-1]-m_DueDates[first_gene];
		//end first job, which it does not have SDST

		//for the second job until the last one, consider the SDST
		for (z=1;z<m_jobs;z++)
		{
			job=gen[z];
			prev_job=gen[z-1];
			//machine 0 is always incremental, so:
			m_timeTable[0]+=m_processingtimes[0][job]+m_SDST[0][prev_job][job];
			prev_machine=m_timeTable[0];
			for (machine=1;machine<m_machines;machine++)
			{
				m_timeTable[machine]+=m_SDST[machine][prev_job][job];
				m_timeTable[machine]= MAX(prev_machine,m_timeTable[machine])+ m_processingtimes[machine][job];
				prev_machine=m_timeTable[machine];
			}
			fitness0+=m_timeTable[m_machines-1];;
		}


		fitness[1]=(int)fitness0;

	//objective 2 = Make Span

		double fitness1=0;
		double timeTable[m_machines];

		for (int j=0;j<m_machines;j++)
		{
			timeTable[j]=0;
		}
		for (int z=0;z<m_jobs;z++)
		{
			job=gen[z];
			//cout<<"Job "<<job<<endl;
			for (int machine=0;machine<m_machines;machine++)
			{
				double processingTime=m_processingtimes[machine][job];
				if (machine==0)
				{
					timeTable[machine]+=processingTime;

					if(z!=0)
					{
						prev_job=gen[z-1];
						timeTable[machine]+=m_SDST[machine][prev_job][job];
					}
				}
				else
				{
					if (z!=0)
					{
						prev_job = gen[z-1];
						timeTable[machine]+=m_SDST[machine][prev_job][job];
					}

					if (timeTable[machine-1]<timeTable[machine])
					{
						timeTable[machine]=timeTable[machine]+processingTime;
					}
					else
					{
						timeTable[machine]= timeTable[machine-1]+processingTime;
					}
				}
				//cout<<"M "<<machine<<" Time "<<timeTable[machine]<<endl;
			}
		}

		fitness1=timeTable[m_machines-1];
		fitness[0]=(int)fitness1;
		//cout << fitness[0] << " " << fitness[1] << endl; getchar();
		return fitness;

}

void SDST_Cmax_TFT::solutionFitness(Subprob& solution) {
	//objective 1 = Evaluates the given solution with the total flow time criterion

		for (int i=0;i<m_machines;i++) m_timeTable[i]=0;
		int j,z, job, prev_job;
		int machine;
		int prev_machine=0;

		int first_gene=solution.getVectorSolution(0);
		m_timeTable[0]=m_processingtimes[0][first_gene];
		for (j=1;j<m_machines;j++)
		{
			m_timeTable[j]=m_timeTable[j-1]+m_processingtimes[j][first_gene];
		}

		double fitness0=m_timeTable[m_machines-1];
		//end first job, which it does not have SDST

		//for the second job until the last one, consider the SDST
		for (z=1;z<m_jobs;z++)
		{
			job=solution.getVectorSolution(z);
			prev_job=solution.getVectorSolution(z-1);
			//machine 0 is always incremental, so:
			m_timeTable[0]+=m_processingtimes[0][job]+m_SDST[0][prev_job][job];
			prev_machine=m_timeTable[0];
			for (machine=1;machine<m_machines;machine++)
			{
				m_timeTable[machine]+=m_SDST[machine][prev_job][job];
				m_timeTable[machine]= MAX(prev_machine,m_timeTable[machine])+ m_processingtimes[machine][job];
				prev_machine=m_timeTable[machine];
			}
			fitness0+= m_timeTable[m_machines-1];

		}

		solution.setFunction((int)fitness0,1);

	//objective 2 = Make Span

		double fitness1=0;
		double timeTable[m_machines];

		for (int j=0;j<m_machines;j++)
		{
			timeTable[j]=0;
		}
		for (int z=0;z<m_jobs;z++)
		{
			job=solution.getVectorSolution(z);
			//cout<<"Job "<<job<<endl;
			for (int machine=0;machine<m_machines;machine++)
			{
				double processingTime=m_processingtimes[machine][job];
				if (machine==0)
				{
					timeTable[machine]+=processingTime;

					if(z!=0)
					{
						prev_job=solution.getVectorSolution(z-1);
						timeTable[machine]+=m_SDST[machine][prev_job][job];
					}
				}
				else
				{
					if (z!=0)
					{
						prev_job = solution.getVectorSolution(z-1);
						timeTable[machine]+=m_SDST[machine][prev_job][job];
					}

					if (timeTable[machine-1]<timeTable[machine])
					{
						timeTable[machine]=timeTable[machine]+processingTime;
					}
					else
					{
						timeTable[machine]= timeTable[machine-1]+processingTime;
					}
				}
				//cout<<"M "<<machine<<" Time "<<timeTable[machine]<<endl;
			}
		}


		fitness1=timeTable[m_machines-1];
		solution.setFunction((int)fitness1,0);
	//	cout << solution.getFunction(0) << " " << solution.getFunction(1) << endl; getchar();
}

int SDST_Cmax_TFT::GetProblemSize()
{
    return m_jobs;
}


int SDST_Cmax_TFT::EvalCmax(int * genes, int size)
{
	int fitness;
	int job, prev_job;
//	int prev_machine=0;
	int partition=size;

		//objective Make Span SDST <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		double timeTable[m_machines];

		for (int j=0;j<m_machines;j++)
		{
			timeTable[j]=0;
		}
		for (int z=0;z<partition;z++)
		{
			job=genes[z];
			//cout<<"Job "<<job<<endl;
			for (int machine=0;machine<m_machines;machine++)
			{
				double processingTime=m_processingtimes[machine][job];
				if (machine==0)
				{
					timeTable[machine]+=processingTime;

					if(z!=0)
					{
						prev_job=genes[z-1];
						timeTable[machine]+=m_SDST[machine][prev_job][job];
					}
				}
				else
				{
					if (z!=0)
					{
						prev_job = genes[z-1];
						timeTable[machine]+=m_SDST[machine][prev_job][job];
					}

					if (timeTable[machine-1]<timeTable[machine])
					{
						timeTable[machine]=timeTable[machine]+processingTime;
					}
					else
					{
						timeTable[machine]= timeTable[machine-1]+processingTime;
					}
				}
				//cout<<"M "<<machine<<" Time "<<timeTable[machine]<<endl;
			}
		}

		fitness=timeTable[m_machines-1];

		//cout << fitness << " " << endl; getchar();

    return fitness;
}

int SDST_Cmax_TFT::EvalTFT(int * genes, int size)
{
	//objective 1 = Evaluates the given solution with the total tardiness criterion
	int fitness;
	for (int i=0;i<m_machines;i++) m_timeTable[i]=0;
	int j,z, job, prev_job;
	int machine;
	int prev_machine=0;

	//first job, which it does not have SDST
	int first_gene=genes[0];

	m_timeTable[0]=m_processingtimes[0][first_gene]; //first machine
	for (j=1;j<m_machines;j++)
	{
		m_timeTable[j]=m_timeTable[j-1]+m_processingtimes[j][first_gene];
	}

	fitness=m_timeTable[m_machines-1];
	//end first job, which it does not have SDST

	//for the second job until the last one, consider the SDST
	for (z=1;z<size;z++)
	{
		job=genes[z];
		prev_job=genes[z-1];
		//machine 0 is always incremental, so:
		m_timeTable[0]+=m_processingtimes[0][job]+m_SDST[0][prev_job][job];
		prev_machine=m_timeTable[0];
		for (machine=1;machine<m_machines;machine++)
		{
			m_timeTable[machine]+=m_SDST[machine][prev_job][job];
			m_timeTable[machine]= MAX(prev_machine,m_timeTable[machine])+ m_processingtimes[machine][job];
			prev_machine=m_timeTable[machine];
		}
		fitness+=m_timeTable[m_machines-1];
	}
	return fitness;
}
