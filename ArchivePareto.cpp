#include"ArchivePareto.h"
#include"Global.h"


ArchivePareto::ArchivePareto() {
	sizeEP=0; sizeAux=0; add=false;

}

ArchivePareto::~ArchivePareto() {
}

void ArchivePareto::UpdateExternalPareto(Subprob solution) {
	add=false;
	int auxSol[objMax];
	 //add the new solutions to EP
    	for(int j=0;j<numbObj;j++){
			auxSol[j]=solution.getFunction(j);
		}
//		sizeEP++;

		//EP update
			for(int i=0;i<sizeEP;i++){

					int c1=0;
					int c2=numbObj;
					for(int k=0;k<numbObj;k++){
						if(EP[i][k]<auxSol[k]){
							c1++;
						}else{
							if(EP[i][k]==auxSol[k]){
								c2--;
							}
						}
					}
					//cout<< "values of c1 e c2 " << c1 << " " << c2 << endl;
					if(c1==0){
						if(c2==0){ 					//they are equal, thus, delete the first, set a flag that determine that the first is dominated
							for(int k=0;k<numbObj;k++){
								EP[i][k]=intMax;
								add=true;
							}
						}else{ 						//dominates, j dominates i, thus, set a flag of i that determine that i is dominated
							for(int k=0;k<numbObj;k++){
								EP[i][k]=intMax;
								add=true;
							}
						}
					}else{
						if(c1>0 and c1<c2){ 		//non dominated, do not delete anything
								add=true;
						}else{ 						//dominated, the i dominated the j, set a flag in j that determine that j is dominated
							for(int k=0;k<numbObj;k++){
								auxSol[k]=intMax;
								add=false;
							}
						}
					}

			}
			sizeAux=0;
			for(int i=0;i<sizeEP;i++){
				if(EP[i][0]!=intMax){ 				//if the solution is non dominated add it in the auxEP
					for(int k=0;k<numbObj;k++){
						auxEP[sizeAux][k]=EP[i][k];
					}
					sizeAux++;
				}
			}

			if(auxSol[0]!=intMax){

				for(int k=0;k<numbObj;k++){
					auxEP[sizeAux][k]=auxSol[k];
				}
				sizeAux++;
			}

		for(int i=0;i<sizeAux;i++){ 			//reorganize the external pareto (EP) add senquantialy the solution of auxEP in EP
			for(int k=0;k<numbObj;k++){
				EP[i][k]=auxEP[i][k];
			}
		}
		sizeEP=sizeAux; //set the size of EP
//		cout << "current size EP " << sizeAux << endl;



}





void ArchivePareto::ShowExternalPareto(int gen) {

		string g = static_cast<ostringstream*>( &(ostringstream() << gen) )->str();
		ofstream pareto;
		pareto.open(("PARETO/"+test + "_gen_" + g +".dat").c_str());
		for(int pop=0;pop<this->sizeEP;pop++){
			for(int k=0;k<numbObj;k++){
				pareto << this->EP[pop][k] << " ";
//				cout << this->EP[pop][k] << " ";
			}
			pareto << endl;
//			cout << endl;
		}
		pareto.close();

}
