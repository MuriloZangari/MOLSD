/*
 *  Tools.cpp
 *  RankingEDAsCEC
 *
 *  Created by Josu Ceberio Uribe on 11/21/11.
 *  Copyright 2011 University of the Basque Country. All rights reserved.
 *
 */

#include "NTools.h"
#include <math.h>
#include <sys/time.h>
#include "NVariables.h"

/*
 * It determines if the given int sequecen if it is indeed a permutation or not.
 */
bool isPermutation(int * permutation, int size)
{
	int flags[size];
	//int * flags=new int[size];
	for (int i=0;i<size;i++) flags[i]=1;
	
	for (int i=0;i<size;i++)
	{
		int value=permutation[i];
		flags[value]=0;
	}
	
	int result,sum=0;
	for(int i=0;i<size;i++)
		sum+=flags[i];
	if (sum==0) result=true;
	else result=false;
	//delete [] flags;
	return result;
}

/*
 * Generates a random permutation of size 'n' in the given array.
 */
void GenerateRandomPermutation(int * permutation, int n)
{
	for (int i = 0; i < n; ++i)
	{
		int j = rand() % (i + 1);
		permutation[i] = permutation[j];
		permutation[j] = i;
	}
}

/*
 * Determines if a given string contains a certain substring.
 */
bool strContains(const string inputStr, const string searchStr)
{
	size_t contains;
	
	contains = inputStr.find(searchStr);
	
	if(contains != string::npos)
		return true;
	else
		return false;
}

/*
 * Prints in standard output 'length' integer elements of a given array.
 */
void PrintArray(int* array, int length, string text)
{
	cout<<text;
	for (int i=0;i<length;i++){
		cout<<array[i]<<" ";
	}
	cout<<" "<<endl;
}


/*
 * Prints in standard output 'length' long double elements of a given array.
 */
void PrintArray(long double* array, int length, string text)
{
	cout<<text;
	for (int i=0;i<length;i++){
		cout<<array[i]<<" ";
	}
	cout<<" "<<endl;
}

/*
 * Prints in the standard output the 'size' elements of the given array..
 */
void PrintArray(int * array, int size)
{
	for (int i=0;i<size;i++){
		cout<<array[i]<<" ";
	}
	cout<<" "<<endl;
}


/*
 * Prints in standard output 'length' double elements of a given array.
 */
void PrintArray(float* array, int length, string text)
{
	cout<<text<<": ";
	for (int i=0;i<length;i++){
        printf(" %3.4f,",array[i]);
	}
	printf("\n ");

}

/*
 * Prints in standard output 'length' double elements of a given array.
 */
void PrintArray(double* array, int length, string text)
{
    int i;
	cout<<text;
	for (i=0;i<length;i++)
		printf(" %3.4f ",array[i]);
	printf("\n ");
}

/*
 * Prints in the standard output given matrix.
 */
void PrintMatrix(int** matrix, int length, int length2, string text)
{
    int i,j;
	cout<<text;
	for (i=0;i<length;i++)
	{
		cout<<""<<endl;
		for (j=0;j<length2;j++)
		{
			cout<<matrix[i][j]<<" ";
		}
	}
	cout<<" "<<endl;
}

void PrintMatrix(float** matrix, int length, int length2, string text)
{
    int i,j;
	cout<<text;
	for (i=0;i<length;i++)
	{
		cout<<""<<endl;
		for (j=0;j<length2;j++)
		{
			cout<<matrix[i][j]<<" ";
		}
	}
	cout<<" "<<endl;
}


/*
 * Prints in standard output 'lengthxlength' double elements of a given matrix.
 */
void PrintMatrixDouble(double** matrix, int length, int length2, string text)
{
    int i,j;
	cout<<text<<endl;
	for (i=0;i<length;i++)
	{
		for (j=0;j<length2;j++)
			printf("%3.9f, ",matrix[i][j]);
		printf("\n");
	}
}

/*
 * Calculates the tau Kendall distance between 2 permutations.
 */
int Kendall(int* permutationA, int*permutationB, int size)
{
    int i,dist;
    int * v,*composition,*invertedB,*m_aux;
    
    dist=0;
    v=new int[size-1];
    composition=new int[size];
    invertedB=new int[size];
    m_aux=new int[size];
    
    Invert(permutationB, size, invertedB);
    Compose(permutationA,invertedB,composition,size);
    vVector_Fast(v,composition,size,m_aux);
    
    for (i = 0; i < size-1; i++)
        dist += v[i];
    
    delete [] composition;
    delete [] invertedB;
    delete [] v;
    delete [] m_aux;
    
    return dist;
}

/*
 * Calculates the tau Kendall distance between 2 permutations.
 */
int Kendall(int* permutationA, int*permutationB, int size, int * m_aux)
{
    int i,dist;
    int * v,*composition,*invertedB;
    
    dist=0;
    v=new int[size-1];
    composition=new int[size];
    invertedB=new int[size];
    
    Invert(permutationB, size, invertedB);
    Compose(permutationA,invertedB,composition,size);
    vVector_Fast(v,composition,size,m_aux);

    for (i = 0; i < size-1; i++)
        dist += v[i];
    
    delete [] composition;
    delete [] invertedB;
    delete [] v;
    
    return dist;
}

/*
 * Calculates the Kendall tau distance between 2 permutations.
 * Auxiliary parameters are used for multiple continuous executions.
 */
int Kendall(int* permutationA, int*permutationB, int size, int * m_aux, int * invertedB, int *  composition, int * v)
{
    int i,dist;
    dist=0;
    Invert(permutationB, size, invertedB);
    Compose(permutationA,invertedB,composition,size);
    vVector_Fast(v,composition,size,m_aux);
    for (i = 0; i < size-1; i++)
        dist += v[i];

    return dist;
}

/*
 * Calculates the Cayley distance between 2 permutations.
 */
int Cayley(int * permutationA, int * permutationB, int size)
{
    int * invertedB= new int[size];
    int * composition= new int[size];
    int * elemsToCycles= new int[size];
    int * maxPosInCycle= new int[size];
    int * freeCycle= new int[size];

    Invert(permutationB, size, invertedB);
    Compose(permutationA,invertedB,composition,size);
    
    int index,cycle,distance;
    
    for(int i=0;i<size;i++)
    {
        elemsToCycles[i]=-1;
        maxPosInCycle[i]=-1;
        freeCycle[i]=1;
    }
    
    while((index=NextUnasignedElem(elemsToCycles,size))!=-1)
    {
        cycle=FindNewCycle(freeCycle,size);
        freeCycle[cycle]=0;
        do
        {
            elemsToCycles[index]=cycle;
            index = composition[index];//para permus de 1..n =>index = sigma[index]-1;
        }
        while(elemsToCycles[index] == -1);
    }
    distance=size-FindNewCycle(freeCycle,size);
    
    delete [] invertedB;
    delete [] composition;
    delete [] elemsToCycles;
    delete [] maxPosInCycle;
    delete [] freeCycle;
    
    return distance;
}

/*
 * Calculates the Cayley distance between 2 permutations.
 */
int Cayley(int * permutationA, int * permutationB, int size,int * invertedB, int *  composition, int * elemsToCycles, int * maxPosInCycle, int * freeCycle)
{
    
    Invert(permutationB, size, invertedB);
    Compose(permutationA,invertedB,composition,size);
    
    int index,cycle,distance;
    
    for(int i=0;i<size;i++)
    {
        elemsToCycles[i]=-1;
        maxPosInCycle[i]=-1;
        freeCycle[i]=1;
    }
    
    while((index=NextUnasignedElem(elemsToCycles,size))!=-1)
    {
        cycle=FindNewCycle(freeCycle,size);
        freeCycle[cycle]=0;
        do
        {
            elemsToCycles[index]=cycle;
            index = composition[index];//para permus de 1..n =>index = sigma[index]-1;
        }
        while(elemsToCycles[index] == -1);
    }
    distance=size-FindNewCycle(freeCycle,size);
    
    return distance;
}


int FindNewCycle(int * freeCycle, int size){
    int i;
    
    
    for(i=0;i<size;i++)
        if(freeCycle[i])
            return i;
    return size;
}

int NextUnasignedElem(int * elemsToCycles, int size)
{
    int i;
    for(i=0;i<size;i++)
        if(elemsToCycles[i]==-1)
            return i;
    return -1;
}

/*
 * Calculates the Ulam distance between 2 permutations.
 */
int Ulam(int * permutationA, int * permutationB, int size){

    int * composition=new int[size];
    int * invertedB=new int[size];
    
    Invert(permutationB, size, invertedB);
    Compose(permutationA,invertedB,composition,size);

    int dist=(size- getLISLength(composition,size));
    
    delete [] composition;
    delete [] invertedB;
    
    return dist;
}

/*
 * Calculates the length of the longest increasing subsequence in the given array of ints.
 */
int getLISLength(int*sigma,int size){
   
    // O(n log k)
    
    int i;
    vector<int> vc(1,sigma[0]);
    vector<int>::iterator vk;
    
    for (i=1;i<size;i++){
        for (vk=vc.begin(); vk != vc.end(); vk++)
            if (*vk>=sigma[i]) break;
        if (vk==vc.end())
            vc.push_back(sigma[i]);
        else *vk=sigma[i];
    }
    
    return (int)vc.size();
    
}
/*
 * Implements the compose of 2 permutations of size n.
 */
void Compose(int*s1, int*s2, int*res, int n)
{
    int i;
    for(i=0;i<n;i++)
        res[i]=s1[s2[i]];
}


/*
 * Calculates V_j-s vector.
 */
void vVector(int*v, int*permutation, int n)
{
    
    int i,j;
    for(i=0;i<n-1;i++)
        v[i]=0;
    
    for (i = n-2; i >= 0; i--)
        for (j = i+1; j < n; j++)
            if(permutation[i] > permutation[j])
                v[i]++;
}

/*
 *  Optimized version proposed by Leti for the calculation of the V_j-s vector.
 */
void vVector_Fast(int*v, int*permutation, int n, int * m_aux)
{

    int i,j, index;

    for(i=0;i<n-1;i++){

        v[i]=0;
        m_aux[i]=0;
    }

    m_aux[n-1]=0;
    for (j=0; j<n-1; j++){
        index=permutation[j];
        v[j]=index-m_aux[index];
        for (i=index; i<n; i++)
            m_aux[i]++;
    }


}

void vVector_Fast2(int*v, int*permutation, int n)
{
	int * m_aux;
	m_aux = new int[n];
    int i,j, index;

    for(i=0;i<n-1;i++){

        v[i]=0;
        m_aux[i]=0;
    }

    m_aux[n-1]=0;
    for (j=0; j<n-1; j++){
        index=permutation[j];
        v[j]=index-m_aux[index];
        for (i=index; i<n; i++)
            m_aux[i]++;
    }

    delete [] m_aux;
}

/*
 * Inverts a permutation.
 */
void Invert(int*permu, int n, int* inverted)
{
    int i;
    for(i=0; i<n; i++)
        inverted[permu[i]]=i;
}

/*
 * Inverts a permutation.
 */
void NoInvert(int*permu, int n, int* inverted)
{
    int i;
    for(i=0; i<n; i++)
        inverted[i]=permu[i];
}

/*
 * Applies the random keys sorting strategy to the vector of doubles
 */
void RandomKeys( int * a, float * criteriaValues, int size)
{
	bool * fixedValues= new bool[size];
	float criteria, min;
	int i, j;
	for (i=0;i<size;i++){
		fixedValues[i]=false;
		a[i]=0;
	}
	int minPos=0;
	for (i=0;i<size;i++)
	{
		min=MAX_INTEGER;
		for (j=0;j<size;j++)
		{
			criteria=criteriaValues[j];
			if (!fixedValues[j] && min>criteria )
			{
				min=criteria;
				minPos=j;
			}
		}
        
		fixedValues[minPos]=true;
		//a[i]=minPos;// modification por el asunto ordering /ranking
    	a[minPos]=i;// original.
	}
	delete [] fixedValues;
}



/*
 * Calculates the factorial of a solution.
 */
long double factorial(int val) {
    if(val <= 0) return 1;
    //long  N, b, c, p; // use int for fast calculation and small range of calculation..
    long   b, c;
    long double p, N;
    N=(long double)val;
    c = (long)N - 1;
    p = 1;
    while (c > 0) {
        p = 0;
        b = c;
        while (b > 0) {
            if (b & 1) {
                p += N; // p = p + N;
            }
            // if you would like to use double choose the alternative forms instead shifts
            // the code is fast even!
            // you can use the same tips on double or 64 bit int etc.... but you must... ;-)
            //b >>= 1; // b/=2; (b = b / 2;) ( b >> 1; a.s.r. is more efficent for int or long..!)
            b/=2;
            //N <<= 1; // N += N; N = N + N; N = N * 2; (N <<=1; a.s.l. is more efficent for int or long..!)
            N += N;
        } // end of: while(b>0)
        N = p;
        c--; // c = c - 1;
    } // end of: while(c > 0)
    //printf("[%d] is the factorial! \n", p);
    return p;
}

/*
 * This method applies a swap of the given i,j positions in the array.
 */
void Swap(int * array, int i, int j)
{
	int aux=array[i];
	array[i]=array[j];
	array[j]=aux;
}

/*
 * PBI scalarizing functions methods
 */

double Norma(vector<double> &nambda){
	double sum=0.0;
	for(int k=0;k<numbObj;k++){
		sum = sum + nambda[k]*nambda[k];
	}
	return(sqrt(sum));

}

double innerProduct(vector <double>&vec1, vector <double>&vec2){
	double sum =0.0;
	for(int i=0; i<numbObj; i++){
		sum+= vec1[i]*vec2[i];
	}
	return sum;
}

int CalculateDistanceAndX2(int * sigma, int *x){
    //also updates the x vector if it isnot null
	bool * m_visited;
	m_visited = new bool[numbVar];
    if(x!=NULL)for (int i = 0 ; i < numbVar; i ++ )x[ i ] =1;

    int num_cycles=0, num_visited=0, item= 0;

    for (int i = 0 ; i < numbVar; i ++ )
        m_visited[ i ] =false;
    while(num_visited<numbVar){
        item=num_cycles;
        while(m_visited[item])item++;
        num_cycles++;
        int maxItemInCycle= 0;
        do{
            if(item>maxItemInCycle)maxItemInCycle=item;
            m_visited[item] =true;
            num_visited++;
            item=sigma[item];
        }while(!m_visited[item]);
        if(x!=NULL)x[maxItemInCycle] = 0;
    }
    delete [] m_visited;
    return (numbVar-num_cycles);
}

/*
 * This method moves the value in position i to the position j.
 */
void InsertAt(int * array, int i, int j, int n)
{
    if (i!=j)
    {
        int res[n];
        int val=array[i];
        if (i<j)
        {
            memcpy(res,array,sizeof(int)*i);

            for (int k=i+1;k<=j;k++)
                res[k-1]=array[k];

            res[j]=val;

            for (int k=j+1;k<n;k++)
                res[k]=array[k];
        }
        else if (i>j)
        {
            memcpy(res,array,sizeof(int)*j);

            res[j]=val;

            for (int k=j;k<i;k++)
                res[k+1]=array[k];

            for (int k=i+1;k<n;k++)
                res[k]=array[k];
        }
        memcpy(array,res,sizeof(int)*n);
    }
}
/*
 * Returns the min value's position in the given array.
 */
int MinPosition(int*array, int n)
{
    int min_value=intMax;
    int min_value_position=0;
    for(int i=0;i<n;i++)
    {
        if (array[i]<min_value)
        {
            min_value=array[i];
            min_value_position=i;
        }
    }
    return min_value_position;
}
/*
 * Sorts the doubles array in the ascending sequence.
 */
void AscendingSequence(int * array, int* indexes, int* resultingSequence, int size)
{
    int * aux= new int[size];
    //memcpy(aux, array, sizeof(int)*size);
    for (int i=0;i<size;i++)
    {
        int pos=MinPosition(array, size);
        resultingSequence[i]=indexes[pos];
        aux[i]=array[pos];
        array[pos]=intMax;

    }
    memcpy(array,aux, sizeof(int)*size);
    delete[] aux;
}
/*
 * Removes a position in a given array of ints and swap the solutions on the left.
 */
void RemoveAt(int*array,int positionIndex, int n)
{
    for (int i=positionIndex;i<n-1;i++)
    {
        array[i]=array[i+1];
    }
    array[n-1]=0;
}

void Shake_Insert(int * genes,  int shake_power)//, int orbit_size)
{
 //   int * inverted= new int[numbVar];
 //    Invert(genes, numbVar, inverted);
    int i,j;

    int min_range,max_range,val;
    int half_orbit=5;
    for (int iter=0;iter<shake_power;iter++)
    {
        //permute randomly genes in position i and j to scape the stackeness in both neighborhood.
        i = rand() % numbVar;
        //j = rand() % m_size;

        min_range=MAX(i-half_orbit,0);
        max_range=MIN(i+half_orbit,numbVar-1);
        val = rand() % (max_range-min_range);
        j= min_range+val;

        InsertAt(genes,i,j,numbVar);
        //PrintArray(genes,numbVar,"In ");
    }

//    Invert(inverted, numbVar,genes);
//    PrintArray(genes,numbVar,"test ");
//   m_value=MIN_LONG_INTEGER;
//   m_generations=0;
//   m_better_samples=0;
//    delete[] inverted;
}

void Shake_Swap(int * genes, int shake_power)
{
    //shake_power= rand() % shake_power;

    int i,j,aux;
    for (int iter=0;iter<shake_power;iter++)
    {
        //permute randomly genes in position i and j to scape the stackeness in both neighborhood.
        i = rand() % numbVar;
        j = rand() % numbVar;
        aux=genes[j];
        genes[j]=genes[i];
        genes[i]=aux;
        //cout<<individual<<endl;
    }
    //individual->SetGenes(genes);
}
