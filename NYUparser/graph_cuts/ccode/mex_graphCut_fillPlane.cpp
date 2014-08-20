#include <math.h>
#include <stdlib.h>
#include <queue>

#include "mex.h"
#include "graph.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
 	
	if (nrhs != 3)
		mexErrMsgTxt("3 input arguments are required");
	if (nlhs>4)
		mexErrMsgTxt("Too many outputs");

	const int *sizeDimUnary = mxGetDimensions(prhs[0]);

	////////// Get the dimensions of the input image  /////////////////////////////
	int nrows = sizeDimUnary[0]; 
	int ncols = sizeDimUnary[1]; 
	const int nodeNum = nrows*ncols; // number of nodes

	/////////// Get the pointers to the input Data ///////////////////////////////
	double *UnaryPtr = mxGetPr(prhs[0]);
	
    double *PairPtr = mxGetPr(prhs[1]);
	
	double *NeighSystemPtr = mxGetPr(prhs[2]);
	int NeighSystem = (int) *NeighSystemPtr;

	
	/////////////////////// create graph //////////////////////
	typedef Graph<double,double,double> GraphType;
	GraphType *g = new GraphType( nodeNum,nodeNum*4 + nodeNum*4*NeighSystem); 
	g -> add_node(nodeNum);

	///////////// Add data term ////////////////////////////
	for (int i = 0; i < nodeNum; i++){
		g -> add_tweights(i, UnaryPtr[i + nodeNum] ,UnaryPtr[i]);
	}

	///////////// Add edge term /////////////////////////////
	double term01; 
    double term10;
	assert(lambda1 >= 0);
	assert(lambda2 >= 0);
	for (int i = 0; i < nrows; i++){
		for (int j = 0; j < ncols; j++){
			// vertical edges (i,j) <-> (i+1,j)
			if( i < nrows-1){				
                term01 = PairPtr[i+j*nrows];
                term10 = PairPtr[(i+1)+j*nrows];                        
				assert(term01 >= 0);  
                assert(term10 >= 0);
                g -> add_edge(i+j*nrows,(i+1)+j*nrows,term01,term10);
			}
			//horizontal edges (i,j) <-> (i,j+1)
			if( j<ncols-1){
                term01 = PairPtr[i+j*nrows];
                term10 = PairPtr[i+(j+1)*nrows];                        
				assert(term01 >= 0);  
                assert(term10 >= 0);
                g -> add_edge(i+j*nrows,i+(j+1)*nrows,term01,term10);				
			}	
			if(NeighSystem){
				//diagonal edges (i,j) <-> (i+1,j+1) 
				if( i<nrows-1 && j<ncols-1 ){
					term01 = PairPtr[i+j*nrows];
                    term10 = PairPtr[(i+1)+(j+1)*nrows];                        
                    assert(term01 >= 0);  
                    assert(term10 >= 0);
                    g -> add_edge(i+j*nrows,(i+1)+(j+1)*nrows,term01/sqrt(2.0),term10/sqrt(2.0));                    
				}
				//diagonal edges (i,j) <-> (i-1,j+1)
				if ( i > 0 && j<ncols-1){
					term01 = PairPtr[i+j*nrows];
                    term10 = PairPtr[(i-1)+(j+1)*nrows];                        
                    assert(term01 >= 0);  
                    assert(term10 >= 0);
                    g -> add_edge(i+j*nrows,(i-1)+(j+1)*nrows,term01/sqrt(2.0),term10/sqrt(2.0));
				}
			}
		}
	}

	double Emin = g -> maxflow();
                   
	//////////////////////////////// ASSIGN OUTPUT //////////////////////////
	plhs[0] = mxCreateDoubleMatrix(nrows, ncols,  mxREAL);
	double *labelOutPtr = mxGetPr(plhs[0]); 

	for (int j=0; j<nodeNum; j++)
	{
		labelOutPtr[j] = (double) g->what_segment(j);
	}


	plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
	double *energyPtr = mxGetPr(plhs[1]); 
	energyPtr[0] = Emin;
	
	if (nlhs==4 && nrhs < 7){
		plhs[2] = mxCreateDoubleMatrix(nrows, ncols,  mxREAL);
		double *minMarginal0Ptr = mxGetPr(plhs[2]); 

		plhs[3] = mxCreateDoubleMatrix(nrows, ncols,  mxREAL);
		double *minMarginal1Ptr = mxGetPr(plhs[3]); 

		double hard_constraint = ((((unsigned)-1)/2));
		for (int j=0; j<nodeNum; j++)
		{
			if (labelOutPtr[j] == 1)
			{
				g -> mark_node(j);
				g -> add_tweights(j, hard_constraint,0);
				
				minMarginal0Ptr[j] = g -> maxflow(1);
				minMarginal1Ptr[j] = Emin; 
				
				g -> add_tweights(j, -hard_constraint,0);
				g -> mark_node(j);
			}
			else
			{
				g -> mark_node(j);
				g -> add_tweights(j,0, hard_constraint);
				
				minMarginal1Ptr[j] = g -> maxflow(1);
				minMarginal0Ptr[j] = Emin; 
				
				g -> mark_node(j);
				g -> add_tweights(j,0, -hard_constraint);		
			}
		}       
        
	}
	
    if (nlhs==4 && nrhs == 7 ){
    	double *CoordinatesPtr = mxGetPr(prhs[6]);// assumes is a pointer to a coordinates_length * 2 matrix
        const int *sizeDimCoordinates = mxGetDimensions(prhs[6]);
        int coordinates_length = sizeDimCoordinates[0]; 
        
        plhs[2] = mxCreateDoubleMatrix(nrows, ncols,  mxREAL);
		double *minMarginal0Ptr = mxGetPr(plhs[2]); 

		plhs[3] = mxCreateDoubleMatrix(nrows, ncols,  mxREAL);
		double *minMarginal1Ptr = mxGetPr(plhs[3]); 
        
        int ind_i, ind_j;

        double hard_constraint = ((((unsigned)-1)/2));
        
        //0 min-marginal
        for (int i = 0; i < nrows; i++)
        {
            for (int j = 0; j < ncols; j++)
            {
                for (int k=0; k<coordinates_length; k++)
                {               
                    ind_i = i + (int) CoordinatesPtr[k];
                    ind_j = j + (int) CoordinatesPtr[k+coordinates_length];
                    if(ind_i < 0 || ind_i >=nrows || ind_j <0 || ind_j >= ncols)
                        continue;
                    g -> mark_node(ind_i + ind_j*nrows);
                    g -> add_tweights(ind_i + ind_j*nrows, hard_constraint,0);
                }    
                minMarginal0Ptr[i + j*nrows] = g -> maxflow(1);

                for (int k=0; k<coordinates_length; k++)
                {               
                    ind_i = i + (int) CoordinatesPtr[k];
                    ind_j = j + (int) CoordinatesPtr[k+coordinates_length];
                    if(ind_i < 0 || ind_i >=nrows || ind_j <0 || ind_j >= ncols)
                        continue;
                    g -> mark_node(ind_i + ind_j*nrows);
                    g -> add_tweights(ind_i + ind_j*nrows, -hard_constraint,0);
                } 
            }
        }
            
        //1 min-marginal
        for (int i = 0; i < nrows; i++)
        {
            for (int j = 0; j < ncols; j++)
            {
                for (int k=0; k<coordinates_length; k++)
                {               
                    ind_i = i + (int) CoordinatesPtr[k];
                    ind_j = j + (int) CoordinatesPtr[k+coordinates_length];
                    if(ind_i < 0 || ind_i >=nrows || ind_j <0 || ind_j >= ncols)
                        continue;
                    g -> mark_node(ind_i + ind_j*nrows);
                    g -> add_tweights(ind_i + ind_j*nrows, 0, hard_constraint);
                }    
                minMarginal1Ptr[i + j*nrows] = g -> maxflow(1);

                for (int k=0; k<coordinates_length; k++)
                {               
                    ind_i = i + (int) CoordinatesPtr[k];
                    ind_j = j + (int) CoordinatesPtr[k+coordinates_length];
                    if(ind_i < 0 || ind_i >=nrows || ind_j <0 || ind_j >= ncols)
                        continue;
                    g -> mark_node(ind_i + ind_j*nrows);
                    g -> add_tweights(ind_i + ind_j*nrows, 0, -hard_constraint);
                } 
            }
        }
    }
  
	delete g;

}

