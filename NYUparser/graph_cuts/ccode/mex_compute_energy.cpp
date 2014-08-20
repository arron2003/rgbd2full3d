#include <math.h>
#include <stdlib.h>
#include <queue>

#include "mex.h"
#include "graph.h"



double getEdgeTerm(int i, int j, int nodeNum, double* ImagePtr, double EdgeBeta)
{
        double dist1 = pow( ((double) ImagePtr[i]-ImagePtr[j]), 2);
        dist1 += pow( ((double) ImagePtr[i + nodeNum]-ImagePtr[j + nodeNum]), 2);
        dist1 += pow( ((double) ImagePtr[i + nodeNum*2]-ImagePtr[j + nodeNum*2]), 2);

        return exp(-EdgeBeta*dist1);                 
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
 	
	if (nrhs > 9)
		mexErrMsgTxt("6/7 input arguments are required");
	if (nlhs>4)
		mexErrMsgTxt("Too many outputs");

	const int *sizeDimUnary = mxGetDimensions(prhs[0]);

	////////// Get the dimensions of the input image  /////////////////////////////
	int nrows = sizeDimUnary[0]; 
	int ncols = sizeDimUnary[1]; 
	const int nodeNum = nrows*ncols; // number of nodes

	/////////// Get the pointers to the input Data ///////////////////////////////
	double *UnaryPtr = mxGetPr(prhs[0]);
	
	double *lambda1Ptr = mxGetPr(prhs[1]);
	double lambda1 = (double) *lambda1Ptr;
  
	double *lambda2Ptr = mxGetPr(prhs[2]);
	double lambda2 = (double) *lambda2Ptr;
	
	double *NeighSystemPtr = mxGetPr(prhs[3]);
	int NeighSystem = (int) *NeighSystemPtr;

	double *ImagePtr = mxGetPr(prhs[4]);
	
	double *EdgeBetaPtr = mxGetPr(prhs[5]);
	double EdgeBeta = (double) *EdgeBetaPtr;

    double *activeNodes = mxGetPr(prhs[6]);
    
    double *numActiveNodesPtr = mxGetPr(prhs[7]);
	int numActiveNodes = (int) *numActiveNodesPtr;
 
    double constEnergy=0;
	double *tempUnary_0 = new double[nodeNum];
    double *tempUnary_1 = new double[nodeNum];
    int *nodeIndex = new int[nodeNum];
    for (int i=0;i<nodeNum;i++)
    {
        tempUnary_0[i] = 0;
        tempUnary_1[i] = 0;
    }
    
    /////////////////////// create graph //////////////////////
	typedef Graph<double,double,double> GraphType;
	GraphType *g = new GraphType(numActiveNodes,numActiveNodes*4*NeighSystem); 
	g -> add_node(numActiveNodes);

    int nodeCount=0;
	for (int i = 0; i < nodeNum; i++)
    {
		if (activeNodes[i]==1) 
        {
            nodeIndex[i] = nodeCount;
            nodeCount++;
        }
        else nodeIndex[i] = -1;
	}

    	//////////// Add edge term /////////////////////////////
	double term01;
	assert(lambda1 >= 0);
	assert(lambda2 >= 0);
    
    int temp_f, temp_t;
    
	for (int i = 0; i < nrows; i++){
		for (int j = 0; j < ncols; j++){
           
			// vertical edges (i,j) <-> (i+1,j)
			if( i < nrows-1){
				term01 = lambda1+lambda2*getEdgeTerm(i+j*nrows, (i+1)+j*nrows, nodeNum, ImagePtr, EdgeBeta);
				assert(term01 >= 0);
                
                temp_f = i+j*nrows;
                temp_t = (i+1)+j*nrows;
                
                if (activeNodes[temp_f]==1 && activeNodes[temp_t]==1)
                    g -> add_edge(nodeIndex[temp_f],nodeIndex[temp_t],term01,term01);
                else if ((activeNodes[temp_f]!=1) && (activeNodes[temp_t]!=1) && (activeNodes[temp_f]!=activeNodes[temp_f]))
                    constEnergy += term01;
                else if (activeNodes[temp_f]!=1)
                {
                    if (activeNodes[temp_f]==0) tempUnary_1[temp_t] += term01;
                    else tempUnary_0[temp_t] += term01;
                }
                else 
                {
                    if (activeNodes[temp_t]==0) tempUnary_1[temp_f] += term01;
                    else tempUnary_0[temp_f] += term01;
                }
			}
			//horizontal edges (i,j) <-> (i,j+1)
			if( j<ncols-1){
				term01 = lambda1+lambda2*getEdgeTerm(i+j*nrows, i+(j+1)*nrows, nodeNum, ImagePtr, EdgeBeta);
				assert(term01 >= 0);
    
                temp_f = i+j*nrows;
                temp_t = i+ (j+1)*nrows;
                
                if (activeNodes[temp_f]==1 && activeNodes[temp_t]==1)
                    g -> add_edge(nodeIndex[temp_f],nodeIndex[temp_t],term01,term01);
                else if ((activeNodes[temp_f]!=1) && (activeNodes[temp_t]!=1) && (activeNodes[temp_f]!=activeNodes[temp_f]))
                    constEnergy += term01;
                else if (activeNodes[temp_f]!=1)
                {
                    if (activeNodes[temp_f]==0) tempUnary_1[temp_t] += term01;
                    else tempUnary_0[temp_t] += term01;
                }
                else 
                {
                    if (activeNodes[temp_t]==0) tempUnary_1[temp_f] += term01;
                    else tempUnary_0[temp_f] += term01;
                }
			}	
			if(NeighSystem){
				//diagonal edges (i,j) <-> (i+1,j+1) 
				if( i<nrows-1 && j<ncols-1 ){
					term01 = (lambda1+lambda2*getEdgeTerm(i+j*nrows,(i+1)+(j+1)*nrows, nodeNum, ImagePtr, EdgeBeta))/sqrt(2.0);

                    assert(term01 >= 0);

                    temp_f = i+j*nrows;
                    temp_t = (i+1) + (j+1)*nrows;

                    if (activeNodes[temp_f]==1 && activeNodes[temp_t]==1)
                        g -> add_edge(nodeIndex[temp_f],nodeIndex[temp_t],term01,term01);
                    else if ((activeNodes[temp_f]!=1) && (activeNodes[temp_t]!=1) && (activeNodes[temp_f]!=activeNodes[temp_f]))
                        constEnergy += term01;
                    else if (activeNodes[temp_f]!=1)
                    {
                        if (activeNodes[temp_f]==0) tempUnary_1[temp_t] += term01;
                        else tempUnary_0[temp_t] += term01;
                    }
                    else 
                    {
                        if (activeNodes[temp_t]==0) tempUnary_1[temp_f] += term01;
                        else tempUnary_0[temp_f] += term01;
                    }
                }
				//diagonal edges (i,j) <-> (i-1,j+1)
				if ( i > 0 && j<ncols-1){
					term01 = (lambda1+lambda2*getEdgeTerm(i+j*nrows,(i-1)+(j+1)*nrows, nodeNum, ImagePtr, EdgeBeta))/sqrt(2.0);
					assert(term01 >= 0);

                    temp_f = i+j*nrows;
                    temp_t = (i-1) + (j+1)*nrows;

                    if (activeNodes[temp_f]==1 && activeNodes[temp_t]==1)
                        g -> add_edge(nodeIndex[temp_f],nodeIndex[temp_t],term01,term01);
                    else if ((activeNodes[temp_f]!=1) && (activeNodes[temp_t]!=1) && (activeNodes[temp_f]!=activeNodes[temp_f]))
                        constEnergy += term01;
                    else if (activeNodes[temp_f]!=1)
                    {
                        if (activeNodes[temp_f]==0) tempUnary_1[temp_t] += term01;
                        else tempUnary_0[temp_t] += term01;
                    }
                    else 
                    {
                        if (activeNodes[temp_t]==0) tempUnary_1[temp_f] += term01;
                        else tempUnary_0[temp_f] += term01;
                    }
                }
			}
		}
	}

    ///////////// Add data term ////////////////////////////
	for (int i = 0; i < nodeNum; i++){
		if (activeNodes[i]==1)
//            g->add_tweights(nodeIndex[i],tempUnary_1[i],tempUnary_0[i]);
            g->add_tweights(nodeIndex[i], UnaryPtr[i+nodeNum]+tempUnary_1[i],UnaryPtr[i]+tempUnary_0[i]);
        else if (activeNodes[i] == 0)
            constEnergy += UnaryPtr[i+nodeNum]+tempUnary_1[i];
        else if (activeNodes[i] == 2)
            constEnergy += UnaryPtr[i]+tempUnary_0[i];
	}

    double Emin = g -> maxflow();
                   
	//////////////////////////////// ASSIGN OUTPUT //////////////////////////
	plhs[0] = mxCreateDoubleMatrix(nrows, ncols,  mxREAL);
	double *labelOutPtr = mxGetPr(plhs[0]); 

	for (int j=0; j<nodeNum; j++)
	{
		if (activeNodes[j]==1) 
            labelOutPtr[j] = (double) g->what_segment(nodeIndex[j]);
        else if (activeNodes[j]==0) 
            labelOutPtr[j] = (double) 0;
        else
            labelOutPtr[j] = (double) 1;
	}

	plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
	double *energyPtr = mxGetPr(plhs[1]); 
	energyPtr[0] = Emin;
	delete g;
}

