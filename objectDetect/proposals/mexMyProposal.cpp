#include "mex.h"
#include "int_tree.h"
#include <map>
#if !defined(MAX)
#define    MAX(A, B)    ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define    MIN(A, B)    ((A) < (B) ? (A) : (B))
#endif

double Unif01(){
  return ((double) rand() / (double) (RAND_MAX+1.0));
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
  // input 0 edge list (n x 2)
  // input 1 edge weights (n x 1)
  // input 2 list of region size (n x 1)
  // input 3 list of initial seeds (sx1)

  if (nrhs != 4) {
    mexErrMsgTxt("3 input argument required: edgeList, edgeWeight, sizeList, Seeds");
  } 
  if ( !mxIsInt16(prhs[0]) )
    mexErrMsgTxt("input 0 (edgeList must be INT16");
  if ( !mxIsDouble(prhs[1]) )
    mexErrMsgTxt("input 1 (edgeWeight must be DOUBLE");
  if ( !mxIsDouble(prhs[2]) )
    mexErrMsgTxt("input 2 (sizeList must be DOUBLE");
  if ( !mxIsInt16(prhs[3]) )
    mexErrMsgTxt("input 3 (seedList must be INT16");
  if ( mxGetN(prhs[0])!=2 )
    mexErrMsgTxt("Dim of input 1 must be n x 2");
  if ( mxGetM(prhs[0])!=mxGetM(prhs[1]) )
    mexErrMsgTxt("Dim of input 1 and input 2 mismatch");

  const size_t nEdge = mxGetM(prhs[0]);
  const size_t nProposals = mxGetM(prhs[3]);
  plhs[0]=mxCreateCellMatrix(nProposals, 1);

  int nSps = -1;
  std::vector< std::pair<short, short> > edgeList;
  std::vector<double> edgeWeight;
  std::vector<short> seedList;

  const short* tPtr = (short*) mxGetData(prhs[0]);
  const double* wPtr = (double*) mxGetData(prhs[1]);
  const double* rPtr = (double*) mxGetData(prhs[2]);
  const short* sPtr = (short*) mxGetData(prhs[3]);
  for (int i=0; i<nEdge; ++i) {
    nSps = MAX(nSps, tPtr[i]);
    nSps = MAX(nSps, tPtr[i+nEdge]);
	edgeList.push_back(std::pair<short, short>(tPtr[i], tPtr[i+nEdge]) );
	edgeWeight.push_back(wPtr[i]);
	//mexPrintf("%d %d %.3f\n", edgeList[i].first, edgeList[i].second, edgeWeight[i]);
  }
  for (int i=0; i<nProposals; ++i)
    seedList.push_back(sPtr[i]);
  //mexPrintf("nSps = %d\n", nSps);
  std::vector< std::vector< std::pair<short, double> > > neighList(nSps);
  for (int i=0; i<nEdge; ++i) {
    uint nodei=MAX(edgeList[i].second, edgeList[i].first)-1;
    uint nodej=MIN(edgeList[i].second, edgeList[i].first)-1;
    neighList[nodei].push_back( std::pair<short, double>(nodej, edgeWeight[i]) );
    neighList[nodej].push_back( std::pair<short, double>(nodei, edgeWeight[i]) );
  }
  /* -- Start to do the work --*/
  IntTree T(nSps);
  std::vector<std::vector<bool> > spGroups(nProposals, std::vector<bool>(nSps,0));

  uint nextSp=0, oSp=0;
  uint n=0, nSpsInGroup=0;
  double E0=0.0, Ea = 0.0, E=0.0, nextS=0.0;
  uint xmin=0, ymin=0, xmax=0, ymax=0;
  for( int k=0; k<nProposals; k++) {
    int nextSp = seedList[k]-1;
    assert(nextSp<nSps && "Bad seed");
    assert(k < spGroups.size());
    assert(nextSp < spGroups.at(k).size());
    spGroups.at(k).at(nextSp)=1;

    nSpsInGroup=1;

    if(nSpsInGroup==nSps)
      continue;
	//mexPrintf("initial seed = %d\n", nextSp);
    std::vector<std::pair<short, double> > neighs= neighList[nextSp];
    std::vector<IntTree::Node> ns;     
    for(n=0; n<neighs.size(); ++n){
      assert(n < neighs.size());
      double e=neighs.at(n).second;
      uint i=MAX(neighs.at(n).first,nextSp);
      uint j=MIN(neighs.at(n).first,nextSp);

      ns.push_back(IntTree::Node(e,i,j,std::vector<double>::iterator(NULL)));
      T.AddNode(ns.back());
	  //mexPrintf("Add edge = <%d, %d> -> %.3f\n", i, j, e);
    }

    assert(T.AreWConsistent());

    E0=Unif01();E0=.5*E0;
    Ea=1.0*Unif01(); Ea=Ea*Ea*Ea;
    E=rPtr[nextSp];
    while(1){
      assert(T.AreWConsistent());
      IntTree::Node nextNode(T.SampleNode(Unif01()));

      assert(k < spGroups.size());
      assert(nextNode.i() < spGroups.at(k).size());
      if(spGroups.at(k).at(nextNode.i())){
        nextSp=nextNode.j();
        oSp=nextNode.i();
      }else{
        assert(k < spGroups.size());
        assert(nextNode.j() < spGroups.at(k).size());
        assert(spGroups.at(k).at(nextNode.j()));
        nextSp=nextNode.i();
        oSp=nextNode.j();
      }
      nextS=nextNode.e();

      //E=sqrt(1-(1-double(nSpsInGroup)/nSps)*(1-double(nSpsInGroup)/nSps));
	  //E = double(nSpsInGroup)/nSps;
      //E = (double(nSpsInGroup)/nSps)*(double(nSpsInGroup)/nSps);
      E += rPtr[nextSp];
      assert(E>=0 && E<=2.0);

      if( E>Ea || nextS<E0 ){
	    //mexPrintf("Break! E=%.3f<E0=%.3f\n", E, E0);
        break;
      }else{
        assert(k < spGroups.size());
        assert(nextSp < spGroups.at(k).size());
        spGroups.at(k).at(nextSp)=1;
        nSpsInGroup++;

        if(nSpsInGroup==nSps)
          break;

        std::vector<std::pair<short, double> > N = neighList[nextSp];

	    //mexPrintf("Node osp = %d, sampled = %d, nextS=%.3f\n", oSp, nextSp, nextNode.e());
        for(n=0;n<N.size();n++){
          assert(n < N.size());
          uint j=N.at(n).first;
          assert(k < spGroups.size());
          assert(j < spGroups.at(k).size());
          if(spGroups.at(k).at(j)){
	        //mexPrintf("Removing node = %d %d\n", MAX(nextSp,j),MIN(nextSp,j));
            T.RemoveNode(MAX(nextSp,j),MIN(nextSp,j));

            assert(T.AreWConsistent());
          }else{
            assert(n < N.size());
            T.AddNode(IntTree::Node(N.at(n).second,MAX(nextSp,j),MIN(nextSp,j),std::vector<double>::iterator(NULL)));
            assert(T.AreWConsistent());
          }
        }
      }
    }
    /*-- Output the region proposals --*/
	std::vector<short> members;
	for(uint s=0;s<nSps;s++)
      if(spGroups.at(k).at(s))
	    members.push_back(s);
    mxArray* matMembers = mxCreateNumericMatrix(members.size() , 1, mxINT16_CLASS, mxREAL);
	mxSetCell(plhs[0], k, matMembers);
	short* memberList = (short*) mxGetPr(matMembers);
    for (int i=0; i<members.size(); i++)
	  memberList[i] = members[i];

    T.Reset();
  }
}
