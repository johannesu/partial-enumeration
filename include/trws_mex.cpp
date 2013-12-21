// Johannes Ulén, 2013
// MATLAB wrapper for Convergent Tree-reweighted Message Passing for Energy Minimization (TRW-S) by Vladimir Kolmogorov
// http://pub.ist.ac.at/~vnk/papers/TRW-S.html
//
// The wrapper is inspired by Olivier Woodfoord's imrender
// and uses it's modification to TRW-S to make it compile with MATLAB.
// http://www.robots.ox.ac.uk/~ojw/software.htm
#include "mex.h"
#include "utils/mexutils.h"
#include "utils/cppmatrix.h"
#include "trw-s/MRFEnergy.h"
#include <cstring>

// Catch errors
static void erfunc(char *err) {
	mexErrMsgTxt(err);
}

template<typename TYPE> 
static void solve_mrf(std::string type, 
					  int nlhs, 
					  mxArray *plhs[], 
					  int nrhs, 
					  const mxArray  *prhs[])
{  
  unsigned int curarg = 1; // type skipped

  // Pairwise costs are defined by pairwise ids and cost
  matrix<int>  patch_size(prhs[curarg++]);
  matrix<double> unary(prhs[curarg++]);
  matrix<unsigned int> pairwise_ids(prhs[curarg++]);
  matrix<double> pairwise_costs(prhs[curarg++]);
  
  // Extra data needed for the sparse solver
  matrix<int> labels(prhs[curarg++]);
  

  MexParams params(nrhs-curarg, prhs+curarg); //Structure to hold and parse additional parameters
 
  const int max_iter      = params.get<int>("max_iter", 1000);
  const double max_relgap = params.get<double>("max_relgap", 0);
  const bool verbose   = params.get<bool>("verbose", false);
  const bool automatic_ordering = params.get<bool>("automatic_ordering", false);
    
  ASSERT(pairwise_ids.M == 2);
  ASSERT(pairwise_ids.N == pairwise_costs.N);

  int n_labels = unary.M; // Number of labels
  int n_nodes = unary.N; // Number of unary terms
  int n_edges = pairwise_ids.N; // Number of pairwise terms

  if (verbose)
  {
    mexPrintf("Type: %s \n", type.c_str());
    mexPrintf("Maximum number of iterations %d \n", max_iter);
    mexPrintf("Maximum relative duality gap %g \n", max_relgap);
    mexPrintf("Number of labels: %d \n", n_labels);
    mexPrintf("Number of variables: %d \n", n_nodes);
    mexPrintf("Number of edges: %d \n", n_edges);
  }

  // Create MRF instance
  typedef typename TYPE::REAL REAL;
  MRFEnergy<TYPE> *mrf;

  // Allocate memory destructor of TYPE is not invoked.
  matrix<int> hor_source(n_labels);
  matrix<int> hor_sink(n_labels);
  matrix<int> vert_source(n_labels);
  matrix<int> vert_sink(n_labels);

  mrf =  new MRFEnergy<TYPE>( construct_globalsize( mrf,
                                                    n_labels,
                                                    patch_size(0),
                                                    &labels(0),
                                                    &hor_source(0),
                                                    &hor_sink(0),
                                                    &vert_source(0),
                                                    &vert_sink(0)
                                                    ),
                                                  erfunc);    

  typename MRFEnergy<TYPE>::NodeId * nodes =  new typename MRFEnergy<TYPE>::NodeId[n_nodes];

	// Add unary terms
  for (int u = 0; u < n_nodes; u++)
  {
    nodes[u] = mrf->AddNode(typename TYPE::LocalSize(n_labels), 
                            typename TYPE::NodeData(&unary(0,u)));
  }

  if (verbose)
    mexPrintf("Added all unary costs \n");

  for (int p = 0; p < n_edges; p++)
  {
    add_edge( mrf, 
              nodes[pairwise_ids(0,p)],
              nodes[pairwise_ids(1,p)],
              &pairwise_costs(0,p)
           );
  }

  if (verbose)
    mexPrintf("Added all pairwise costs \n");

  if (automatic_ordering)
    mrf->SetAutomaticOrdering();

  // Options struct
  typename MRFEnergy<TYPE>::Options mrf_options;
  mrf_options.m_iterMax = max_iter;
  mrf_options.m_relgapMax = max_relgap; 

  if (verbose)
    mrf_options.m_printMinIter = true;
  else
    mrf_options.m_printMinIter = false;

  // Solve

  if (verbose)
    mexPrintf("Starting the optimization. \n");

  REAL energy, lowerBound = 0;
  int iterations = mrf->Minimize_TRW_S(mrf_options, lowerBound, energy);

  if (verbose)
    mexPrintf("Optimization, done parsing output. \n");

  // Output
  matrix<double> labelling(n_nodes);

  // Read solution
  for (int u = 0; u < n_nodes; u++ )
    labelling(u) = mrf->GetSolution(nodes[u])+1;

  plhs[0] = labelling;
  plhs[1] = mxCreateDoubleScalar((double)energy);
  plhs[2] = mxCreateDoubleScalar((double)lowerBound);
  plhs[3] = mxCreateDoubleScalar((double)iterations);

  delete nodes;
  delete mrf;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,  const mxArray  *prhs[])
{
  ASSERT(nrhs == 6 || nrhs == 7);
  ASSERT(nlhs == 4 || nlhs == 5);
  
  char buffer[1024];
  if (mxGetString(prhs[0], buffer, 1024)) 
    throw std::runtime_error("Expected string parameter value");
  std::string type(buffer);

  if (!type.compare("general"))
  	solve_mrf<TypeGeneral>(type, nlhs, plhs, nrhs, prhs);
  else if (!type.compare("partial_enumeration"))
    solve_mrf<TypePartialEnumeration>(type, nlhs, plhs, nrhs, prhs);
  else
   throw std::runtime_error( "Unknown energy type" );
}

// Template specializations
//
// TypeGeneral
//
static inline TypeGeneral::GlobalSize 
construct_globalsize(MRFEnergy<TypeGeneral> *mrf,
                     int n_labels, 
                     int patch_size,
                     int * labels,
                     int * hor_source,
                     int * hor_sink,
                     int * vert_source,
                     int * vert_sink
                     )
{
  return TypeGeneral::GlobalSize(n_labels);
}

static inline void add_edge(MRFEnergy<TypeGeneral> *mrf, 
                            MRFEnergy<TypeGeneral>::NodeId node1, 
                            MRFEnergy<TypeGeneral>::NodeId node2, 
                            TypeGeneral::REAL *data)
{
  mrf->AddEdge(node1, node2, TypeGeneral::EdgeData(TypeGeneral::GENERAL, data));
}

//
//  typePatchConsSparseSorted
//

static inline TypePartialEnumeration::GlobalSize 
construct_globalsize(MRFEnergy<TypePartialEnumeration> *mrf,
                     int n_labels, 
                     int patch_size,
                     int * labels,
                     int * hor_source,
                     int * hor_sink,
                     int * vert_source,
                     int * vert_sink
                     )
{
  return TypePartialEnumeration::GlobalSize(n_labels, 
                                               patch_size,
                                               labels,
                                               hor_source,
                                               hor_sink,
                                               vert_source,
                                               vert_sink,
                                               true
                                            );
}

static inline void add_edge(MRFEnergy<TypePartialEnumeration> *mrf, 
                            MRFEnergy<TypePartialEnumeration>::NodeId node1, 
                            MRFEnergy<TypePartialEnumeration>::NodeId node2, 
                            TypePartialEnumeration::REAL *data)
{
  if (data[0] == 1)
    mrf->AddEdge(node1, 
                 node2, 
                 TypePartialEnumeration::EdgeData(TypePartialEnumeration::PatchConsSparseSortedVert));
  else
    mrf->AddEdge(node1, 
                 node2, 
                 TypePartialEnumeration::EdgeData(TypePartialEnumeration::PatchConsSparseSortedHor)); 
}