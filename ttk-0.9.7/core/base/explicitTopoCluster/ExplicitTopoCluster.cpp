#include                  <ExplicitTopoCluster.h>


using namespace std;
using namespace ttk;

ExplicitTopoCluster::ExplicitTopoCluster(){
  clear();
  /* LFU cache */
  // caches_.resize(1);
  // freqMaps_.resize(1);
  // minFrequency_.resize(1);
  // #ifdef TTK_ENABLE_OPENMP
  // caches_.resize(threadNumber_);
  // freqMaps_.resize(threadNumber_);
  // minFrequency_.resize(threadNumber_);
  // #endif

  /* LRU cache */
  caches_.resize(1);
  cacheMaps_.resize(1);
  #ifdef TTK_ENABLE_OPENMP
  caches_.resize(threadNumber_);
  cacheMaps_.resize(threadNumber_);
  #endif
}

ExplicitTopoCluster::~ExplicitTopoCluster(){
}

int ExplicitTopoCluster::clear(){
  vertexNumber_ = 0;
  cellNumber_ = 0;
  nodeNumber_ = 0;
  doublePrecision_ = false;

  {
    stringstream msg;
    msg << "[ExplicitTopoCluster] Triangulation cleared." << endl;
    dMsg(cout, msg.str(), detailedInfoMsg);
  }

  return AbstractTriangulation::clear();
}