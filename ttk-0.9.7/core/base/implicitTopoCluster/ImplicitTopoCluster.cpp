#include                  <ImplicitTopoCluster.h>


using namespace std;
using namespace ttk;

ImplicitTopoCluster::ImplicitTopoCluster(){
  clear();
  caches_.resize(1);
  cacheMaps_.resize(1);
  #ifdef TTK_ENABLE_OPENMP
  caches_.resize(threadNumber_);
  cacheMaps_.resize(threadNumber_);
  #endif
}

ImplicitTopoCluster::~ImplicitTopoCluster(){
}

int ImplicitTopoCluster::clear(){
  vertexNumber_ = 0;
  cellNumber_ = 0;
  nodeNumber_ = 0;
  doublePrecision_ = false;

  {
    stringstream msg;
    msg << "[ImplicitTopoCluster] Triangulation cleared." << endl;
    dMsg(cout, msg.str(), detailedInfoMsg);
  }

  return AbstractTriangulation::clear();
}