#include                  <ExplicitTopoCluster.h>


using namespace std;
using namespace ttk;

ExplicitTopoCluster::ExplicitTopoCluster(){
    clear();
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