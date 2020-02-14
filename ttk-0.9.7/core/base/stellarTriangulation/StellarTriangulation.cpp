#include                  <StellarTriangulation.h>


using namespace std;
using namespace ttk;

StellarTriangulation::StellarTriangulation(){
    clear();
}

StellarTriangulation::~StellarTriangulation(){
  
}

int StellarTriangulation::clear(){
  vertexNumber_ = 0;

  cellNumber_ = 0;
  doublePrecision_ = false;

  {
    stringstream msg;
    msg << "[ExplicitTriangulation] Triangulation cleared." << endl;
    dMsg(cout, msg.str(), detailedInfoMsg);
  }

  return AbstractTriangulation::clear();
}