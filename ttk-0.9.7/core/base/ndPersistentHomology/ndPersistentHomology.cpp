#include                  <ndPersistentHomology.h>

using namespace std;
using namespace ttk;

ndPersistentHomology::ndPersistentHomology(){

  inputData_ = NULL;
  triangulation_ = NULL;
  outputData = NULL;
}

ndPersistentHomology::~ndPersistentHomology(){
  delete ph; 
}

