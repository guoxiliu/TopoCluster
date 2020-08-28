/// \author Guoxi Liu <guoxil@g.clemson.edu>.
/// \date Jan. 2020.
///

// include the local headers
#include                  <ttkTestTopoCluster.h>
#include                  <ttkProgramBase.h>

using namespace std;
using namespace ttk;

int main(int argc, char **argv) {

  vtkProgram<ttkTestTopoCluster> program;

  double ratio = 0.1;
  program.parser_.setArgument("r", &ratio, "Cache ratio", true);
  
  int ret = 0;
  ret = program.init(argc, argv);
 
  if(ret != 0)
    return ret;

  // execute data processing

  program.ttkObject_->SetCacheRatio(ratio);
  ret = program.run();
  
  if(ret != 0)
    return ret;
 
  ret = program.save();
  
  return ret;
}
