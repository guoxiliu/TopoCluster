/// \author Guoxi Liu <guoxil@g.clemson.edu>.
/// \date Jan. 2020.
///

// include the local headers
#include                  <ttkTestStellar.h>
#include                  <ttkProgramBase.h>

using namespace std;
using namespace ttk;

int main(int argc, char **argv) {

  vtkProgram<ttkTestStellar> program;
  
  int ret = 0;
  ret = program.init(argc, argv);
 
  if(ret != 0)
    return ret;

  // execute data processing
  ret = program.run();
  
  if(ret != 0)
    return ret;
 
  ret = program.save();
  
  return ret;
}
