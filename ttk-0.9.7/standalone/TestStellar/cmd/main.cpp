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
  
  int size = 100;
  program.parser_.setArgument("s", &size, "Cache size", true);

  int ret = 0;
  ret = program.init(argc, argv);
 
  if(ret != 0)
    return ret;

  // execute data processing
  program.ttkObject_->SetCacheSize(size);
  ret = program.run();
  
  if(ret != 0)
    return ret;

  ret = program.save();
  
  return ret;
}
