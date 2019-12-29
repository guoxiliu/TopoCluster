/// \author Your Name Here <Your Email Address Here>.
/// \date The Date Here.
///
/// \brief dummy program example.

// include the local headers
#include <ttkPreprocessStellar.h>
#include <ttkProgramBase.h>

using namespace std;
using namespace ttk;

int main(int argc, char **argv) {

  vtkProgram<ttkPreprocessStellar> program;
  
  int capacity = 1000;

  program.parser_.setArgument("c", &capacity, "Bucket capacity", true);
  
  int ret = 0;
  ret = program.init(argc, argv);
 
  if(ret != 0)
    return ret;

  program.ttkObject_->SetThreshold(capacity);
  
  ret = program.run();
  
  if(ret != 0)
    return ret;
 
  ret = program.save();
  
  return ret;
}
