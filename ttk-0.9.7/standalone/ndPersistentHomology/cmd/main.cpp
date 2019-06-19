/// \author Your Name Here <Your Email Address Here>.
/// \date The Date Here.
///
/// \brief dummy program example.

// include the local headers
#include                  <ttkndPersistentHomology.h>
#include                  <ttkProgramBase.h>

using namespace std;
using namespace ttk;

int main(int argc, char **argv) {

  vtkProgram<ttkndPersistentHomology> program;

  
  int ret = 0;
  ret = program.init(argc, argv);
 
  if(ret != 0)
    return ret;

  program.ttkObject_->SetSelectedFields("func1");
  program.ttkObject_->SetSelectedFields("func2");
  
  // execute data processing
  ret = program.run();
  
  if(ret != 0)
    return ret;
 
  // save the output
  // optional TODO-4:
  // if you want a different kind of output, re-implement the function save().
  ret = program.save();
  /// end of optional TODO-4
  
  return ret;
}
