/// \author Your Name Here <Your Email Address Here>.
/// \date The Date Here.
///
/// \brief dummy program example.

// include the local headers
#include                  <ttkPersistentHomology.h>
#include                  <ttkProgramBase.h>

using namespace std;
using namespace ttk;

int main(int argc, char **argv) {

  vtkProgram<ttkPersistentHomology> program;
  
  int ret = 0;
  ret = program.init(argc, argv);
 
  program.ttkObject_->SetScalarField("func1");
  
  //program.ttkObject_->SetScalarField("Elevation");
  ret = program.run();
  
  // if(ret != 0)
  //   return ret;
 
  // // save the output
  // // optional TODO-4:
  // // if you want a different kind of output, re-implement the function save().
  // ret = program.save();
  // /// end of optional TODO-4
  
  return ret;
}
