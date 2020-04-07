/// \author Your Name Here <Your Email Address Here>.
/// \date The Date Here.
///
/// \brief dummy program example.

// include the local headers
#include                  <ttkFG_PersistentHomology.h>
#include                  <ttkProgramBase.h>

using namespace std;
using namespace ttk;

int main(int argc, char **argv) {

    vtkProgram<ttkFG_PersistentHomology> program;

    int ret = 0;
    ret = program.init(argc, argv);

    program.ttkObject_->SetScalarField("field");
    program.ttkObject_->SetMinPers(0);
    program.ttkObject_->SetMaxPers(1);
    program.ttkObject_->Setcycles1(false);
    program.ttkObject_->Setcycles2(false);
    program.ttkObject_->SetformanCycles(false);
  //program.ttkObject_->SetScalarField("Elevation");
  ret = program.run();
  
  return 0;
}
