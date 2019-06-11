/// \author Your Name Here <Your Email Address Here>.
/// \date The Date Here.
///
/// \brief dummy GUI program example.

// include the local headers
#include                  <ttkPersistentHomology.h>
#include                  <ttkUserInterfaceBase.h>

using namespace std;
using namespace ttk;

vtkUserInterface<ttkPersistentHomology> program;

class myKeyHandler : public ttkKeyHandler{
  
  public:
  
    int OnKeyPress(vtkRenderWindowInteractor *interactor, string &key){
     
      stringstream msg;
      msg << "[myKeyHandler] The user pressed the key `" << key << "'." << endl;
      dMsg(cout, msg.str(), infoMsg);
      
      // TODO-4
      // depending on the value of "key", trigger the right functions on the
      // program object (or its contained ttkObject_).
      // end of TODO-4
      
      return 0;
    }
};

int main(int argc, char **argv) {
  

  int ret = 0;
  ret = program.init(argc, argv);
 

  program.ttkObject_->SetScalarField("func1");
  ret = program.run();
  
  myKeyHandler myHandler;
  program.setKeyHandler(&myHandler); 



  program.run();
  
  return 0;
}
