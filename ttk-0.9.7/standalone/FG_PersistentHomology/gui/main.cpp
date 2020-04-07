/// \author Your Name Here <Your Email Address Here>.
/// \date The Date Here.
///
/// \brief dummy GUI program example.

// include the local headers
#include                  <ttkFG_PersistentHomology.h>
#include                  <ttkUserInterfaceBase.h>

using namespace std;
using namespace ttk;

vtkUserInterface<ttkFG_PersistentHomology> program;

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
  
  // // TODO-1: 
  // // specify local parameters to the TTK module with default values.
  // bool someOption = false;
  // int someIntegerArgument = -1;
  // double someDoubleArgument = -1.0;
  // // end of TODO-1

  // // TODO-2:
  // // register these arguments to the command line parser
  // program.parser_.setArgument("D", &someDoubleArgument,
  //   "Some optional double argument", true);
  // program.parser_.setArgument("I", &someIntegerArgument,
  //   "Some optional integer argument", true);
  // program.parser_.setOption("O", &someOption,
  //   "Some option to enable or disable");
  // // end of TODO-2
  
  // int ret = program.init(argc, argv);
 
  // if(ret != 0)
  //   return ret;
  
  // myKeyHandler myHandler;
  // program.setKeyHandler(&myHandler); 

  // program.run();
  
  return 0;
}
