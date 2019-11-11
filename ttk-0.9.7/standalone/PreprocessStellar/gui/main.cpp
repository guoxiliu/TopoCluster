/// \author Your Name Here <Your Email Address Here>.
/// \date The Date Here.
///
/// \brief dummy GUI program example.

// include the local headers
#include                  <ttkPreprocessStellar.h>
#include                  <ttkUserInterfaceBase.h>

using namespace std;
using namespace ttk;

vtkUserInterface<ttkPreprocessStellar> program;

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
  
  // TODO-1: 
  // specify local parameters to the TTK module with default values.
  int threshold = -1;
  // end of TODO-1

  // TODO-2:
  // register these arguments to the command line parser
  program.parser_.setArgument("I", &threshold,
    "Bucket threshold", true);
  // end of TODO-2
  
  int ret = program.init(argc, argv);
 
  if(ret != 0)
    return ret;

  // TODO-3:
  // change here the arguments of the vtkWrapper that you want to update prior
  // to execution.
  program.ttkObject_->SetThreshold(threshold);
  // end of TODO-3
  
  myKeyHandler myHandler;
  program.setKeyHandler(&myHandler); 


  // TODO-5
  // specify if the display of certain outputs should be disabled
  // vector<int> hiddenOutputs = {0, 2};
  // program.hideOutputs(hiddenOutputs);
  // end of TODO-5

  program.run();
  
  return 0;
}
