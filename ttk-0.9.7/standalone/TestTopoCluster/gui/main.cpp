/// \author Guoxi Liu <guoxil@g.clemson.edu>.
/// \date Jan. 2020.
///

// include the local headers
#include                  <ttkTestTopoCluster.h>
#include                  <ttkUserInterfaceBase.h>

using namespace std;
using namespace ttk;

vtkUserInterface<ttkTestTopoCluster> program;

class myKeyHandler : public ttkKeyHandler{
  
  public:
  
    int OnKeyPress(vtkRenderWindowInteractor *interactor, string &key){
     
      stringstream msg;
      msg << "[myKeyHandler] The user pressed the key `" << key << "'." << endl;
      dMsg(cout, msg.str(), infoMsg);

      return 0;
    }
};

int main(int argc, char **argv) {
  
  
  int ret = program.init(argc, argv);
 
  if(ret != 0)
    return ret;

  myKeyHandler myHandler;
  program.setKeyHandler(&myHandler); 

  program.run();
  
  return 0;
}
