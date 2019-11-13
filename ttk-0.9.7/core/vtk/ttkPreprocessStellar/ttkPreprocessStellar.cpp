#include <ttkPreprocessStellar.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkPreprocessStellar)

int ttkPreprocessStellar::doIt(vector<vtkDataSet *> &inputs, vector<vtkDataSet *> &outputs){

  Memory m;
  
  // vtkUnstructuredGrid *input = vtkUnstructuredGrid::SafeDownCast(inputs[0]);
  vtkDataSet *input = inputs[0];
  // vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(outputs[0]);;
  
  Triangulation *triangulation = ttkTriangulation::getTriangulation(input);
 
  if(!triangulation)
    return -1;
  
  triangulation->setWrapper(this);
  preprocessStellar_.setupTriangulation(triangulation);
  preprocessStellar_.setWrapper(this);

  // pass necessary parameters for execution
  vector<SimplexId> *vertexArray = new vector<SimplexId>();
  vector<SimplexId> *nodeArray = new vector<SimplexId>();
  vector<SimplexId> *cellArray = new vector<SimplexId>();

  preprocessStellar_.setVerticesPointer(vertexArray);
  preprocessStellar_.setNodesPointer(nodeArray);
  preprocessStellar_.setCellsPointer(cellArray);
  
  preprocessStellar_.execute(Threshold);


  {
    stringstream msg;
    msg << "[ttkPreprocessStellar] Memory usage: " << m.getElapsedUsage() 
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }
  
  return 0;
}
