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
  
  if(preprocessStellar_.execute(Threshold)){
    return -1;
  }

  vector<SimplexId> vertexMap(vertexArray->size());
  vtkUnstructuredGrid* outputMesh = vtkUnstructuredGrid::SafeDownCast(outputs[0]);

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkIntArray> indices = vtkSmartPointer<vtkIntArray>::New();

  indices->SetNumberOfComponents(1);
  indices->SetName("_index");

  // insert the vertices in the output mesh
  for(unsigned int i = 0; i < vertexArray->size(); i++){
    float x, y, z;
    triangulation->getVertexPoint(vertexArray->at(i), x, y, z);
    points->InsertNextPoint(x, y, z);
    vertexMap[vertexArray->at(i)] = i;
    indices->InsertNextTuple1(nodeArray->at(i));
  }
  outputMesh->SetPoints(points);

  vtkPointData *pointData = outputMesh->GetPointData();
  pointData->AddArray(indices);

  // insert the cells in the output mesh
  outputMesh->Allocate(cellArray->size());
  int dimension = triangulation->getCellVertexNumber(0);

  for(unsigned int i = 0; i < cellArray->size(); i++){
    vtkIdType cell[dimension];
    for(int j = 0; j < dimension; j++){
        SimplexId  vertexId;
        triangulation->getCellVertex(cellArray->at(i), j, vertexId);
        cell[j] = vertexMap[vertexId];
    }
    if(dimension == 2){
        outputMesh->InsertNextCell(VTK_LINE, 2, cell);
    }else if(dimension == 3){
        outputMesh->InsertNextCell(VTK_TRIANGLE, 3, cell);
    }else if(dimension == 4){
        outputMesh->InsertNextCell(VTK_TETRA, 4, cell);
    }else{
        cerr << "[ttkPreprocessStellar] Should not get here!\n";
    }
  }


  {
    stringstream msg;
    msg << "[ttkPreprocessStellar] Memory usage: " << m.getElapsedUsage() 
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }
  
  return 0;
}
