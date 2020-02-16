#include <ttkPreprocessStellar.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkPreprocessStellar)

int ttkPreprocessStellar::doIt(vector<vtkDataSet *> &inputs, vector<vtkDataSet *> &outputs){

  Memory m;

  vtkDataSet *input = inputs[0];
  
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
  for(size_t i = 0; i < vertexArray->size(); i++){
    float x, y, z;
    triangulation->getVertexPoint(vertexArray->at(i), x, y, z);
    points->InsertNextPoint(x, y, z);
    vertexMap[vertexArray->at(i)] = i;
    indices->InsertNextTuple1(nodeArray->at(i));
  }
  outputMesh->SetPoints(points);

  vtkPointData *pointData = outputMesh->GetPointData();
  pointData->AddArray(indices);



  for(auto name : scalarFields){
    vtkDataArray* inputScalars_ = input->GetPointData()->GetArray(name.data());
    vtkDataArray* newField = nullptr;

      switch(inputScalars_->GetDataType()){

        case VTK_CHAR:
          newField = vtkCharArray::New();
          break;

        case VTK_DOUBLE:
          newField = vtkDoubleArray::New();
          break;

        case VTK_FLOAT:
          newField = vtkFloatArray::New();
          break;

        case VTK_INT:
          newField = vtkIntArray::New();
          break;

        case VTK_ID_TYPE:
          newField = vtkIdTypeArray::New();
          break;

      }

    newField->DeepCopy(inputScalars_);
    for(size_t i=0; i<vertexArray->size(); i++){
      newField->SetTuple(i, inputScalars_->GetTuple(vertexArray->at(i)));
    }

    pointData->AddArray(newField);
  }

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
    sort(cell, cell+dimension);
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

  scalarFields.clear();


  {
    stringstream msg;
    msg << "[ttkPreprocessStellar] Memory usage: " << m.getElapsedUsage() 
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }
  
  return 0;
}
