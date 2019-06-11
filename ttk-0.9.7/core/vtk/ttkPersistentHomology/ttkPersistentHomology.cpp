#include                  <ttkPersistentHomology.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkPersistentHomology)

int ttkPersistentHomology::doIt(vector<vtkDataSet *> &inputs, vector<vtkDataSet *> &outputs){

  // Memory m;
  vtkDataSet *input = vtkDataSet::SafeDownCast(inputs[0]);
  
  vtkDataArray *inputScalarField;
  if(ScalarField.length()){
    inputScalarField=input->GetPointData()->GetArray(ScalarField.data());
  }
  else{
    inputScalarField=input->GetPointData()->GetArray(0);
  }

  Triangulation *triangulation = ttkTriangulation::getTriangulation(input);

  if(!triangulation)
    return -1;
  
  triangulation->setWrapper(this);
  persistentHomology_.setWrapper(this);
  persistentHomology_.setupTriangulation(triangulation);
  persistentHomology_.setInputDataPointer(inputScalarField->GetVoidPointer(0));
  
  switch(inputScalarField->GetDataType()){
    ttkTemplateMacro(persistentHomology_.execute<VTK_TT>(true));
  }

  outPersistencePairs(outputs[0]);
  outHomology(outputs[1]);

  return 0;
}

void ttkPersistentHomology::outPersistencePairs(vtkDataSet* output){
  
  vtkUnstructuredGrid* outputCriticalPoints=vtkUnstructuredGrid::SafeDownCast(output);

  vector<float> criticalPoints = vector<float>();
  vector<char> criticalPointsCellDimension = vector<char>();
  vector<double> filtration = vector<double>();

  persistentHomology_.readPersistencePairs(criticalPoints,criticalPointsCellDimension,filtration);

  vtkSmartPointer<vtkPoints> points=vtkSmartPointer<vtkPoints>::New();
  //prepare array of critical points dimension
  vtkSmartPointer<vtkCharArray> cellDimensions=vtkSmartPointer<vtkCharArray>::New();
  cellDimensions->SetNumberOfComponents(1);
  cellDimensions->SetName("CellDimension");
 
  vtkSmartPointer<vtkFloatArray> filtrPoints=vtkSmartPointer<vtkFloatArray>::New();
  filtrPoints->SetNumberOfComponents(1);
  filtrPoints->SetName("Filtration");

  int index=0;
  for(int i=0; i<criticalPointsCellDimension.size(); i++){
    //add new point to list
    points->InsertNextPoint(criticalPoints[3*i],criticalPoints[3*i+1],criticalPoints[3*i+2]);

    //add new point dimension (corresponding to the original cell originating it)
    cellDimensions->InsertNextTuple1(criticalPointsCellDimension[i]);
    filtrPoints->InsertNextTuple1(filtration[index]);

    if(i%2!=0) index++;

  }
  outputCriticalPoints->SetPoints(points);
  vtkPointData* pointData=outputCriticalPoints->GetPointData();
  pointData->AddArray(cellDimensions);
  pointData->AddArray(filtrPoints);

  // vtkSmartPointer<vtkFloatArray> cellFiltration=vtkSmartPointer<vtkFloatArray>::New();
  // cellFiltration->SetNumberOfComponents(1);
  // cellFiltration->SetName("Filtration");

  outputCriticalPoints->Allocate(criticalPointsCellDimension.size()/2.0);
  for(int i=0; i<criticalPointsCellDimension.size(); i=i+2){
    vtkIdType line[2];
    line[0]=i;
    line[1]=i+1;
    outputCriticalPoints->InsertNextCell(VTK_LINE, 2, line);
    // cellFiltration->InsertNextTuple1(filtration[i/2]);
  }

  // vtkCellData* cellData = outputCriticalPoints->GetCellData();
  // cellData->AddArray(cellFiltration);

  cout << criticalPointsCellDimension.size()/2.0 << " persistence pairs" << endl;
}

void ttkPersistentHomology::outHomology(vtkDataSet* output){

  vtkUnstructuredGrid* outputCriticalPoints=vtkUnstructuredGrid::SafeDownCast(output);

  vector<float> criticalPoints = vector<float>();
  vector<char> criticalPointsCellDimension = vector<char>();

  persistentHomology_.readHomology(criticalPoints,criticalPointsCellDimension);

  vtkSmartPointer<vtkPoints> points=vtkSmartPointer<vtkPoints>::New();
  //prepare array of critical points dimension
  vtkSmartPointer<vtkCharArray> cellDimensions=vtkSmartPointer<vtkCharArray>::New();
  cellDimensions->SetNumberOfComponents(1);
  cellDimensions->SetName("CellDimension");
 
  for(int i=0; i<criticalPointsCellDimension.size(); i++){
    //add new point to list
    points->InsertNextPoint(criticalPoints[3*i],criticalPoints[3*i+1],criticalPoints[3*i+2]);

    //add new point dimension (corresponding to the original cell originating it)
    cellDimensions->InsertNextTuple1(criticalPointsCellDimension[i]);

  }
  outputCriticalPoints->SetPoints(points);

  vtkPointData* pointData=outputCriticalPoints->GetPointData();
  pointData->AddArray(cellDimensions);

  cout << criticalPointsCellDimension.size() << " homology classes" << endl;
}