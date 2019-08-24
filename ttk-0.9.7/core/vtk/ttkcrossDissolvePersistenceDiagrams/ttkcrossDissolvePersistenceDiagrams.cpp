#include                  <ttkcrossDissolvePersistenceDiagrams.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkcrossDissolvePersistenceDiagrams)

int ttkcrossDissolvePersistenceDiagrams::doIt(vector<vtkDataSet *> &inputs, vector<vtkDataSet *> &outputs){

  if(crossDissolvePersistenceDiagrams_ == NULL)
    crossDissolvePersistenceDiagrams_ = new ttk::crossDissolvePersistenceDiagrams();

  Memory m;
  
  vtkDataSet *input = inputs[0];
  vtkDataSet *output = outputs[0];
  
  Triangulation *triangulation = ttkTriangulation::getTriangulation(input);

  if(!triangulation)
    return -1;
  
  triangulation->setWrapper(this);
  crossDissolvePersistenceDiagrams_->setupTriangulation(triangulation);
  crossDissolvePersistenceDiagrams_->setWrapper(this);
 
  vector<void*>* inputScalarFields = new vector<void*>(SelectedFields.size());
  
  cout << SelectedFields.size() << endl;
  for(int i=0; i<SelectedFields.size(); i++)
    (*inputScalarFields)[i] = (input->GetPointData()->GetArray(SelectedFields[i].data())->GetVoidPointer(0));
    
  if(inputScalarFields == 0)
      return -2;

  output->ShallowCopy(input);

  vtkDataArray* outputScalarField_ = NULL;
  if(!outputScalarField_){
    switch(input->GetPointData()->GetArray(SelectedFields[0].data())->GetDataType()){
      
      case VTK_CHAR:
        outputScalarField_ = vtkCharArray::New();
        break;
        
      case VTK_DOUBLE:
        outputScalarField_ = vtkDoubleArray::New();
        break;

      case VTK_FLOAT:
        outputScalarField_ = vtkFloatArray::New();
        break;
       
      case VTK_INT:
        outputScalarField_ = vtkIntArray::New();
        break;

      case VTK_ID_TYPE:
        outputScalarField_ = vtkIdTypeArray::New();
        break;
        
      stringstream msg;
      msg << "[ttkprova] Unsupported data type :(" << endl;
      dMsg(cerr, msg.str(), fatalMsg);
    }
  }

  outputScalarField_->SetNumberOfTuples(input->GetNumberOfPoints());
  outputScalarField_->SetName("Filtration 1D");
  output->GetPointData()->AddArray(outputScalarField_);
  

  // calling the executing package
  crossDissolvePersistenceDiagrams_->setInputDataPointer(inputScalarFields);
  crossDissolvePersistenceDiagrams_->setOutputDataPointer(outputScalarField_->GetVoidPointer(0));
  switch(input->GetPointData()->GetArray(SelectedFields[0].data())->GetDataType()){
    ttkTemplateMacro(crossDissolvePersistenceDiagrams_->execute<VTK_TT>(Alpha));
  }

  outputs[0] = input;
  outHomology(crossDissolvePersistenceDiagrams_->getPH(), outputs[1]);
  outPersistencePairs(crossDissolvePersistenceDiagrams_->getPH(), outputs[2]);
  //print results here
  
  return 0;
}

void ttkcrossDissolvePersistenceDiagrams::outPersistencePairs(PersistentHomology* ph, vtkDataSet* output){
  
  vtkUnstructuredGrid* outputCriticalPoints=vtkUnstructuredGrid::SafeDownCast(output);

  vector<float> criticalPoints = vector<float>();
  vector<char> criticalPointsCellDimension = vector<char>();
  vector<double> filtration = vector<double>();

  ph->readPersistencePairs(criticalPoints,criticalPointsCellDimension,filtration);

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

void ttkcrossDissolvePersistenceDiagrams::outHomology(PersistentHomology* ph, vtkDataSet* output){

  vtkUnstructuredGrid* outputCriticalPoints=vtkUnstructuredGrid::SafeDownCast(output);

  vector<float> criticalPoints = vector<float>();
  vector<char> criticalPointsCellDimension = vector<char>();

  ph->readHomology(criticalPoints,criticalPointsCellDimension);

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