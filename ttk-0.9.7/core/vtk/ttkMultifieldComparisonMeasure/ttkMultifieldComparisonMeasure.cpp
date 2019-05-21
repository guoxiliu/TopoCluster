#include                  <ttkMultifieldComparisonMeasure.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkMultifieldComparisonMeasure)

int ttkMultifieldComparisonMeasure::doIt(vector<vtkDataSet *> &inputs, vector<vtkDataSet *> &outputs){

  Memory m;
  
  vtkDataSet *input = inputs[0];
  vtkDataSet *output = outputs[0];
  
  Triangulation *triangulation = ttkTriangulation::getTriangulation(input);
 
  if(!triangulation)
    return -1;

  output->ShallowCopy(input);
  
  triangulation->setWrapper(this);
  multifieldComparisonMeasure_.setupTriangulation(triangulation);
  multifieldComparisonMeasure_.setWrapper(this);
 
  
  // in the following, the target scalar field of the input is replaced in the 
  // variable 'output' with the result of the computation.
  // if your wrapper produces an output of the same type of the input, you 
  // should proceed in the same way.
  vector<void *>* inputScalarFields = new vector<void *>(SelectedFields.size());

  if(SelectedFields.size() > 0){
    for(int i=0; i<SelectedFields.size(); i++)
        (*inputScalarFields)[i] = (input->GetPointData()->GetArray(SelectedFields[i].data())->GetVoidPointer(0));
  }
  else{
    inputScalarFields->push_back(input->GetPointData()->GetArray(0));
  }
  
  cout << "Are we here? " << endl;

  if(inputScalarFields->size() == 0)
    return -2;
  
  // allocate the memory for the output scalar field
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
      msg << "[ttkMultifieldComparisonMeasure] Unsupported data type :(" << endl;
      dMsg(cerr, msg.str(), fatalMsg);
    }
  }

  cout << "Are we here? 2" << endl;
  // calling the executing package
  multifieldComparisonMeasure_.setInputDataPointer(inputScalarFields);

  cout << "Are we here? 3 " << input->GetNumberOfPoints() << endl;
  
  outputScalarField_->SetNumberOfTuples(input->GetNumberOfPoints());
  outputScalarField_->SetName("CompMeasure");
  multifieldComparisonMeasure_.setOutputDataPointer(outputScalarField_->GetVoidPointer(0));



  switch(input->GetPointData()->GetArray(SelectedFields[0].data())->GetDataType()){
    ttkTemplateMacro(multifieldComparisonMeasure_.execute<VTK_TT>());
  }

  cout << "Are we here? 5" << endl;  

  output->GetPointData()->AddArray(outputScalarField_);

  {
    stringstream msg;
    msg << "[ttkMultifieldComparisonMeasure] Memory usage: " << m.getElapsedUsage() 
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }
  
  return 0;
}
