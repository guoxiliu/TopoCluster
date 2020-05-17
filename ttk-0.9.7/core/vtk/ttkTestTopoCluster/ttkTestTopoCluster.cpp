#include                  <ttkTestTopoCluster.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkTestTopoCluster)

int ttkTestTopoCluster::doIt(vector<vtkDataSet *> &inputs, vector<vtkDataSet *> &outputs){

  MemoryUsage m;
  
  vtkDataSet *input = inputs[0];
  vtkDataSet *output = outputs[0];
  
  Triangulation *triangulation = ttkTriangulation::getTriangulation(input);
 
  if(!triangulation)
    return -1;
  
  triangulation->setWrapper(this);
  testTopoCluster_.setupTriangulation(triangulation, CacheSize);
  testTopoCluster_.setWrapper(this);
 
  // use a pointer-base copy for the input data -- to adapt if your wrapper does
  // not produce an output of the type of the input.
  output->ShallowCopy(input);
  
  // in the following, the target scalar field of the input is replaced in the 
  // variable 'output' with the result of the computation.
  // if your wrapper produces an output of the same type of the input, you 
  // should proceed in the same way.
  vtkDataArray *inputScalarField = NULL;
  
  if(ScalarField.length()){
    inputScalarField = input->GetPointData()->GetArray(ScalarField.data());
  }
  else{
    inputScalarField = input->GetPointData()->GetArray(0);
  }
  
  if(!inputScalarField)
    return -2;
  
  // allocate the memory for the output scalar field
  if(!outputScalarField_){
    switch(inputScalarField->GetDataType()){
      
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
      msg << "[ttkTestTopoCluster] Unsupported data type :(" << endl;
      dMsg(cerr, msg.str(), fatalMsg);
    }
  }
  outputScalarField_->SetNumberOfTuples(input->GetNumberOfPoints());
  outputScalarField_->SetName(inputScalarField->GetName());
  
  
  // on the output, replace the field array by a pointer to its processed
  // version
  if(ScalarField.length()){
    output->GetPointData()->RemoveArray(ScalarField.data());
  }
  else{
    output->GetPointData()->RemoveArray(0);
  }
  output->GetPointData()->AddArray(outputScalarField_);
  
  // calling the executing package
  testTopoCluster_.setInputDataPointer(inputScalarField->GetVoidPointer(0));
  testTopoCluster_.setOutputDataPointer(outputScalarField_->GetVoidPointer(0));
  switch(inputScalarField->GetDataType()){
    ttkTemplateMacro(testTopoCluster_.execute<VTK_TT>());
  }
  
  {
    stringstream msg;
    msg << "[ttkTestTopoCluster] Memory usage: " << m.getValue_in_MB(false) 
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }
  
  return 0;
}
