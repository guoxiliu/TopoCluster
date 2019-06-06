#include                  <ttkPersistentHomology.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkPersistentHomology)

int ttkPersistentHomology::doIt(vector<vtkDataSet *> &inputs, vector<vtkDataSet *> &outputs){

  Memory m;
  
  vtkDataSet *input = inputs[0];
  vtkDataSet *output = outputs[0];
  
  Triangulation *triangulation = ttkTriangulation::getTriangulation(input);
 
  if(!triangulation)
    return -1;
  
  triangulation->setWrapper(this);
  persistentHomology_.setupTriangulation(triangulation);
  persistentHomology_.setWrapper(this);
 

  vtkDataArray *inputScalarField = NULL;
  
  if(ScalarField.length()){
    inputScalarField = input->GetPointData()->GetArray(ScalarField.data());
  }
  else{
    inputScalarField = input->GetPointData()->GetArray(0);
  }
  
  persistentHomology_.setInputDataPointer(inputScalarField);
  
  switch(inputScalarField->GetDataType()){
    ttkTemplateMacro(persistentHomology_.execute<VTK_TT>());
  }
  
  return 0;
}
