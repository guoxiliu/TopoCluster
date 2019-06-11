#include                  <ttkmapFunctionTo3D.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkmapFunctionTo3D)

int ttkmapFunctionTo3D::doIt(vector<vtkDataSet *> &inputs, vector<vtkDataSet *> &outputs){

  Memory m;
  
  vtkDataSet *input = vtkDataSet::SafeDownCast(inputs[0]);

  vtkPointSet *output = vtkPointSet::SafeDownCast(outputs[0]);
  
 
  vtkDataArray *inputScalarField = input->GetPointData()->GetArray(ScalarField.data());
  output->ShallowCopy(input);

  vtkPoints * points = output->GetPoints();

  for(int i=0; i<output->GetNumberOfPoints(); i++){
    std::cout << inputScalarField->GetTuple1(i) << std::endl;
    double* coords = points->GetPoint(i);
    coords[2] = inputScalarField->GetTuple1(i)*Mult;
    points->InsertPoint(i,coords);
  }

  return 0;
}
