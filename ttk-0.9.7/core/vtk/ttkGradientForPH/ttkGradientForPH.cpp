#include                  <ttkGradientForPH.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkGradientForPH)

int ttkGradientForPH::doIt(vector<vtkDataSet *> &inputs, vector<vtkDataSet *> &outputs){

  //set input
  vtkDataSet *input = inputs[0];
  Triangulation *triangulation;

  gradientForPH_ = new GradientForPH();

  triangulation = ttkTriangulation::getTriangulation(input);

  if(!triangulation)
    return -1;

  gradientForPH_->setupTriangulation(triangulation);
  gradientForPH_->setWrapper(this);


  void* field = input->GetPointData()->GetArray(SelectedField.data())->GetVoidPointer(0);

    // cout << "Triangulation initialized " << SelectedFields.size() << endl;
    // vector<void *>* inputScalarFields = new vector<void *>(SelectedFields.size());

    // if(SelectedFields.size() > 0){
    //   for(int i=0; i<SelectedFields.size(); i++)
    //       (*inputScalarFields)[i] = (input->GetPointData()->GetArray(SelectedFields[i].data())->GetVoidPointer(0));
    // }
    // else{
    //   inputScalarFields->push_back(input->GetPointData()->GetArray(0));
    // }

    // cout << "Scalar field reordererd " << m.getElapsedUsage() << " " << t.getElapsedTime() << endl;

    // if(inputScalarFields->size() == 0)
    //   return -2;

    // calling the executing package
    gradientForPH_->setInputDataPointer(field);
    switch(input->GetPointData()->GetArray(SelectedField.data())->GetDataType()){
      ttkTemplateMacro(gradientForPH_->computeGradient<VTK_TT>());
    }

  outputCriticalCells(outputs);
  if(extractGradient)
    outputGradient(outputs);
    

  return 0;
}


void ttkGradientForPH::outputGradient(vector<vtkDataSet *> &outputs){

  vtkUnstructuredGrid* outputGradientGlyphs=vtkUnstructuredGrid::SafeDownCast(outputs[1]);

  vector<float> pairedPoints;
  gradientForPH_->readOutputGradientPairs(pairedPoints);


  vtkSmartPointer<vtkPoints> pairOrigins=vtkSmartPointer<vtkPoints>::New();
  for(int i=0; i<pairedPoints.size()/3; i++){
    pairOrigins->InsertNextPoint(pairedPoints[3*i],pairedPoints[3*i+1],pairedPoints[3*i+2]);
  }
  outputGradientGlyphs->SetPoints(pairOrigins);
  outputGradientGlyphs->Allocate(pairedPoints.size()/6);

  int vertex =0;
   for(int i=0; i<pairedPoints.size(); i=i+6){
    vtkIdType line[2];
    line[0]=vertex++;
    line[1]=vertex++;
    outputGradientGlyphs->InsertNextCell(VTK_LINE, 2, line);
  }

}

void ttkGradientForPH::outputCriticalCells(vector<vtkDataSet *> &outputs){

  vtkUnstructuredGrid* outputCriticalPoints=vtkUnstructuredGrid::SafeDownCast(outputs[0]);

  vector<float> criticalPoints;
  vector<char> criticalPointsCellDimension;

  gradientForPH_->readOutputCriticalPoints(-1, criticalPoints, criticalPointsCellDimension);

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
}
