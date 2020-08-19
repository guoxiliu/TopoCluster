#include                  <ttkFG_Segmentation.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkFG_Segmentation)

int ttkFG_Segmentation::doIt(vector<vtkDataSet *> &inputs, vector<vtkDataSet *> &outputs){

  Memory m;
  
  vtkDataSet *input = inputs[0];
  
  Triangulation *triangulation = ttkTriangulation::getTriangulation(input);
 
  if(!triangulation)
    return -1;
  
  triangulation->setWrapper(this);
  fG_Segmentation_->setupTriangulation(triangulation, CacheSize);
  fG_Segmentation_->setWrapper(this);
  

  vtkDataArray *inputScalarField = NULL;
  
  if(ScalarField.length()){
    inputScalarField = input->GetPointData()->GetArray(ScalarField.data());
  }
  else{
    inputScalarField = input->GetPointData()->GetArray(0);
  }

  int vertices = triangulation->getNumberOfVertices();
  
  if(!inputScalarField)
    return -2;

    fG_Segmentation_->setInputDataPointer(inputScalarField->GetVoidPointer(0));

    cout << "Start computing gradient " << endl;
    Timer t;
    t.reStart();
    switch(inputScalarField->GetDataType()){
        ttkTemplateMacro(fG_Segmentation_->computeIndexing<VTK_TT>());
    }
    cout << "Done computing gradient " << t.getElapsedTime() << endl;

    cout << "Start computing MS " << endl;
    t.reStart();
    fG_Segmentation_->computeMorseSmale();
    cout << "Done computing MS " << t.getElapsedTime() << endl;


//  outputCriticalCells(outputs, triangulation);
//  cout << "Morse complex extracted " << endl;

  return 0;
}

void ttkFG_Segmentation::outputCriticalCells(vector<vtkDataSet *> outputs, Triangulation* triangulation){

  vtkUnstructuredGrid* outputCriticalPoints=vtkUnstructuredGrid::SafeDownCast(outputs[0]);

  vector<float> criticalPoints;
  vector<char> criticalPointsCellDimension;
  vector<vector<int> > indexes;

  fG_Segmentation_->readCriticalPoints(criticalPoints, criticalPointsCellDimension, indexes);
  
  vtkSmartPointer<vtkPoints> points=vtkSmartPointer<vtkPoints>::New();
  //prepare array of critical points dimension
  vtkSmartPointer<vtkCharArray> cellDimensions=vtkSmartPointer<vtkCharArray>::New();
  cellDimensions->SetNumberOfComponents(1);
  cellDimensions->SetName("CellDimension");

  vtkSmartPointer<vtkCharArray> index_labels=vtkSmartPointer<vtkCharArray>::New();
  index_labels->SetNumberOfComponents(1);
  index_labels->SetName("CellId");

  int level=0;
  int count=0;
  for(int i=0; i<criticalPointsCellDimension.size(); i++){
    //add new point to list
    points->InsertNextPoint(criticalPoints[3*i],criticalPoints[3*i+1],criticalPoints[3*i+2]);

    //add new point dimension (corresponding to the original cell originating it)
    cellDimensions->InsertNextTuple1(criticalPointsCellDimension[i]);
    index_labels->InsertNextTuple1(indexes[level][count++]);
    if(count > indexes[level].size()){
      level++;
      count=0;
    }
  }
  outputCriticalPoints->SetPoints(points);

  vtkPointData* pointData=outputCriticalPoints->GetPointData();
  pointData->AddArray(cellDimensions);
  pointData->AddArray(index_labels);
  //--- HERE DONE PREPARING THE CRITICAL POINTS



  vtkUnstructuredGrid* outputDescendingCells=vtkUnstructuredGrid::SafeDownCast(outputs[1]);

  list<Simplex> simplices;
  set<SimplexId> vertices;
  list<int> simplicesPerCell;

  int last_size=0;
  for(int d=1; d<indexes.size(); d++){
    for(int i=0; i<indexes[d].size(); i++){
      fG_Segmentation_->extractDescendingCell(Simplex(d,indexes[d][i]), simplices, vertices);

      simplicesPerCell.push_back(simplices.size()-last_size);
      last_size = simplices.size();
    }
  }

  vtkSmartPointer<vtkPoints> pointsDesc=vtkSmartPointer<vtkPoints>::New();
  vector<SimplexId> updateIndex = vector<SimplexId>(triangulation->getNumberOfVertices(),-1);
  
  count=0;
  for(auto v : vertices){
    updateIndex[v] = count++;
    float x,y,z;
    triangulation->getVertexPoint(v,x,y,z);
    pointsDesc->InsertNextPoint(x,y,z);
  }
  outputDescendingCells->SetPoints(pointsDesc);
  //setting the points for the descending cells

  outputDescendingCells->Allocate(simplices.size());
  for(auto simplex : simplices){
    vector<SimplexId> vSimplex;
    fG_Segmentation_->simplexToVertices(simplex,vSimplex);
    vtkIdType simpl[vSimplex.size()];

    for(int i=0; i<vSimplex.size(); i++){
      simpl[i] = updateIndex[vSimplex[i]];
  }

    if(vSimplex.size() == 2)
      outputDescendingCells->InsertNextCell(VTK_LINE, 2, simpl);
    else if(vSimplex.size() == 3)
      outputDescendingCells->InsertNextCell(VTK_TRIANGLE, 3, simpl);
    else if(vSimplex.size() == 4)
      outputDescendingCells->InsertNextCell(VTK_TETRA, 4, simpl);
  }

  vtkSmartPointer<vtkCharArray> index_labels_desc=vtkSmartPointer<vtkCharArray>::New();
  index_labels_desc->SetNumberOfComponents(1);
  index_labels_desc->SetName("CellIdDesc");

  count=0;
  int sum=0;
  for(auto val : simplicesPerCell){
    for(int i=0; i<val; i++){
      index_labels_desc->InsertNextTuple1(count);
      sum++;
    }
    count++;
  }

  cout << "Compare " << count << " " << sum <<" " << simplices.size() << endl;
  vtkCellData* cellData=outputDescendingCells->GetCellData();
  cellData->AddArray(index_labels_desc);
  //setting the geometry for the descending cells



  vtkUnstructuredGrid* outputAscendingCells=vtkUnstructuredGrid::SafeDownCast(outputs[2]);

  simplices.clear();
  vertices.clear();
  simplicesPerCell.clear();

  last_size=0;
  for(int d=0; d<indexes.size()-1; d++){
    for(int i=0; i<indexes[d].size(); i++){
      fG_Segmentation_->extractAscendingCell(Simplex(d,indexes[d][i]), simplices, vertices);

      simplicesPerCell.push_back(simplices.size()-last_size);
      last_size = simplices.size();
    }
  }

  vtkSmartPointer<vtkPoints> pointsAsc=vtkSmartPointer<vtkPoints>::New();
  updateIndex = vector<SimplexId>(triangulation->getNumberOfVertices(),-1);
  
  count=0;
  for(auto v : vertices){
    updateIndex[v] = count++;
    float x,y,z;
    triangulation->getVertexPoint(v,x,y,z);
    pointsAsc->InsertNextPoint(x,y,z);
  }
  outputAscendingCells->SetPoints(pointsAsc);
  //setting the points for the descending cells

  outputAscendingCells->Allocate(simplices.size());
  for(auto simplex : simplices){
    vector<SimplexId> vSimplex;
    fG_Segmentation_->simplexToVertices(simplex,vSimplex);
    vtkIdType simpl[vSimplex.size()];

    for(int i=0; i<vSimplex.size(); i++){
      simpl[i] = updateIndex[vSimplex[i]];
  }

    if(vSimplex.size() == 2)
      outputAscendingCells->InsertNextCell(VTK_LINE, 2, simpl);
    else if(vSimplex.size() == 3)
      outputAscendingCells->InsertNextCell(VTK_TRIANGLE, 3, simpl);
    else if(vSimplex.size() == 4)
      outputAscendingCells->InsertNextCell(VTK_TETRA, 4, simpl);
  }

  vtkSmartPointer<vtkCharArray> index_labels_asc=vtkSmartPointer<vtkCharArray>::New();
  index_labels_asc->SetNumberOfComponents(1);
  index_labels_asc->SetName("CellIdDesc");

  count=0;
  for(auto val : simplicesPerCell){
    for(int i=0; i<val; i++)
      index_labels_asc->InsertNextTuple1(count);
    count++;
  }
  vtkCellData* cellDataAsc=outputAscendingCells->GetCellData();
  cellDataAsc->AddArray(index_labels_asc);
  //--- HERE DONE PREPARING THE ASCENDING MORSE COMPLEXES

  // cout << "Done preparing " << endl;
}
