#include                  <ttkMultivariateGradient.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkMultivariateGradient)

int ttkMultivariateGradient::doIt(vector<vtkDataSet *> &inputs, vector<vtkDataSet *> &outputs){

  Memory m;
  Timer t;

  //set input
  vtkDataSet *input = inputs[0];
  Triangulation *triangulation;

  if(multivariateGradient_ == NULL){
    triangulation = ttkTriangulation::getTriangulation(input);

    if(!triangulation)
      return -1;


    triangulation->setWrapper(this);
    multivariateGradient_ = new MultivariateGradient();
    multivariateGradient_->setupTriangulation(triangulation);
    multivariateGradient_->setWrapper(this);
    multivariateGradient_->setThreadNumber(UseAllCores,ThreadNumber);
    cout << "Triangulation initialized " << SelectedFields.size() << endl;
    vector<void *>* inputScalarFields = new vector<void *>(SelectedFields.size());

    if(SelectedFields.size() > 0){
      for(int i=0; i<SelectedFields.size(); i++)
          (*inputScalarFields)[i] = (input->GetPointData()->GetArray(SelectedFields[i].data())->GetVoidPointer(0));
    }
    else{
      inputScalarFields->push_back(input->GetPointData()->GetArray(0));
    }

    cout << "Scalar field reordererd " << m.getElapsedUsage() << " " << t.getElapsedTime() << endl;

    if(inputScalarFields->size() == 0)
      return -2;

    // calling the executing package
    multivariateGradient_->setInputDataPointer(inputScalarFields);
    switch(input->GetPointData()->GetArray(SelectedFields[0].data())->GetDataType()){

      ttkTemplateMacro(multivariateGradient_->computeGradient<VTK_TT>());
    }

    cout << "Gradient computed " << m.getElapsedUsage() << " " << t.getElapsedTime() << endl;
    //multivariateGradient_->computeCriticalClusters();

  }

  multivariateGradient_->computeCriticalClustersAlt();
  cout << "Critical clusters computed " << m.getElapsedUsage()
  << " " << t.getElapsedTime() << endl;


  multivariateGradient_->simplifyGraph(Path);
  float val = multivariateGradient_->simplifyClusters(CL_size);
  multivariateGradient_->sortClustersBySize();
  cout << "Graph simplified " << Path << " " << val << endl;

  if(!extractInfluence)
  {
    outputCriticalCells(outputs);
    cout << "Done with critical cells" << endl;
    if(extractGradient)
      outputGradient(outputs);
    cout << "Done with gradient" << endl;
    outputClustersGraph(outputs);
    cout << "Done with clusters" << endl;
  }

  else{
    cout << "Computing here " << endl;
    outputOneCluster(outputs);
    if(extractAscending){
      cout << "Computing Ascending" << endl;
      outputAscendingRegionOfInfluence(outputs);
      }
    if(extractDescending){
      cout << "Computing Descending" << endl;
      outputDescendingRegionOfInfluence(outputs);
    }
  }


  {
    stringstream msg;
    msg << "[ttkMultivariateGradient] Memory usage: " << m.getElapsedUsage()
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}


void ttkMultivariateGradient::outputGradient(vector<vtkDataSet *> &outputs){

  vtkUnstructuredGrid* outputGradientGlyphs=vtkUnstructuredGrid::SafeDownCast(outputs[1]);

  vector<float> pairedPoints;
  multivariateGradient_->readOutputGradientPairs(pairedPoints);


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

void ttkMultivariateGradient::outputCriticalCells(vector<vtkDataSet *> &outputs){

  vtkUnstructuredGrid* outputCriticalPoints=vtkUnstructuredGrid::SafeDownCast(outputs[0]);

  vector<float> criticalPoints;
  vector<char> criticalPointsCellDimension;
  vector<char> criticalPointsClusterType;
  vector<int> criticalPointsCluster;

  multivariateGradient_->readOutputCriticalPoints(-1, criticalPoints, criticalPointsCellDimension, criticalPointsClusterType, criticalPointsCluster);

  vtkSmartPointer<vtkPoints> points=vtkSmartPointer<vtkPoints>::New();
  //prepare array of critical points dimension
  vtkSmartPointer<vtkCharArray> cellDimensions=vtkSmartPointer<vtkCharArray>::New();
  cellDimensions->SetNumberOfComponents(1);
  cellDimensions->SetName("CellDimension");

  vtkSmartPointer<vtkCharArray> cluster_labels=vtkSmartPointer<vtkCharArray>::New();
  cluster_labels->SetNumberOfComponents(1);
  cluster_labels->SetName("ClusterId");

  vtkSmartPointer<vtkCharArray> cluster_type=vtkSmartPointer<vtkCharArray>::New();
  cluster_type->SetNumberOfComponents(1);
  cluster_type->SetName("ClusterType");

  for(int i=0; i<criticalPointsCellDimension.size(); i++){
    //add new point to list
    points->InsertNextPoint(criticalPoints[3*i],criticalPoints[3*i+1],criticalPoints[3*i+2]);

    //add new point dimension (corresponding to the original cell originating it)
    cellDimensions->InsertNextTuple1(criticalPointsCellDimension[i]);
    cluster_labels->InsertNextTuple1(criticalPointsCluster[i]);
    cluster_type->InsertNextTuple1(criticalPointsClusterType[i]);
  }
  outputCriticalPoints->SetPoints(points);

  vtkPointData* pointData=outputCriticalPoints->GetPointData();
  pointData->AddArray(cellDimensions);
  pointData->AddArray(cluster_labels);
  pointData->AddArray(cluster_type);
}

void ttkMultivariateGradient::outputClustersGraph(vector<vtkDataSet *> &outputs){

  vtkUnstructuredGrid* outputCriticalClusters=vtkUnstructuredGrid::SafeDownCast(outputs[2]);

  list<float> points;
  list<pair<int,int> > cells;
  vector<int> criticalPointsSize;
  vector<char> criticalPointsClusterType;
  vector<int> criticalPointsCluster;


  //here I should read the clusters and transform everything in a graph
  multivariateGradient_->getOutputClusterGraph(points, cells, criticalPointsSize, criticalPointsClusterType, criticalPointsCluster);

  vtkSmartPointer<vtkPoints> point_list=vtkSmartPointer<vtkPoints>::New();

  vtkSmartPointer<vtkIntArray> cellDimensions=vtkSmartPointer<vtkIntArray>::New();
  cellDimensions->SetNumberOfComponents(1);
  cellDimensions->SetName("CellSize");

  vtkSmartPointer<vtkCharArray> cluster_labels=vtkSmartPointer<vtkCharArray>::New();
  cluster_labels->SetNumberOfComponents(1);
  cluster_labels->SetName("ClusterId");

  vtkSmartPointer<vtkCharArray> cluster_type=vtkSmartPointer<vtkCharArray>::New();
  cluster_type->SetNumberOfComponents(1);
  cluster_type->SetName("ClusterType");

  list<float>::iterator it = points.begin();
  int i=0;
  while(it != points.end()){
    float x = *it; it++;
    float y = *it; it++;
    float z = *it; it++;

    point_list->InsertNextPoint(x,y,z);

    cellDimensions->InsertNextTuple1(criticalPointsSize[i]);
    cluster_labels->InsertNextTuple1(criticalPointsCluster[i]);
    cluster_type->InsertNextTuple1(criticalPointsClusterType[i++]);
  }

  outputCriticalClusters->SetPoints(point_list);

  vtkPointData* pointData=outputCriticalClusters->GetPointData();
  pointData->AddArray(cellDimensions);
  pointData->AddArray(cluster_labels);
  pointData->AddArray(cluster_type);

  outputCriticalClusters->Allocate(cells.size());

  for(auto cell : cells){
    vtkIdType line[2];
    line[0]=cell.first;
    line[1]=cell.second;
    outputCriticalClusters->InsertNextCell(VTK_LINE, 2, line);
  }

}

void ttkMultivariateGradient::outputOneCluster(vector<vtkDataSet *> &outputs){

  int selectedCluster = clusterInfluence;

  vtkUnstructuredGrid* outputCriticalPoints=vtkUnstructuredGrid::SafeDownCast(outputs[3]);

  vector<float> criticalPoints;
  vector<char> criticalPointsCellDimension;
  vector<char> criticalPointsClusterType;
  vector<int> criticalPointsCluster;

  multivariateGradient_->readOutputCriticalPoints(selectedCluster, criticalPoints, criticalPointsCellDimension, criticalPointsClusterType, criticalPointsCluster);

  vtkSmartPointer<vtkPoints> points=vtkSmartPointer<vtkPoints>::New();
  //prepare array of critical points dimension
  vtkSmartPointer<vtkCharArray> cellDimensions=vtkSmartPointer<vtkCharArray>::New();
  cellDimensions->SetNumberOfComponents(1);
  cellDimensions->SetName("CellDimension");

  vtkSmartPointer<vtkCharArray> cluster_labels=vtkSmartPointer<vtkCharArray>::New();
  cluster_labels->SetNumberOfComponents(1);
  cluster_labels->SetName("ClusterId");

  vtkSmartPointer<vtkCharArray> cluster_type=vtkSmartPointer<vtkCharArray>::New();
  cluster_type->SetNumberOfComponents(1);
  cluster_type->SetName("ClusterType");

  for(int i=0; i<criticalPointsCellDimension.size(); i++){
    //add new point to list
    points->InsertNextPoint(criticalPoints[3*i],criticalPoints[3*i+1],criticalPoints[3*i+2]);

    //add new point dimension (corresponding to the original cell originating it)
    cellDimensions->InsertNextTuple1(criticalPointsCellDimension[i]);
    cluster_labels->InsertNextTuple1(criticalPointsCluster[i]);
    cluster_type->InsertNextTuple1(criticalPointsClusterType[i]);
  }
  outputCriticalPoints->SetPoints(points);

  vtkPointData* pointData=outputCriticalPoints->GetPointData();
  pointData->AddArray(cellDimensions);
  pointData->AddArray(cluster_labels);
  pointData->AddArray(cluster_type);

  cout << "Points after i've added them " << criticalPoints.size() << endl;
}


void ttkMultivariateGradient::outputDescendingRegionOfInfluence(vector<vtkDataSet *> &outputs){

  int selectedCluster = clusterInfluence;

  vtkUnstructuredGrid* simplicialSoup = vtkUnstructuredGrid::SafeDownCast(outputs[4]);

  list<float> points;
  list<vector<int> > cells;
  vector<int> nCells(3,0); //(nEdges, nTriangles, nTetrahedra)

  //extract the minimum number of points, edges, triangles and tetrahedra do define the descending region
  multivariateGradient_->readOutputDescendingCell(selectedCluster, points, cells, nCells);

  vtkSmartPointer<vtkPoints> point_list=vtkSmartPointer<vtkPoints>::New();

  //set points
  list<float>::iterator it = points.begin();
  while(it != points.end()){
    float x = *it; it++;
    float y = *it; it++;
    float z = *it; it++;

    point_list->InsertNextPoint(x,y,z);
  }
  simplicialSoup->SetPoints(point_list);
  simplicialSoup->Allocate(nCells[0]);

  for(auto cell : cells){
    vtkIdType simplex[cell.size()];
    for(int i=0; i<cell.size(); i++){
      simplex[i] = cell[i];
  }

    if(cell.size() == 2)
      simplicialSoup->InsertNextCell(VTK_LINE, 2, simplex);
    else if(cell.size() == 3)
      simplicialSoup->InsertNextCell(VTK_TRIANGLE, 3, simplex);
    else if(cell.size() == 4)
      simplicialSoup->InsertNextCell(VTK_TETRA, 4, simplex);
  }

  cout << "Went through here " << cells.size() << " " << points.size() << endl;
}

void ttkMultivariateGradient::outputAscendingRegionOfInfluence(vector<vtkDataSet *> &outputs){

  int selectedCluster = clusterInfluence;
  vtkUnstructuredGrid* simplicialSoup = vtkUnstructuredGrid::SafeDownCast(outputs[5]);

  list<float> points;
  list<vector<int> > cells;
  vector<int> nCells(3,0); //(nEdges, nTriangles, nTetrahedra)

  //extract the minimum number of points, edges, triangles and tetrahedra do define the descending region
  multivariateGradient_->readOutputAscendingCell(selectedCluster, points, cells, nCells);

  vtkSmartPointer<vtkPoints> point_list=vtkSmartPointer<vtkPoints>::New();

  //set points
  list<float>::iterator it = points.begin();
  while(it != points.end()){
    float x = *it; it++;
    float y = *it; it++;
    float z = *it; it++;

    point_list->InsertNextPoint(x,y,z);
  }
  simplicialSoup->SetPoints(point_list);
  simplicialSoup->Allocate(nCells[0]);

  for(auto cell : cells){
    vtkIdType simplex[cell.size()];
    for(int i=0; i<cell.size(); i++){
      simplex[i] = cell[i];
  }

    if(cell.size() == 2)
      simplicialSoup->InsertNextCell(VTK_LINE, 2, simplex);
    else if(cell.size() == 3)
      simplicialSoup->InsertNextCell(VTK_TRIANGLE, 3, simplex);
    else if(cell.size() == 4)
      simplicialSoup->InsertNextCell(VTK_TETRA, 4, simplex);
  }

  cout << "Went through here " << cells.size() << " " << points.size() << endl;


}
