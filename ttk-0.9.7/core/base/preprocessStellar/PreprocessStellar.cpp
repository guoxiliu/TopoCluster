#include <PreprocessStellar.h>

using namespace std;
using namespace ttk;

PreprocessStellar::PreprocessStellar(){
  vertices = nullptr;
  nodes = nullptr;
  cells = nullptr;
}

PreprocessStellar::~PreprocessStellar(){
  
}



// execute function
int ttk::PreprocessStellar::execute(
  const int &argument) const{

  Timer t;
  
  // check the consistency of the variables -- to adapt
#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation_)
    return -1;
  if(!vertices)
    return -2;
  if(!nodes)
    return -3;
  if(!cells)
    return -4;
#endif

  SimplexId vertexNumber = triangulation_->getNumberOfVertices();
  SimplexId cellNumber = triangulation_->getNumberOfCells();

  // the following open-mp processing is only relevant for embarrassingly 
  // parallel algorithms (such as smoothing) -- to adapt
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) 
#endif

  // create the octree
  Octree preOctree(triangulation_, argument);
  for(SimplexId i = 0; i < vertexNumber; i++){
    preOctree.insertVertex(i);
  }

  preOctree.verifyTree(vertexNumber);

  // SimplexId testCell = 0;
  // preOctree.insertCell(testCell);
  for(SimplexId i = 0; i < cellNumber; i++){
    preOctree.insertCell(i);
  }

  std::vector<SimplexId> *vertexVec = static_cast<std::vector<SimplexId>*>(vertices);
  std::vector<SimplexId> *nodeVec = static_cast<std::vector<SimplexId>*>(nodes);
  std::vector<SimplexId> *cellVec = static_cast<std::vector<SimplexId>*>(cells);
  preOctree.reindex(vertexVec, nodeVec, cellVec);
  std::cout << "[PreprocessStellar] Size of vertexVec: " << vertexVec->size();
  std::cout << "; Size of cellVec: " << cellVec->size() << std::endl;

  {
    std::stringstream msg;
    msg << "[PreprocessStellar] Data-set (" << vertexNumber
      << " points) processed in "
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }
  
  return 0;
}
