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

  // create the octree
  Octree preOctree(triangulation_, argument);
  for(SimplexId i = 0; i < vertexNumber; i++){
    preOctree.insertVertex(i);
  }

  for(SimplexId i = 0; i < cellNumber; i++){
    preOctree.insertCell(i);
  }

  if(preOctree.verifyTree(vertexNumber)){
    cerr << "[PreprocessStellar] The construction of the tree failed!\n";
    return -1;
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
