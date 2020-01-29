/// \ingroup base
/// \class ttk::TestStellar 
/// \author Guoxi Liu <guoxil@g.clemson.edu>
/// \date Jan. 2020.
///
/// \brief TTK %testStellar processing package.
///
/// %TestStellar is a TTK processing package that takes a scalar field on the input 
/// and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa ttkTestStellar.cpp %for a usage example.

#pragma once

// base code includes
#include                  <Triangulation.h>
#include                  <Wrapper.h>
#include                  <cstdlib>     // for random generator


namespace ttk{
  
  class TestStellar : public Debug{

    public:
        
      TestStellar();
      
      ~TestStellar();

      template <class dataType>
        int execute() const;
    
      inline int setInputDataPointer(void *data){
        inputData_ = data;
        return 0;
      }

      inline int setOutputDataPointer(void *data){
        outputData_ = data;
        return 0;
      }
     
      inline int setupTriangulation(Triangulation *triangulation){
        triangulation_ = triangulation;
        
        if(triangulation_){
          // build edges and triangles
          triangulation_->preprocessEdges();
          triangulation_->preprocessTriangles();
          // vertex related relationships
          triangulation_->preprocessVertexEdges();
          triangulation_->preprocessVertexStars();
          triangulation_->preprocessVertexNeighbors();
          triangulation_->preprocessVertexTriangles();
          // edge related relationships
          triangulation_->preprocessEdgeStars();
          triangulation_->preprocessEdgeTriangles();
          // triangle related relationships
          triangulation_->preprocessTriangleEdges();
          triangulation_->preprocessTriangleStars();
          // cell related relationships
          triangulation_->preprocessCellEdges();
          triangulation_->preprocessCellTriangles();
        }
        
        return 0;
      }
    
    protected:
    
      void                  *inputData_, *outputData_;
      Triangulation         *triangulation_;
  };
}

// template functions
template <class dataType> int ttk::TestStellar::execute() const{

  Timer t;
  
  // check the consistency of the variables -- to adapt
#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation_)
    return -1;
  if(!inputData_)
    return -2;
  if(!outputData_)
    return -3;
#endif

  dataType *outputData = (dataType *) outputData_;
  dataType *inputData = (dataType *) inputData_;

  srand(1);   // initialize the seed
  
  SimplexId vertexNumber = triangulation_->getNumberOfVertices();
  SimplexId edgeNumber = triangulation_->getNumberOfEdges();
  SimplexId triangleNumber = triangulation_->getNumberOfTriangles();
  SimplexId cellNumber = triangulation_->getNumberOfCells();

  std::cout << "[TestStellar] vertex num: " << vertexNumber << ", edge num: "
    << edgeNumber << ", triangle num: " << triangleNumber << ", cell num: " << 
    cellNumber << std::endl;
  
  // test vertex edge relationships
  // for(SimplexId vertexId = 0; vertexId < vertexNumber; vertexId++){
  // // for each vertex 
  //   SimplexId edgeNum = triangulation_->getVertexEdgeNumber(vertexId);
  //   for(SimplexId j = 0; j < edgeNum; j++){
  //   // for each edge
  //     SimplexId edgeId;
  //     if(!triangulation_->getVertexEdge(vertexId, j, edgeId)){
  //       SimplexId edgeVertexId;
  //       bool hasFound = false;
  //       for(SimplexId k = 0; k < 2; k++){
  //         int result = triangulation_->getEdgeVertex(edgeId, k, edgeVertexId);
  //         if(!result){
  //           if(edgeVertexId == vertexId){
  //             hasFound = true;
  //             break;
  //           }
  //         }
  //         else{
  //           std::cout << "[TestStellar] vertexId " << vertexId << ":  Something wrong in getEdgeVertex()! Error code: " << result << "\n";
  //         }
  //       }
  //       if(!hasFound){
  //         std::cout << "[TestStellar] vertexId " << vertexId << " Cannot find in edge id " << edgeId << "\n";
  //         triangulation_->getEdgeVertex(edgeId, 0, edgeVertexId);
  //         std::cout << "edge id " << edgeId <<": " << edgeVertexId << ", ";
  //         triangulation_->getEdgeVertex(edgeId, 0, edgeVertexId);
  //         std::cout << edgeVertexId << ".\n";
  //       }
  //     }
  //     else{
  //       std::cout << "[TestStellar] vertexId " << vertexId << " Something wrong in getVertexEdge()!\n";
  //     }
  //   }
  // }

  // test vertex triangle relationships
  for(SimplexId vertexId = 0; vertexId < vertexNumber; vertexId++){
  // for each vertex 
    SimplexId triangleNum = triangulation_->getVertexTriangleNumber(vertexId);
    for(SimplexId j = 0; j < triangleNum; j++){
    // for each Triangle
      SimplexId triangleId;
      int result1 = triangulation_->getVertexTriangle(vertexId, j, triangleId);
      if(!result1){
        SimplexId triangleVertexId;
        bool hasFound = false;
        for(SimplexId k = 0; k < 3; k++){
          int result2 = triangulation_->getTriangleVertex(triangleId, k, triangleVertexId);
          if(!result2){
            if(triangleVertexId == vertexId){
              hasFound = true;
              break;
            }
          }
          else{
            std::cout << "[TestStellar] vertexId " << vertexId << ":  Something wrong in getTriangleVertex()! Error code: " << result2 << "\n";
          }
        }
        if(!hasFound){
          std::cout << "[TestStellar] vertexId " << vertexId << " Cannot find in Triangle id " << triangleId << "\n";
          triangulation_->getTriangleVertex(triangleId, 0, triangleVertexId);
          std::cout << "Triangle id " << triangleId <<": " << triangleVertexId << ", ";
          triangulation_->getTriangleVertex(triangleId, 0, triangleVertexId);
          std::cout << triangleVertexId << ".\n";
        }
      }
      else{
        std::cout << "[TestStellar] vertexId " << vertexId << " Something wrong in getVertexTriangle()! Error code: " << result1 << "\n";
      }
    }
  }
  

  // init the output -- to adapt
  for(SimplexId i = 0; i < vertexNumber; i++){
    outputData[i] = inputData[i];
  }
  
  {
    std::stringstream msg;
    msg << "[TestStellar] Data-set (" << vertexNumber
      << " points) processed in "
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }
  
  return 0;
}
