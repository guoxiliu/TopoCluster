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
#include                  <chrono>

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

  SimplexId vertexNumber = triangulation_->getNumberOfVertices();
  SimplexId edgeNumber = triangulation_->getNumberOfEdges();
  SimplexId triangleNumber = triangulation_->getNumberOfTriangles();
  SimplexId cellNumber = triangulation_->getNumberOfCells();

  std::cout << "[TestStellar] vertex num: " << vertexNumber << ", edge num: "
    << edgeNumber << ", triangle num: " << triangleNumber << ", cell num: " << 
    cellNumber << std::endl;


  SimplexId verticesPerCell = triangulation_->getCellVertexNumber(0);

  // test vertex edge relationship
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  for(SimplexId vertexId = 0; vertexId < vertexNumber; vertexId++){
  // for each vertex 
    SimplexId edgeNum = triangulation_->getVertexEdgeNumber(vertexId);
    for(SimplexId j = 0; j < edgeNum; j++){
    // for each edge
      SimplexId edgeId;
      if(!triangulation_->getVertexEdge(vertexId, j, edgeId)){
        SimplexId edgeVertexId;
        bool hasFound = false;
        for(SimplexId k = 0; k < 2; k++){
          int result = triangulation_->getEdgeVertex(edgeId, k, edgeVertexId);
          if(!result){
            if(edgeVertexId == vertexId){
              hasFound = true;
              break;
            }
          }
          else{
            std::cout << "[TestStellar] vertexId " << vertexId << ":  Something wrong in getEdgeVertex()! Error code: " << result << "\n";
          }
        }
        if(!hasFound){
          std::cout << "[TestStellar] vertexId " << vertexId << " Cannot find in edge id " << edgeId << "\n";
          triangulation_->getEdgeVertex(edgeId, 0, edgeVertexId);
          std::cout << "edge id " << edgeId <<": " << edgeVertexId << ", ";
          triangulation_->getEdgeVertex(edgeId, 1, edgeVertexId);
          std::cout << edgeVertexId << ".\n";
        }
      }
      else{
        std::cout << "[TestStellar] vertexId " << vertexId << " Something wrong in getVertexEdge()!\n";
      }
    }
  }

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "[TestStellar] Time usage for VE: " << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << " ms\n";

  // test vertex triangle relationship
  begin = std::chrono::steady_clock::now();
  for(SimplexId vertexId = 0; vertexId < vertexNumber; vertexId++){
  // for each vertex 
    SimplexId triangleNum = triangulation_->getVertexTriangleNumber(vertexId);
    for(SimplexId j = 0; j < triangleNum; j++){
    // for each triangle
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
          std::cout << "[TestStellar] vertexId " << vertexId << " Cannot find in triangle id " << triangleId << "\n";
          triangulation_->getTriangleVertex(triangleId, 0, triangleVertexId);
          std::cout << "Triangle id " << triangleId <<": " << triangleVertexId << ", ";
          triangulation_->getTriangleVertex(triangleId, 1, triangleVertexId);
          std::cout << triangleVertexId << ", ";
          triangulation_->getTriangleVertex(triangleId, 2, triangleVertexId);
          std::cout << triangleVertexId << ".\n";
        }
      }
      else{
        std::cout << "[TestStellar] vertexId " << vertexId << " Something wrong in getVertexTriangle()! Error code: " << result1 << "\n";
      }
    }
  }
  end = std::chrono::steady_clock::now();
  std::cout << "[TestStellar] Time usage for VT: " << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << " ms\n";
  

  // test vertex star relationship
  begin = std::chrono::steady_clock::now();
  for(SimplexId vertexId = 0; vertexId < vertexNumber; vertexId++){
  // for each vertex 
    SimplexId cellNum = triangulation_->getVertexStarNumber(vertexId);
    for(SimplexId j = 0; j < cellNum; j++){
    // for each star
      SimplexId cellId;
      int result1 = triangulation_->getVertexStar(vertexId, j, cellId);
      if(!result1){
        SimplexId cellVertexId;
        bool hasFound = false;
        for(SimplexId k = 0; k < verticesPerCell; k++){
          int result2 = triangulation_->getCellVertex(cellId, k, cellVertexId);
          if(!result2){
            if(cellVertexId == vertexId){
              hasFound = true;
              break;
            }
          }
          else{
            std::cout << "[TestStellar] vertexId " << vertexId << ":  Something wrong in getCellVertex()! Error code: " << result2 << "\n";
          }
        }
        if(!hasFound){
          std::cout << "[TestStellar] vertexId " << vertexId << " Cannot find in cell id " << cellId << "\n";
          triangulation_->getCellVertex(cellId, 0, cellVertexId);
          std::cout << "Cell id " << cellId <<": " << cellVertexId << ", ";
          triangulation_->getCellVertex(cellId, 1, cellVertexId);
          std::cout << cellVertexId << ".\n";
        }
      }
      else{
        std::cout << "[TestStellar] vertexId " << vertexId << " Something wrong in getVertexStar()! Error code: " << result1 << "\n";
      }
    }
  }
  end = std::chrono::steady_clock::now();
  std::cout << "[TestStellar] Time usage for VS: " << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << " ms\n";

  // // test edge triangle relationship
  // begin = std::chrono::steady_clock::now();
  // for(SimplexId edgeId = 0; edgeId < edgeNumber; edgeId++){
  // // for each vertex 
  //   SimplexId triangleNum = triangulation_->getEdgeTriangleNumber(edgeId);
  //   for(SimplexId j = 0; j < triangleNum; j++){
  //   // for each triangle
  //     SimplexId triangleId;
  //     int result1 = triangulation_->getEdgeTriangle(edgeId, j, triangleId);
  //     if(!result1){
  //       SimplexId triangleEdgeId;
  //       bool hasFound = false;
  //       for(SimplexId k = 0; k < 3; k++){
  //         int result2 = triangulation_->getTriangleEdge(triangleId, k, triangleEdgeId);
  //         if(!result2){
  //           if(triangleEdgeId == edgeId){
  //             hasFound = true;
  //             break;
  //           }
  //         }
  //         else{
  //           std::cout << "[TestStellar] edgeId " << edgeId << ":  Something wrong in getTriangleEdge()! Error code: " << result2 << "\n";
  //         }
  //       }
  //       if(!hasFound){
  //         std::cout << "[TestStellar] edgeId " << edgeId << " Cannot find in triangle id " << triangleId << "\n";
  //         triangulation_->getTriangleEdge(triangleId, 0, triangleEdgeId);
  //         std::cout << "Triangle id " << triangleId <<": " << triangleEdgeId << ", ";
  //         triangulation_->getTriangleEdge(triangleId, 1, triangleEdgeId);
  //         std::cout << triangleEdgeId << ", ";
  //         triangulation_->getTriangleEdge(triangleId, 2, triangleEdgeId);
  //         std::cout << triangleEdgeId << ".\n";
  //       }
  //     }
  //     else{
  //       std::cout << "[TestStellar] edgeId " << edgeId << " Something wrong in getEdgeTriangle()! Error code: " << result1 << "\n";
  //     }
  //   }
  // }
  // end = std::chrono::steady_clock::now();
  // std::cout << "[TestStellar] Time usage for ET: " << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << " ms\n";

  // // test edge star relationship
  // begin = std::chrono::steady_clock::now();
  // SimplexId edgesPerCell = triangulation_->getCellEdgeNumber(0);
  // for(SimplexId edgeId = 0; edgeId < edgeNumber; edgeId++){
  // // for each vertex 
  //   SimplexId cellNum = triangulation_->getEdgeStarNumber(edgeId);
  //   for(SimplexId j = 0; j < cellNum; j++){
  //   // for each triangle
  //     SimplexId cellId;
  //     int result1 = triangulation_->getEdgeStar(edgeId, j, cellId);
  //     if(!result1){
  //       SimplexId cellEdgeId;
  //       bool hasFound = false;
  //       for(SimplexId k = 0; k < edgesPerCell; k++){
  //         int result2 = triangulation_->getCellEdge(cellId, k, cellEdgeId);
  //         if(!result2){
  //           if(cellEdgeId == edgeId){
  //             hasFound = true;
  //             break;
  //           }
  //         }
  //         else{
  //           std::cout << "[TestStellar] edgeId " << edgeId << ":  Something wrong in getCellEdge()! Error code: " << result2 << "\n";
  //         }
  //       }
  //       if(!hasFound){
  //         std::cout << "[TestStellar] edgeId " << edgeId << " Cannot find in triangle id " << cellId << "\n";
  //       }
  //     }
  //     else{
  //       std::cout << "[TestStellar] edgeId " << edgeId << " Something wrong in getEdgeStar()! Error code: " << result1 << "\n";
  //     }
  //   }
  // }
  // end = std::chrono::steady_clock::now();
  // std::cout << "[TestStellar] Time usage for ES: " << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << " ms\n";

  // // test triangle star relationship
  // begin = std::chrono::steady_clock::now();
  // SimplexId trianglesPerCell = triangulation_->getCellTriangleNumber(0);
  // for(SimplexId triangleId = 0; triangleId < edgeNumber; triangleId++){
  // // for each vertex 
  //   SimplexId cellNum = triangulation_->getTriangleStarNumber(triangleId);
  //   for(SimplexId j = 0; j < cellNum; j++){
  //   // for each triangle
  //     SimplexId cellId;
  //     int result1 = triangulation_->getTriangleStar(triangleId, j, cellId);
  //     if(!result1){
  //       SimplexId cellTriangleId;
  //       bool hasFound = false;
  //       for(SimplexId k = 0; k < trianglesPerCell; k++){
  //         int result2 = triangulation_->getCellTriangle(cellId, k, cellTriangleId);
  //         if(!result2){
  //           if(cellTriangleId == triangleId){
  //             hasFound = true;
  //             break;
  //           }
  //         }
  //         else{
  //           std::cout << "[TestStellar] triangleId " << triangleId << ":  Something wrong in getCellTriangle()! Error code: " << result2 << "\n";
  //         }
  //       }
  //       if(!hasFound){
  //         std::cout << "[TestStellar] triangleId " << triangleId << " Cannot find in triangle id " << cellId << "\n";
  //       }
  //     }
  //     else{
  //       std::cout << "[TestStellar] triangleId " << triangleId << " Something wrong in getTriangleStar()! Error code: " << result1 << "\n";
  //     }
  //   }
  // }
  // end = std::chrono::steady_clock::now();
  // std::cout << "[TestStellar] Time usage for TS: " << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << " ms\n";


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
