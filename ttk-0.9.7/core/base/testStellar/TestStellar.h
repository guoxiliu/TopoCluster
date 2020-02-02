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
  
  // test vertex edge relationship
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
  //         triangulation_->getEdgeVertex(edgeId, 1, edgeVertexId);
  //         std::cout << edgeVertexId << ".\n";
  //       }
  //     }
  //     else{
  //       std::cout << "[TestStellar] vertexId " << vertexId << " Something wrong in getVertexEdge()!\n";
  //     }
  //   }
  // }

  // // test vertex triangle relationship
  // for(SimplexId vertexId = 0; vertexId < vertexNumber; vertexId++){
  // // for each vertex 
  //   SimplexId triangleNum = triangulation_->getVertexTriangleNumber(vertexId);
  //   for(SimplexId j = 0; j < triangleNum; j++){
  //   // for each triangle
  //     SimplexId triangleId;
  //     int result1 = triangulation_->getVertexTriangle(vertexId, j, triangleId);
  //     if(!result1){
  //       SimplexId triangleVertexId;
  //       bool hasFound = false;
  //       for(SimplexId k = 0; k < 3; k++){
  //         int result2 = triangulation_->getTriangleVertex(triangleId, k, triangleVertexId);
  //         if(!result2){
  //           if(triangleVertexId == vertexId){
  //             hasFound = true;
  //             break;
  //           }
  //         }
  //         else{
  //           std::cout << "[TestStellar] vertexId " << vertexId << ":  Something wrong in getTriangleVertex()! Error code: " << result2 << "\n";
  //         }
  //       }
  //       if(!hasFound){
  //         std::cout << "[TestStellar] vertexId " << vertexId << " Cannot find in triangle id " << triangleId << "\n";
  //         triangulation_->getTriangleVertex(triangleId, 0, triangleVertexId);
  //         std::cout << "Triangle id " << triangleId <<": " << triangleVertexId << ", ";
  //         triangulation_->getTriangleVertex(triangleId, 1, triangleVertexId);
  //         std::cout << triangleVertexId << ", ";
  //         triangulation_->getTriangleVertex(triangleId, 2, triangleVertexId);
  //         std::cout << triangleVertexId << ".\n";
  //       }
  //     }
  //     else{
  //       std::cout << "[TestStellar] vertexId " << vertexId << " Something wrong in getVertexTriangle()! Error code: " << result1 << "\n";
  //     }
  //   }
  // }

  // // test vertex star relationship
  // for(SimplexId vertexId = 0; vertexId < vertexNumber; vertexId++){
  // // for each vertex 
  //   SimplexId cellNum = triangulation_->getVertexStarNumber(vertexId);
  //   SimplexId verticesPerCell = triangulation_->getCellVertexNumber(0);
  //   for(SimplexId j = 0; j < cellNum; j++){
  //   // for each star
  //     SimplexId cellId;
  //     int result1 = triangulation_->getVertexStar(vertexId, j, cellId);
  //     if(!result1){
  //       SimplexId cellVertexId;
  //       bool hasFound = false;
  //       for(SimplexId k = 0; k < verticesPerCell; k++){
  //         int result2 = triangulation_->getCellVertex(cellId, k, cellVertexId);
  //         if(!result2){
  //           if(cellVertexId == vertexId){
  //             hasFound = true;
  //             break;
  //           }
  //         }
  //         else{
  //           std::cout << "[TestStellar] vertexId " << vertexId << ":  Something wrong in getCellVertex()! Error code: " << result2 << "\n";
  //         }
  //       }
  //       if(!hasFound){
  //         std::cout << "[TestStellar] vertexId " << vertexId << " Cannot find in cell id " << cellId << "\n";
  //         triangulation_->getCellVertex(cellId, 0, cellVertexId);
  //         std::cout << "Cell id " << cellId <<": " << cellVertexId << ", ";
  //         triangulation_->getCellVertex(cellId, 1, cellVertexId);
  //         std::cout << cellVertexId << ".\n";
  //       }
  //     }
  //     else{
  //       std::cout << "[TestStellar] vertexId " << vertexId << " Something wrong in getVertexStar()! Error code: " << result1 << "\n";
  //     }
  //   }
  // }

  // // test edge triangle relationship
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
  //           std::cout << "[TestStellar] edgeId " << edgeId << ":  Something wrong in getCellVertex()! Error code: " << result2 << "\n";
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

  // test edge star relationship
  const vector<vector<SimplexId>> *globalEdgeStars = triangulation_->getEdgeStars();
  const vector<vector<SimplexId>> *globalCellEdges = triangulation_->getCellEdges();
  for(SimplexId edgeId = 0; edgeId < edgeNumber; edgeId++){
    for(SimplexId cid : (*globalEdgeStars)[edgeId]){
      bool hasFound = false;
      for(SimplexId eid : (*globalCellEdges)[cid]){
        if(eid == edgeId){
          hasFound = true;
          break;
        }
      }
      if(!hasFound){
        std::cout << "[TestStellar] Edge id " << edgeId << " cannot be found in cell id " << cid << std::endl;
      }
    }
  }

  // test cell edge relationship
  // const vector<vector<SimplexId>> *globalEdgeStars = triangulation_->getEdgeStars();
  // const vector<vector<SimplexId>> *globalCellEdges = triangulation_->getCellEdges();
  for(SimplexId cellId = 0; cellId < cellNumber; cellId++){
    for(SimplexId eid : (*globalCellEdges)[cellId]){
      bool hasFound = false;
      for(SimplexId cid : (*globalEdgeStars)[eid]){
        if(cid == cellId){
          hasFound = true;
          break;
        }
      }
      if(!hasFound){
        std::cout << "[TestStellar] Cell id " << cellId << " cannot be found in edge id " << eid << std::endl;
      }
    }
  }

  // test triangle star relationship
  const vector<vector<SimplexId>> *globalTriangleStars = triangulation_->getTriangleStars();
  const vector<vector<SimplexId>> *globalCellTriangles = triangulation_->getCellTriangles();
  if(globalTriangleStars && globalCellTriangles){
    for(SimplexId triangleId = 0; triangleId < edgeNumber; triangleId++){
      for(SimplexId cid : (*globalTriangleStars)[triangleId]){
        bool hasFound = false;
        for(SimplexId tid : (*globalCellTriangles)[cid]){
          if(tid == triangleId){
            hasFound = true;
            break;
          }
        }
        if(!hasFound){
          std::cout << "[TestStellar] Triangle id " << triangleId << " cannot be found in cell id " << cid << std::endl;
        }
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
