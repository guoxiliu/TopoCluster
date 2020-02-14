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

  SimplexId vertexNumber = triangulation_->getNumberOfVertices();
  SimplexId edgeNumber = triangulation_->getNumberOfEdges();
  SimplexId triangleNumber = triangulation_->getNumberOfTriangles();
  SimplexId cellNumber = triangulation_->getNumberOfCells();

  std::cout << "[TestStellar] vertex num: " << vertexNumber << ", edge num: "
    << edgeNumber << ", triangle num: " << triangleNumber << ", cell num: " << 
    cellNumber << std::endl;


  SimplexId verticesPerCell = triangulation_->getCellVertexNumber(0);
  SimplexId edgesPerCell = triangulation_->getCellEdgeNumber(0);
  SimplexId trianglesPerCell = triangulation_->getCellTriangleNumber(0);

  // /* test vertex edge relationship */
  // t.reStart();
  // for(SimplexId vertexId = 0; vertexId < vertexNumber; vertexId++){
  //   SimplexId edgeNum = triangulation_->getVertexEdgeNumber(vertexId);
  //   SimplexId edgeId;
  //   for(SimplexId j = 0; j < edgeNum; j++){
  //     int result1 = triangulation_->getVertexEdge(vertexId, j, edgeId);
  //     if(result1){
  //       std::cerr << "[TestStellar] vertexId " << vertexId << " Something wrong in getVertexEdge()! Error code: " << result1 << "\n";
  //     }
  //   }
  // }
  // std::cout << "[TestStellar] Time usage for VE: " << t.getElapsedTime() << " s.\n";

  // /* test vertex triangle relationship */
  // t.reStart();
  // for(SimplexId vertexId = 0; vertexId < vertexNumber; vertexId++){
  //   SimplexId triangleNum = triangulation_->getVertexTriangleNumber(vertexId);
  //   SimplexId triangleId;
  //   for(SimplexId j = 0; j < triangleNum; j++){
  //     int result1 = triangulation_->getVertexTriangle(vertexId, j, triangleId);
  //     if(result1){
  //       std::cerr << "[TestStellar] vertexId " << vertexId << " Something wrong in getVertexTriangle()! Error code: " << result1 << "\n";
  //     }
  //   }
  // }
  // std::cout << "[TestStellar] Time usage for VT: " << t.getElapsedTime() << " s.\n";
  
  // /* test vertex star relationship */
  // t.reStart();
  // for(SimplexId vertexId = 0; vertexId < vertexNumber; vertexId++){
  //   SimplexId cellNum = triangulation_->getVertexStarNumber(vertexId);
  //   SimplexId cellId;
  //   for(SimplexId j = 0; j < cellNum; j++){
  //     int result1 = triangulation_->getVertexStar(vertexId, j, cellId);
  //     if(result1){
  //       std::cout << "[TestStellar] vertexId " << vertexId << " Something wrong in getVertexStar()! Error code: " << result1 << "\n";
  //     }
  //   }
  // }
  // std::cout << "[TestStellar] Time usage for VS: " << t.getElapsedTime() << " s.\n";

  // /* test edge vertex relationship */
  // t.reStart();
  // for(SimplexId edgeId = 0; edgeId < edgeNumber; edgeId++){
  //   SimplexId vertexId;
  //   for(SimplexId j = 0; j < 2; j++){
  //     int result1 = triangulation_->getEdgeVertex(edgeId, j, vertexId);
  //     if(result1){
  //       std::cout << "[TestStellar] edgeId " << edgeId << " Something wrong in getEdgeVertex()! Error code: " << result1 << "\n";
  //     }
  //   }
  // }
  // std::cout << "[TestStellar] Time usage for EV: " << t.getElapsedTime() << " s.\n";


  // /* test edge triangle relationship */
  // t.reStart();
  // for(SimplexId edgeId = 0; edgeId < edgeNumber; edgeId++){
  //   SimplexId triangleNum = triangulation_->getEdgeTriangleNumber(edgeId);
  //   SimplexId triangleId;
  //   for(SimplexId j = 0; j < triangleNum; j++){
  //     int result1 = triangulation_->getEdgeTriangle(edgeId, j, triangleId);
  //     if(result1){
  //       std::cout << "[TestStellar] edgeId " << edgeId << " Something wrong in getEdgeTriangle()! Error code: " << result1 << "\n";
  //     }
  //   }
  // }
  // std::cout << "[TestStellar] Time usage for ET: " << t.getElapsedTime() << " s.\n";
  
  // /* test edge star relationship */
  // t.reStart();
  // for(SimplexId edgeId = 0; edgeId < edgeNumber; edgeId++){
  //   SimplexId cellNum = triangulation_->getEdgeStarNumber(edgeId);
  //   SimplexId cellId;
  //   for(SimplexId j = 0; j < cellNum; j++){
  //     int result1 = triangulation_->getEdgeStar(edgeId, j, cellId);
  //     if(result1){
  //       std::cout << "[TestStellar] edgeId " << edgeId << " Something wrong in getEdgeStar()! Error code: " << result1 << "\n";
  //     }
  //   }
  // }
  // std::cout << "[TestStellar] Time usage for ES: " << t.getElapsedTime() << " s.\n";

  // /* test triangle vertex relationship */
  // t.reStart();
  // for(SimplexId triangleId = 0; triangleId < triangleNumber; triangleId++){
  //   SimplexId vertexId;
  //   for(SimplexId j = 0; j < 3; j++){
  //     int result1 = triangulation_->getTriangleVertex(triangleId, j, vertexId);
  //     if(result1){
  //       std::cout << "[TestStellar] triangleId " << triangleId << " Something wrong in getTriangleVertex()! Error code: " << result1 << "\n";
  //     }
  //   }
  // }
  // std::cout << "[TestStellar] Time usage for TV: " << t.getElapsedTime() << " s.\n";

  // /* test triangle edge relationship */
  // t.reStart();
  // for(SimplexId triangleId = 0; triangleId < triangleNumber; triangleId++){
  //   SimplexId edgeId;
  //   for(SimplexId j = 0; j < 3; j++){
  //     int result1 = triangulation_->getTriangleEdge(triangleId, j, edgeId);
  //     if(result1){
  //       std::cout << "[TestStellar] triangleId " << triangleId << " Something wrong in getTriangleEdge()! Error code: " << result1 << "\n";
  //     }
  //   }
  // }
  // std::cout << "[TestStellar] Time usage for TE: " << t.getElapsedTime() << " s.\n";

  // /* test triangle star relationship */
  // t.reStart();
  // for(SimplexId triangleId = 0; triangleId < edgeNumber; triangleId++){
  //   SimplexId cellNum = triangulation_->getTriangleStarNumber(triangleId);
  //   SimplexId cellId;
  //   for(SimplexId j = 0; j < cellNum; j++){
  //     int result1 = triangulation_->getTriangleStar(triangleId, j, cellId);
  //     if(result1){
  //       std::cout << "[TestStellar] triangleId " << triangleId << " Something wrong in getTriangleStar()! Error code: " << result1 << "\n";
  //     }
  //   }
  // }
  // std::cout << "[TestStellar] Time usage for TS: " << t.getElapsedTime() << " s.\n";

  // /* test cell vertex relationship */
  // t.reStart();
  // for(SimplexId cellId = 0; cellId < cellNumber; cellId++){
  //   SimplexId vertexId;
  //   for(SimplexId j = 0; j < verticesPerCell; j++){
  //     int result1 = triangulation_->getCellVertex(cellId, j, vertexId);
  //     if(result1){
  //       std::cout << "[TestStellar] cellId " << cellId << " Something wrong in getCellVertex()! Error code: " << result1 << "\n";
  //     }
  //   }
  // }
  // std::cout << "[TestStellar] Time usage for CV: " << t.getElapsedTime() << " s.\n";

  // /* test cell edge relationship */
  // t.reStart();
  // for(SimplexId cellId = 0; cellId < cellNumber; cellId++){
  //   SimplexId edgeId;
  //   for(SimplexId j = 0; j < edgesPerCell; j++){
  //     int result1 = triangulation_->getCellEdge(cellId, j, edgeId);
  //     if(result1){
  //       std::cout << "[TestStellar] cellId " << cellId << " Something wrong in getCellEdge()! Error code: " << result1 << "\n";
  //     }
  //   }
  // }
  // std::cout << "[TestStellar] Time usage for CE: " << t.getElapsedTime() << " s.\n";

  // /* test cell triangle relationship */
  // t.reStart();
  // for(SimplexId cellId = 0; cellId < cellNumber; cellId++){
  //   SimplexId triangleId;
  //   for(SimplexId j = 0; j < trianglesPerCell; j++){
  //     int result1 = triangulation_->getCellTriangle(cellId, j, triangleId);
  //     if(result1){
  //       std::cout << "[TestStellar] cellId " << cellId << " Something wrong in getCellTriangle()! Error code: " << result1 << "\n";
  //     }
  //   }
  // }
  // std::cout << "[TestStellar] Time usage for CT: " << t.getElapsedTime() << " s.\n";

  // /* To test the vertex edge */
  // for(SimplexId i = 0; i < edgeNumber; i++){
  //   for(SimplexId j = 0; j < 3; j++){
  //     SimplexId vertexId;
  //     triangulation_->getEdgeVertex(i, j, vertexId);
  //     SimplexId vertexEdgeNum = triangulation_->getVertexEdgeNumber(vertexId);
  //     vector<SimplexId> vertexEdges;
  //     bool hasFound = false;
  //     for(SimplexId k = 0; k < vertexEdgeNum; k++){
  //       SimplexId edgeId;
  //       triangulation_->getVertexEdge(vertexId, k, edgeId);
  //       vertexEdges.push_back(edgeId);
  //       if(edgeId == i){
  //         hasFound = true;
  //         break;
  //       }
  //     }
  //     if(!hasFound){
  //       std::cout << "Edge id " << i << " not found in vertex id " << vertexId << std::endl;
  //       std::cout << "Vertex id " << vertexId << ":[ ";
  //       for(SimplexId tid : vertexEdges){
  //         std::cout << tid << " ";
  //       }
  //       std::cout << "]\n";
  //     }
  //   }
  // }

  // /* To test the vertex triangle */
  // for(SimplexId i = 0; i < triangleNumber; i++){
  //   for(SimplexId j = 0; j < 3; j++){
  //     SimplexId vertexId;
  //     triangulation_->getTriangleVertex(i, j, vertexId);
  //     SimplexId vertexTriangleNum = triangulation_->getVertexTriangleNumber(vertexId);
  //     vector<SimplexId> vertexTriangles;
  //     bool hasFound = false;
  //     for(SimplexId k = 0; k < vertexTriangleNum; k++){
  //       SimplexId triangleId;
  //       triangulation_->getVertexTriangle(vertexId, k, triangleId);
  //       vertexTriangles.push_back(triangleId);
  //       if(triangleId == i){
  //         hasFound = true;
  //         break;
  //       }
  //     }
  //     if(!hasFound){
  //       std::cout << "Triangle id " << i << " not found in vertex id " << vertexId << std::endl;
  //       std::cout << "Vertex id " << vertexId << ":[ ";
  //       for(SimplexId tid : vertexTriangles){
  //         std::cout << tid << " ";
  //       }
  //       std::cout << "]\n";
  //     }
  //   }
  // }

  /* To test the edge triangle */
  for(SimplexId i = 0; i < triangleNumber; i++){
    for(SimplexId j = 0; j < 3; j++){
      SimplexId edgeId;
      triangulation_->getTriangleEdge(i, j, edgeId);
      SimplexId edgeTriangleNum = triangulation_->getEdgeTriangleNumber(edgeId);
      vector<SimplexId> edgeTriangles;
      bool hasFound = false;
      for(SimplexId k = 0; k < edgeTriangleNum; k++){
        SimplexId triangleId;
        triangulation_->getEdgeTriangle(edgeId, k, triangleId);
        edgeTriangles.push_back(triangleId);
        if(triangleId == i){
          hasFound = true;
          break;
        }
      }
      if(!hasFound){
        std::cout << "Triangle id " << i << " not found in edge id " << edgeId << std::endl;
        std::cout << "Edge id " << edgeId << ":[ ";
        for(SimplexId tid : edgeTriangles){
          std::cout << tid << " ";
        }
        std::cout << "]\n";
      }
    }
  }

  /* To test the triangle edge */
  // for(SimplexId i = 0; i < edgeNumber; i++){
  //   SimplexId edgeTriangleNum = triangulation_->getEdgeTriangleNumber(i);
  //   for(SimplexId j = 0; j < edgeTriangleNum; j++){
  //     SimplexId triangleId;
  //     triangulation_->getEdgeTriangle(i, j, triangleId);
  //     vector<SimplexId> triangleEdges;
  //     bool hasFound = false;
  //     for(SimplexId k = 0; k < 3; k++){
  //       SimplexId edgeId;
  //       triangulation_->getTriangleEdge(triangleId, k, edgeId);
  //       triangleEdges.push_back(edgeId);
  //       if(edgeId == i){
  //         hasFound = true;
  //         break;
  //       }
  //     }
  //     if(!hasFound){
  //       std::cout << "Edge id " << i << " not found in triangle id " << triangleId << std::endl;
  //       std::cout << "Triangle id " << triangleId << ":[ ";
  //       for(SimplexId eid : triangleEdges){
  //         std::cout << eid << " ";
  //       }
  //       std::cout << "]\n";
  //     }
  //   }
  // }

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
