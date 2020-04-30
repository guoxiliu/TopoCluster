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
#include                  <MemoryUsage.h>

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
     
      inline int setupTriangulation(Triangulation *triangulation, const int &size){
        triangulation_ = triangulation;
        
        if(triangulation_){
          Timer t;
          t.getStartTime();
          triangulation_->setCacheSize(size);
          // build edges and triangles
          triangulation_->preprocessEdges();
          triangulation_->preprocessTriangles();
          // boundary relationships
          triangulation_->preprocessBoundaryVertices();
          triangulation_->preprocessBoundaryEdges();
          triangulation_->preprocessBoundaryTriangles();
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
          std::cout << "[TestStellar] Time usage for preprocessing: " << t.getElapsedTime() << " s.\n";
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

  Timer t, tot;
  
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


  // std::cout << "dimension: " << triangulation_->getDimensionality() << std::endl;
  // t.reStart();
  // int boundaryVertexNum = 0, boundaryEdgeNum = 0, boundaryTriangleNum = 0;
  // for(SimplexId vid = 0; vid < vertexNumber; vid++){
  //   if(triangulation_->isVertexOnBoundary(vid)){
  //     boundaryVertexNum++;
  //   }
  // }
  // std::cout << "[TestStellar] Boundary vertex number: " << boundaryVertexNum << std::endl;

  // for(SimplexId eid = 0; eid < edgeNumber; eid++){
  //   if(triangulation_->isEdgeOnBoundary(eid)){
  //     boundaryEdgeNum++;
  //   }
  // }
  // std::cout << "[TestStellar] Boundary edge number: " << boundaryEdgeNum << std::endl;

  // for(SimplexId tid = 0; tid < triangleNumber; tid++){
  //   if(triangulation_->isTriangleOnBoundary(tid)){
  //     boundaryTriangleNum++;
  //   }
  // }
  // std::cout << "[TestStellar] Boundary triangle number: " << boundaryTriangleNum << std::endl;
  // std::cout << "[TestStellar] Time usage for getting boundaries: " << t.getElapsedTime() << " s.\n";


  // const std::vector<std::vector<SimplexId>> *vertexNeighbors = triangulation_->getVertexNeighbors();
  // for(SimplexId vid = 0; vid < vertexNumber; vid++){
  //   for(SimplexId neighbor : (*vertexNeighbors)[vid]){
  //     auto it = find((*vertexNeighbors)[neighbor].begin(), (*vertexNeighbors)[neighbor].end(), vid);
  //     if(it == (*vertexNeighbors)[neighbor].end()){
  //       std::cout << "[TestStellar] Cannot find " << vid << " in " << neighbor << "'s neighbors.\n";
  //       break;
  //     }
  //   }
  // }
  

  // // VE relation
  // for(SimplexId vid = 0; vid < vertexNumber; vid++){
  //   int edgeNum = triangulation_->getVertexEdgeNumber(vid);
  //   for(int i = 0; i < edgeNum; i++){
  //     SimplexId edgeId;
  //     int res = triangulation_->getVertexEdge(vid, i, edgeId);
  //     if(res){
  //       std::cout << "[TestStellar] Cannot get vertex edge for vertex id " << vid << ", error code: " << res << "\n";
  //       return -4;
  //     }
  //     bool hasFound = false;
  //     for(int j = 0; j < 2; j++){
  //       SimplexId vertexId;
  //       triangulation_->getEdgeVertex(edgeId, j, vertexId);
  //       if(vertexId == vid){
  //         hasFound = true;
  //         break;
  //       }
  //     }
  //     if(!hasFound){
  //       std::cout << "[TestStellar] Cannot find vertex " << vid << " in edge " << edgeId << "\n";
  //       break;
  //     }
  //   }
  // } 


  // // VF relation
  // for(SimplexId vid = 0; vid < vertexNumber; vid++){
  //   int triangleNum = triangulation_->getVertexTriangleNumber(vid);
  //   for(int i = 0; i < triangleNum; i++){
  //     SimplexId triangleId;
  //     int res = triangulation_->getVertexTriangle(vid, i, triangleId);
  //     if(res){
  //       std::cout << "[TestStellar] Cannot get vertex triangle for vertex id " << vid << ", error code: " << res << "\n";
  //       return -4;
  //     }
  //     bool hasFound = false;
  //     for(int j = 0; j < 3; j++){
  //       SimplexId vertexId;
  //       triangulation_->getTriangleVertex(triangleId, j, vertexId);
  //       if(vertexId == vid){
  //         hasFound = true;
  //         break;
  //       }
  //     }
  //     if(!hasFound){
  //       std::cout << "[TestStellar] Cannot find vertex " << vid << " in triangle " << triangleId << "\n";
  //       break;
  //     }
  //   }
  // }

  // // VT relation
  // for(SimplexId vid = 0; vid < vertexNumber; vid++){
  //   int starNum = triangulation_->getVertexStarNumber(vid);
  //   for(int i = 0; i < starNum; i++){
  //     SimplexId cellId;
  //     int res = triangulation_->getVertexStar(vid, i, cellId);
  //     if(res){
  //       std::cout << "[TestStellar] Cannot get vertex edge for vertex id " << vid << ", error code: " << res << "\n";
  //       return -4;
  //     }
  //     bool hasFound = false;
  //     for(int j = 0; j < 4; j++){
  //       SimplexId vertexId;
  //       triangulation_->getCellVertex(cellId, j, vertexId);
  //       if(vertexId == vid){
  //         hasFound = true;
  //         break;
  //       }
  //     }
  //     if(!hasFound){
  //       std::cout << "[TestStellar] Cannot find vertex " << vid << " in cell " << cellId << "\n";
  //       break;
  //     }
  //   }
  // }

  // // EV relation
  // for(SimplexId eid = 0; eid < edgeNumber; eid++){
  //   for(int i = 0; i < 2; i++){
  //     SimplexId vid;
  //     triangulation_->getEdgeVertex(eid, i, vid);
  //     SimplexId edgeNum = triangulation_->getVertexEdgeNumber(vid);
  //     bool hasFound = false;
  //     for(int j = 0; j < edgeNum; j++){
  //       SimplexId edgeid;
  //       int res = triangulation_->getVertexEdge(vid, j, edgeid);
  //       if(res){
  //         std::cout << "[TestStellar] Cannot get vertex edge for vertex id " << vid << ", error code: " << res << "\n";
  //         return -4;
  //       }
  //       if(edgeid == eid){
  //         hasFound = true;
  //         break;
  //       }
  //     }
  //     if(!hasFound){
  //       std::cout << "[TestStellar] Cannot find edge " << eid << " in vertex " << vid << "\n";
  //       break;
  //     }
  //   }
  // }

  // // EF relation
  // for(SimplexId eid = 0; eid < edgeNumber; eid++){
  //   SimplexId triangleNum = triangulation_->getEdgeTriangleNumber(eid);
  //   for(int i = 0; i < triangleNum; i++){
  //     SimplexId triangleId;
  //     int res = triangulation_->getEdgeTriangle(eid, i, triangleId);
  //     if(res){
  //       std::cout << "[TestStellar] Cannot get edge triangle for edge id " << eid << ", error code: " << res << "\n";
  //       return -4;
  //     }
  //     bool findEdge = false;
  //     SimplexId edgeId;
  //     for(int i = 0; i < 3; i++){
  //       triangulation_->getTriangleEdge(triangleId, i, edgeId);
  //       if(edgeId == eid){
  //         findEdge = true;
  //         break;
  //       }
  //     }
  //     if(!findEdge){
  //       std::cout << "[TestStellar] Cannot find edge " << eid << " in triangle " << triangleId << "\n";
  //       break;
  //     }
  //   }
  // }

  // // ET relation
  // for(SimplexId eid = 0; eid < edgeNumber; eid++){
  //   SimplexId starNum = triangulation_->getEdgeStarNumber(eid);
  //   for(int i = 0; i < starNum; i++){
  //     SimplexId cellId;
  //     int res = triangulation_->getEdgeStar(eid, i, cellId);
  //     if(res){
  //       std::cout << "[TestStellar] Cannot get edge star for edge id " << eid << ", error code: " << res << "\n";
  //       return -4;
  //     }
  //     bool findEdge = false;
  //     SimplexId edgeId;
  //     for(int i = 0; i < 6; i++){
  //       triangulation_->getCellEdge(cellId, i, edgeId);
  //       if(edgeId == eid){
  //         findEdge = true;
  //         break;
  //       }
  //     }
  //     if(!findEdge){
  //       std::cout << "[TestStellar] Cannot find edge " << eid << " in cell " << cellId << "\n";
  //       break;
  //     }
  //   }
  // }

  // // FV relation
  // for(SimplexId tid = 0; tid < triangleNumber; tid++){
  //   for(int i = 0; i < 3; i++){
  //     SimplexId vid;
  //     triangulation_->getTriangleVertex(tid, i, vid);
  //     SimplexId triangleNum = triangulation_->getVertexTriangleNumber(vid);
  //     bool hasFound = false;
  //     for(int j = 0; j < triangleNum; j++){
  //       SimplexId triangleid;
  //       int res = triangulation_->getVertexTriangle(vid, j, triangleid);
  //       if(res){
  //         std::cout << "[TestStellar] Cannot get vertex triangle for vertex id " << vid << ", error code: " << res << "\n";
  //         return -4;
  //       }
  //       if(triangleid == tid){
  //         hasFound = true;
  //         break;
  //       }
  //     }
  //     if(!hasFound){
  //       std::cout << "[TestStellar] Cannot find triangle " << tid << " in vertex " << vid << "\n";
  //       break;
  //     }
  //   }
  // }

  // // FE relation
  // for(SimplexId tid = 0; tid < triangleNumber; tid++){
  //   for(int i = 0; i < 3; i++){
  //     SimplexId eid;
  //     triangulation_->getTriangleEdge(tid, i, eid);
  //     SimplexId triangleNum = triangulation_->getEdgeTriangleNumber(eid);
  //     bool hasFound = false;
  //     for(int j = 0; j < triangleNum; j++){
  //       SimplexId triangleid;
  //       int res = triangulation_->getEdgeTriangle(eid, j, triangleid);
  //       if(res){
  //         std::cout << "[TestStellar] Cannot get edge triangle for edge id " << eid << ", error code: " << res << "\n";
  //         return -4;
  //       }
  //       if(triangleid == tid){
  //         hasFound = true;
  //         break;
  //       }
  //     }
  //     if(!hasFound){
  //       std::cout << "[TestStellar] Cannot find triangle " << tid << " in edge " << eid << "\n";
  //       break;
  //     }
  //   }
  // }

  // // FT relation
  // for(SimplexId tid = 0; tid < triangleNumber; tid++){
  //   SimplexId starNum = triangulation_->getTriangleStarNumber(tid);
  //   for(int i = 0; i < starNum; i++){
  //     SimplexId cellId;
  //     int res = triangulation_->getTriangleStar(tid, i, cellId);
  //     if(res){
  //       std::cout << "[TestStellar] Cannot get traingle star for traingle id " << tid << ", error code: " << res << "\n";
  //       return -4;
  //     }
  //     bool findTriangle = false;
  //     SimplexId triangleid;
  //     for(int i = 0; i < 4; i++){
  //       triangulation_->getCellTriangle(cellId, i, triangleid);
  //       if(triangleid == tid){
  //         findTriangle = true;
  //         break;
  //       }
  //     }
  //     if(!findTriangle){
  //       std::cout << "[TestStellar] Cannot find triangle " << tid << " in cell " << cellId << "\n";
  //       break;
  //     }
  //   }
  // }

  // // TV relation
  // for(SimplexId cid = 0; cid < cellNumber; cid++){
  //   for(int i = 0; i < 4; i++){
  //     SimplexId vid;
  //     triangulation_->getCellVertex(cid, i, vid);
  //     SimplexId starNum = triangulation_->getVertexStarNumber(vid);
  //     bool hasFound = false;
  //     for(int j = 0; j < starNum; j++){
  //       SimplexId starid;
  //       int res = triangulation_->getVertexStar(vid, j, starid);
  //       if(res){
  //         std::cout << "[TestStellar] Cannot get vertex star for vertex id " << vid << ", error code: " << res << "\n";
  //         return -4;
  //       }
  //       if(starid == cid){
  //         hasFound = true;
  //         break;
  //       }
  //     }
  //     if(!hasFound){
  //       std::cout << "[TestStellar] Cannot find cell " << cid << " in vertex " << vid << "\n";
  //       break;
  //     }
  //   }
  // }

  // // TE relation
  // for(SimplexId cid = 0; cid < cellNumber; cid++){
  //   for(int i = 0; i < 6; i++){
  //     SimplexId eid;
  //     triangulation_->getCellEdge(cid, i, eid);
  //     SimplexId starNum = triangulation_->getEdgeStarNumber(eid);
  //     bool hasFound = false;
  //     for(int j = 0; j < starNum; j++){
  //       SimplexId starid;
  //       int res = triangulation_->getEdgeStar(eid, j, starid);
  //       if(res){
  //         std::cout << "[TestStellar] Cannot get edge star for edge id " << eid << ", error code: " << res << "\n";
  //         return -4;
  //       }
  //       if(starid == cid){
  //         hasFound = true;
  //         break;
  //       }
  //     }
  //     if(!hasFound){
  //       std::cout << "[TestStellar] Cannot find cell " << cid << " in edge " << eid << "\n";
  //       break;
  //     }
  //   }
  // }

  // // TF relation
  // for(SimplexId cid = 0; cid < cellNumber; cid++){
  //   for(int i = 0; i < 4; i++){
  //     SimplexId tid;
  //     triangulation_->getCellTriangle(cid, i, tid);
  //     SimplexId starNum = triangulation_->getTriangleStarNumber(tid);
  //     bool hasFound = false;
  //     for(int j = 0; j < starNum; j++){
  //       SimplexId starid;
  //       int res = triangulation_->getTriangleStar(tid, j, starid);
  //       if(res){
  //         std::cout << "[TestStellar] Cannot get traingle star for traingle id " << tid << ", error code: " << res << "\n";
  //         return -4;
  //       }
  //       if(starid == cid){
  //         hasFound = true;
  //         break;
  //       }
  //     }
  //     if(!hasFound){
  //       std::cout << "[TestStellar] Cannot find cell " << cid << " in traingle " << tid << "\n";
  //       break;
  //     }
  //   }
  // }


  /* test vertex edge relationship */
  t.reStart();
  for(SimplexId vertexId = 0; vertexId < vertexNumber; vertexId++){
    SimplexId edgeNum = triangulation_->getVertexEdgeNumber(vertexId);
    SimplexId edgeId;
    for(SimplexId j = 0; j < edgeNum; j++){
      int result1 = triangulation_->getVertexEdge(vertexId, j, edgeId);
      if(result1){
        std::cerr << "[TestStellar] vertexId " << vertexId << " Something wrong in getVertexEdge()! Error code: " << result1 << "\n";
      }
    }
  }
  std::cout << "[TestStellar] Time usage for VE: " << t.getElapsedTime() << " s.\n";

  /* test vertex triangle relationship */
  t.reStart();
  for(SimplexId vertexId = 0; vertexId < vertexNumber; vertexId++){
    SimplexId triangleNum = triangulation_->getVertexTriangleNumber(vertexId);
    SimplexId triangleId;
    for(SimplexId j = 0; j < triangleNum; j++){
      int result1 = triangulation_->getVertexTriangle(vertexId, j, triangleId);
      if(result1){
        std::cerr << "[TestStellar] vertexId " << vertexId << " Something wrong in getVertexTriangle()! Error code: " << result1 << "\n";
      }
    }
  }
  std::cout << "[TestStellar] Time usage for VT: " << t.getElapsedTime() << " s.\n";
  
  /* test vertex star relationship */
  t.reStart();
  for(SimplexId vertexId = 0; vertexId < vertexNumber; vertexId++){
    SimplexId cellNum = triangulation_->getVertexStarNumber(vertexId);
    SimplexId cellId;
    for(SimplexId j = 0; j < cellNum; j++){
      int result1 = triangulation_->getVertexStar(vertexId, j, cellId);
      if(result1){
        std::cout << "[TestStellar] vertexId " << vertexId << " Something wrong in getVertexStar()! Error code: " << result1 << "\n";
      }
    }
  }
  std::cout << "[TestStellar] Time usage for VS: " << t.getElapsedTime() << " s.\n";

  /* test edge vertex relationship */
  t.reStart();
  for(SimplexId edgeId = 0; edgeId < edgeNumber; edgeId++){
    SimplexId vertexId;
    for(SimplexId j = 0; j < 2; j++){
      int result1 = triangulation_->getEdgeVertex(edgeId, j, vertexId);
      if(result1){
        std::cout << "[TestStellar] edgeId " << edgeId << " Something wrong in getEdgeVertex()! Error code: " << result1 << "\n";
      }
    }
  }
  std::cout << "[TestStellar] Time usage for EV: " << t.getElapsedTime() << " s.\n";


  /* test edge triangle relationship */
  t.reStart();
  for(SimplexId edgeId = 0; edgeId < edgeNumber; edgeId++){
    SimplexId triangleNum = triangulation_->getEdgeTriangleNumber(edgeId);
    SimplexId triangleId;
    for(SimplexId j = 0; j < triangleNum; j++){
      int result1 = triangulation_->getEdgeTriangle(edgeId, j, triangleId);
      if(result1){
        std::cout << "[TestStellar] edgeId " << edgeId << " Something wrong in getEdgeTriangle()! Error code: " << result1 << "\n";
      }
    }
  }
  std::cout << "[TestStellar] Time usage for ET: " << t.getElapsedTime() << " s.\n";
  
  /* test edge star relationship */
  t.reStart();
  for(SimplexId edgeId = 0; edgeId < edgeNumber; edgeId++){
    SimplexId cellNum = triangulation_->getEdgeStarNumber(edgeId);
    SimplexId cellId;
    for(SimplexId j = 0; j < cellNum; j++){
      int result1 = triangulation_->getEdgeStar(edgeId, j, cellId);
      if(result1){
        std::cout << "[TestStellar] edgeId " << edgeId << " Something wrong in getEdgeStar()! Error code: " << result1 << "\n";
      }
    }
  }
  std::cout << "[TestStellar] Time usage for ES: " << t.getElapsedTime() << " s.\n";

  /* test triangle vertex relationship */
  t.reStart();
  for(SimplexId triangleId = 0; triangleId < triangleNumber; triangleId++){
    SimplexId vertexId;
    for(SimplexId j = 0; j < 3; j++){
      int result1 = triangulation_->getTriangleVertex(triangleId, j, vertexId);
      if(result1){
        std::cout << "[TestStellar] triangleId " << triangleId << " Something wrong in getTriangleVertex()! Error code: " << result1 << "\n";
      }
    }
  }
  std::cout << "[TestStellar] Time usage for TV: " << t.getElapsedTime() << " s.\n";

  /* test triangle edge relationship */
  t.reStart();
  for(SimplexId triangleId = 0; triangleId < triangleNumber; triangleId++){
    SimplexId edgeId;
    for(SimplexId j = 0; j < 3; j++){
      int result1 = triangulation_->getTriangleEdge(triangleId, j, edgeId);
      if(result1){
        std::cout << "[TestStellar] triangleId " << triangleId << " Something wrong in getTriangleEdge()! Error code: " << result1 << "\n";
      }
    }
  }
  std::cout << "[TestStellar] Time usage for TE: " << t.getElapsedTime() << " s.\n";

  /* test triangle star relationship */
  t.reStart();
  for(SimplexId triangleId = 0; triangleId < edgeNumber; triangleId++){
    SimplexId cellNum = triangulation_->getTriangleStarNumber(triangleId);
    SimplexId cellId;
    for(SimplexId j = 0; j < cellNum; j++){
      int result1 = triangulation_->getTriangleStar(triangleId, j, cellId);
      if(result1){
        std::cout << "[TestStellar] triangleId " << triangleId << " Something wrong in getTriangleStar()! Error code: " << result1 << "\n";
      }
    }
  }
  std::cout << "[TestStellar] Time usage for TS: " << t.getElapsedTime() << " s.\n";

  /* test cell vertex relationship */
  t.reStart();
  for(SimplexId cellId = 0; cellId < cellNumber; cellId++){
    SimplexId vertexId;
    for(SimplexId j = 0; j < verticesPerCell; j++){
      int result1 = triangulation_->getCellVertex(cellId, j, vertexId);
      if(result1){
        std::cout << "[TestStellar] cellId " << cellId << " Something wrong in getCellVertex()! Error code: " << result1 << "\n";
      }
    }
  }
  std::cout << "[TestStellar] Time usage for CV: " << t.getElapsedTime() << " s.\n";

  /* test cell edge relationship */
  t.reStart();
  for(SimplexId cellId = 0; cellId < cellNumber; cellId++){
    SimplexId edgeId;
    for(SimplexId j = 0; j < edgesPerCell; j++){
      int result1 = triangulation_->getCellEdge(cellId, j, edgeId);
      if(result1){
        std::cout << "[TestStellar] cellId " << cellId << " Something wrong in getCellEdge()! Error code: " << result1 << "\n";
      }
    }
  }
  std::cout << "[TestStellar] Time usage for CE: " << t.getElapsedTime() << " s.\n";

  /* test cell triangle relationship */
  t.reStart();
  for(SimplexId cellId = 0; cellId < cellNumber; cellId++){
    SimplexId triangleId;
    for(SimplexId j = 0; j < trianglesPerCell; j++){
      int result1 = triangulation_->getCellTriangle(cellId, j, triangleId);
      if(result1){
        std::cout << "[TestStellar] cellId " << cellId << " Something wrong in getCellTriangle()! Error code: " << result1 << "\n";
      }
    }
  }
  std::cout << "[TestStellar] Time usage for CT: " << t.getElapsedTime() << " s.\n";

  // init the output -- to adapt
  for(SimplexId i = 0; i < vertexNumber; i++){
    outputData[i] = inputData[i];
  }
  
  {
    std::stringstream msg;
    msg << "[TestStellar] Data-set (" << vertexNumber
      << " points) processed in "
      << tot.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }
  
  return 0;
}
