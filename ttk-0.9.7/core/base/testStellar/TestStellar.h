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
