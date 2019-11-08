/// \ingroup base
/// \class ttk::PreprocessStellar 
/// \author Guoxi Liu <guoxil@g.clemson.edu>
/// \date Nov. 2019.


#pragma once

// base code includes
#include                  <Triangulation.h>
#include                  <Wrapper.h>
#include                  <Octree.h>


namespace ttk{
  
  class PreprocessStellar : public Debug{

    public:
        
      PreprocessStellar();
      
      ~PreprocessStellar();

      /// Execute the package.
      template <class dataType>
        int execute(const int &argument) const;
    
      /// Pass a pointer to an input array representing a scalarfield.
      inline int setInputDataPointer(void *data){
        inputData_ = data;
        return 0;
      }

      /// Pass a pointer to an output array representing a scalar field.
      inline int setOutputDataPointer(void *data){
        outputData_ = data;
        return 0;
      }
     
      /// Setup a (valid) triangulation object for this TTK base object.
      inline int setupTriangulation(Triangulation *triangulation){
        triangulation_ = triangulation;
       
        if(triangulation_){
          // Pre-condition functions.
          triangulation_->preprocessVertexNeighbors();
        }
        
        return 0;
      }
    
    protected:
    
      void                  *inputData_, *outputData_;
      Triangulation         *triangulation_;
  };
}


// template functions
template <class dataType> int ttk::PreprocessStellar::execute(
  const int &argument) const{

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
  SimplexId cellNumber = triangulation_->getNumberOfCells();

  cout << "Number of vertices: " << vertexNumber << endl;
  // init the output -- to adapt
  for(SimplexId i = 0; i < vertexNumber; i++){
    outputData[i] = inputData[i];
  }
  
  // the following open-mp processing is only relevant for embarrassingly 
  // parallel algorithms (such as smoothing) -- to adapt
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) 
#endif

  Octree preOctree(triangulation_, 1000);
  for(SimplexId i = 0; i < vertexNumber; i++){
    preOctree.insertVertex(i);
  }
  for(SimplexId i = 0; i < cellNumber; i++){
    preOctree.insertCell(i);
  }


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
