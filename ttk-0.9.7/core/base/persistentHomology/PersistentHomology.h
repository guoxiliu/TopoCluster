/// \ingroup base
/// \class ttk::PersistentHomology 
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK %persistentHomology processing package.
///
/// %PersistentHomology is a TTK processing package that takes a scalar field on the input 
/// and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa ttkPersistentHomology.cpp %for a usage example.

#pragma once

// base code includes
#include                  <Triangulation.h>
#include                  <Wrapper.h>

#include <boundarymatrix.h>

namespace ttk{
  
  class PersistentHomology : public Debug{

    public:
        
      PersistentHomology();
      
      ~PersistentHomology();

      template <class dataType>
        int execute() const;
    
      inline int setInputDataPointer(void *data){
        inputData_ = data;
        return 0;
      }

      inline int setupTriangulation(Triangulation *triangulation){
        triangulation_ = triangulation;
       
        if(triangulation_){
          
          triangulation_->preprocessVertexNeighbors();

          
        }
        
        return 0;
      }
    
    protected:
    
      void                  *inputData_;
      Triangulation         *triangulation_;
  };
}

// if the package is a pure template class, uncomment the following line
// #include                  <PersistentHomology.cpp>

// template functions
template <class dataType> int ttk::PersistentHomology::execute() const{
 
  // check the consistency of the variables -- to adapt
#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation_)
    return -1;
  if(!inputData_)
    return -2;
#endif

  dataType *inputData = (dataType *) inputData_;
  
  SimplexId vertexNumber = triangulation_->getNumberOfVertices();

  // // init the output -- to adapt
  // for(SimplexId i = 0; i < vertexNumber; i++){
  //   outputData[i] = inputData[i];
  // }
  
  // the following open-mp processing is only relevant for embarrassingly 
  // parallel algorithms (such as smoothing) -- to adapt
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) 
#endif
  for(SimplexId i = 0; i < vertexNumber; i++){
    // TODO-2
    // processing here!
    // end of TODO-2
  }
   

  
  return 0;
}
