/// \ingroup base
/// \class ttk::mapFunctionTo3D 
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK %mapFunctionTo3D processing package.
///
/// %mapFunctionTo3D is a TTK processing package that takes a scalar field on the input 
/// and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa ttkmapFunctionTo3D.cpp %for a usage example.

#pragma once

// base code includes
#include                  <Triangulation.h>
#include                  <Wrapper.h>


namespace ttk{
  
  class mapFunctionTo3D : public Debug{

    public:
        
      mapFunctionTo3D();
      
      ~mapFunctionTo3D();

      /// Execute the package.
      /// \pre If this TTK package uses ttk::Triangulation for fast mesh 
      /// traversals, the function setupTriangulation() must be called on this 
      /// object prior to this function, in a clearly distinct pre-processing 
      /// steps. An error will be returned otherwise.
      /// \note In such a case, it is recommended to exclude 
      /// setupTriangulation() from any time performance measurement.
      /// \param argment Dummy integer argument.
      /// \return Returns 0 upon success, negative values otherwise.
      template <class dataType>
        int execute(const int &argument) const;
    
      /// Pass a pointer to an input array representing a scalarfield.
      /// The expected format for the array is the following:
      /// <vertex0-component0> <vertex0-component1> ... <vertex0-componentN>
      /// <vertex1-component0> <vertex1-component1> ... <vertex1-componentN>
      /// <vertexM-component0> <vertexM-component1> ... <vertexM-componentN>.
      /// The array is expected to be correctly allocated. 
      /// \param data Pointer to the data array.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa setVertexNumber() and setDimensionNumber().
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
    
      void                  *inputData_, *outputData_;
      Triangulation         *triangulation_;
  };
}

// if the package is a pure template class, uncomment the following line
// #include                  <mapFunctionTo3D.cpp>

// template functions
template <class dataType> int ttk::mapFunctionTo3D::execute(
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

  // init the output -- to adapt
  for(SimplexId i = 0; i < vertexNumber; i++){
    outputData[i] = inputData[i];
  }
  
  // the following open-mp processing is only relevant for embarrassingly 
  // parallel algorithms (such as smoothing) -- to adapt
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) 
#endif
  for(SimplexId i = 0; i < vertexNumber; i++){
    
  }
   
  
  return 0;
}
