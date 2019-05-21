/// \ingroup base
/// \class ttk::MultifieldComparisonMeasure 
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK %multifieldComparisonMeasure processing package.
///
/// %MultifieldComparisonMeasure is a TTK processing package that takes a scalar field on the input 
/// and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa ttkMultifieldComparisonMeasure.cpp %for a usage example.

#pragma once

// base code includes
#include                  <Triangulation.h>
#include                  <Wrapper.h>
#include                  <Eigen/Eigenvalues>

using namespace Eigen;
using namespace std;

namespace ttk{
  
  class MultifieldComparisonMeasure : public Debug{

    public:
        
      MultifieldComparisonMeasure();
      
      ~MultifieldComparisonMeasure();

      /// Execute the package.
      /// \pre If this TTK package uses ttk::Triangulation for fast mesh 
      /// traversals, the function setupTriangulation() must be called on this 
      /// object prior to this function, in a clearly distinct pre-processing 
      /// steps. An error will be returned otherwise.
      /// \note In such a case, it is recommended to exclude 
      /// setupTriangulation() from any time performance measurement.
      /// \param argment Dummy integer argument.
      /// \return Returns 0 upon success, negative values otherwise.
      template <class dataType> int execute() const;
    
      /// Pass a pointer to an input array representing a scalarfield.
      /// The expected format for the array is the following:
      /// <vertex0-component0> <vertex0-component1> ... <vertex0-componentN>
      /// <vertex1-component0> <vertex1-component1> ... <vertex1-componentN>
      /// <vertexM-component0> <vertexM-component1> ... <vertexM-componentN>.
      /// The array is expected to be correctly allocated. 
      /// \param data Pointer to the data array.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa setVertexNumber() and setDimensionNumber().
      inline int setInputDataPointer(std::vector<void*> *data){
        inputData_ = data;
        return 0;
      }

      /// Pass a pointer to an output array representing a scalar field.
      /// The expected format for the array is the following:
      /// <vertex0-component0> <vertex0-component1> ... <vertex0-componentN>
      /// <vertex1-component0> <vertex1-component1> ... <vertex1-componentN>
      /// <vertexM-component0> <vertexM-component1> ... <vertexM-componentN>.
      /// The array is expected to be correctly allocated.
      /// \param data Pointer to the data array.
      /// \return Returns 0 upon success, negative values otherwise. 
      /// \sa setVertexNumber() and setDimensionNumber().
      inline int setOutputDataPointer(void *data){
        outputData_ = data;
        return 0;
      }
     
      // General documentation info:
      //
      /// Setup a (valid) triangulation object for this TTK base object.
      ///
      /// \pre This function should be called prior to any usage of this TTK 
      /// object, in a clearly distinct pre-processing step that involves no 
      /// traversal or computation at all. An error will be returned otherwise.
      ///
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement. Therefore, it is recommended to 
      /// call this function ONLY in the pre-processing steps of your program. 
      /// Note however, that your triangulation object must be valid when 
      /// calling this function (i.e. you should have filled it at this point, 
      /// see the setInput*() functions of ttk::Triangulation). See ttkMultifieldComparisonMeasure 
      /// for further examples.
      ///
      /// \param triangulation Pointer to a valid triangulation.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa ttk::Triangulation
      //
      //
      // Developer info:
      // ttk::Triangulation is a generic triangulation representation that 
      // enables fast mesh traversal, either on explicit triangulations (i.e.
      // tet-meshes) or implicit triangulations (i.e. low-memory footprint 
      // implicit triangulations obtained from regular grids).
      //
      // Not all TTK packages need such mesh traversal features. If your 
      // TTK package needs any mesh traversal procedure, we recommend to use 
      // ttk::Triangulation as described here.
      //
      // Each call to a traversal procedure of ttk::Triangulation 
      // must satisfy some pre-condition (see ttk::Triangulation for more 
      // details). Such pre-condition functions are typically called from this
      // function. 
      inline int setupTriangulation(Triangulation *triangulation){
        triangulation_ = triangulation;
       
        if(triangulation_){
          
          // TODO-1
          // Pre-condition functions.
          // Call all the required pre-condition functions here!
          // for example:
          triangulation_->preprocessVertexNeighbors();
          // end of TODO-1
          
        }
        
        return 0;
      }
    
    protected:
    
      std::vector<void *>*        inputData_;
      void                        *outputData_;
      Triangulation         *triangulation_;
  };
}

// if the package is a pure template class, uncomment the following line
// #include                  <MultifieldComparisonMeasure.cpp>

// template functions
template <class dataType> int ttk::MultifieldComparisonMeasure::execute() const{

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
  std::vector<dataType*>* inputData = (std::vector<dataType*>*) inputData_;
  
  SimplexId vertexNumber = triangulation_->getNumberOfVertices();

  std::cout << "Are we here? 4" << std::endl;

  // the following open-mp processing is only relevant for embarrassingly 
  // parallel algorithms (such as smoothing) -- to adapt
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) 
#endif
  for(SimplexId i = 0; i < vertexNumber; i++){

    float x,y,z;
    triangulation_->getVertexPoint(i,x,y,z);

    MatrixXd A = MatrixXd::Zero(2, 2);

    for(int n=0; n<triangulation_->getVertexNeighborNumber(i); n++){
      SimplexId neighb;
      triangulation_->getVertexNeighbor(i,n,neighb);

      float xi,yi,zi;
      triangulation_->getVertexPoint(neighb,xi,yi,zi);

      if(y == yi && xi > x){
        std::cout << "Modified y" << std::endl;
        A(0,0) = (inputData->at(0)[i] - inputData->at(0)[neighb])/2*(xi-x);
        A(0,1) = (inputData->at(1)[i] - inputData->at(1)[neighb])/2*(xi-x); 
      }

      if(x == xi && yi > y){
        std::cout << "Modified x" << std::endl;
        A(1,0) = (inputData->at(0)[i] - inputData->at(0)[neighb])/2*(yi-y);
        A(1,1) = (inputData->at(1)[i] - inputData->at(1)[neighb])/2*(yi-y);                  
      }
    }

	  MatrixXd B = A.transpose() * A;
	  
    double max=0;
    VectorXcd eivals = B.eigenvalues();
    for(int e=0; e<eivals.size(); e++){

      if (max < eivals(e).real() )
        max = eivals(e).real();
    }
    outputData[i] = max;

  }
   
  {
    std::stringstream msg;
    msg << "[MultifieldComparisonMeasure] Data-set (" << vertexNumber
      << " points) processed in "
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }
  
  return 0;
}


  