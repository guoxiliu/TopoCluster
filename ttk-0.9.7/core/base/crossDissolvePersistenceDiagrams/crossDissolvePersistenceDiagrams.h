/// \ingroup base
/// \class ttk::crossDissolvePersistenceDiagrams 
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK %crossDissolvePersistenceDiagrams processing package.
///
/// %crossDissolvePersistenceDiagrams is a TTK processing package that takes a scalar field on the input 
/// and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa ttkcrossDissolvePersistenceDiagrams.cpp %for a usage example.

#pragma once

// base code includes
#include                  <Triangulation.h>
#include                  <Wrapper.h>

#include                  <PersistentHomology.h>

namespace ttk{
  
  class crossDissolvePersistenceDiagrams : public Debug{

    public:
        
      crossDissolvePersistenceDiagrams();
      
      ~crossDissolvePersistenceDiagrams();

      template <class dataType>
        int execute(double alpha);
    
      inline int setInputDataPointer(vector<void *>* data){
        inputData_ = data;
        return 0;
      }

      inline void setOutputDataPointer(void* data){
        outputData = data;
      }

      inline PersistentHomology* getPH(){return ph;}

      inline int setupTriangulation(Triangulation *triangulation){
        triangulation_ = triangulation;
        dimensionality_ = triangulation_->getCellVertexNumber(0) - 1;

        if(triangulation_){
          
          triangulation_->preprocessVertexEdges();

          if (dimensionality_ >= 2) {
            triangulation_->preprocessTriangleEdges();
          }

          if (dimensionality_ == 3) {
            triangulation_->preprocessCellTriangles();
          }

        }
        
        return 0;
      }
    
    protected:
    
      vector<void*>         *inputData_;
      void                  *outputData;
      Triangulation         *triangulation_;

      int                   dimensionality_;
      PersistentHomology    *ph; 
  };
}

// if the package is a pure template class, uncomment the following line
// #include                  <crossDissolvePersistenceDiagrams.cpp>

// template functions
template <class dataType> int ttk::crossDissolvePersistenceDiagrams::execute(double alpha){

  Timer t;
  
  // check the consistency of the variables -- to adapt
#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation_)
    return -1;
  if(!inputData_)
    return -2;
#endif

  cout << alpha << endl;

  vector<dataType*>* inputData = (vector<dataType*> *) inputData_;
  SimplexId vertexNumber = triangulation_->getNumberOfVertices();

  //define first filtration with integers

  //define second filtration with integers

  //cross-dissolve based on alpha

  dataType* oneFiltration = (dataType*) outputData;
  for(int j=0; j<vertexNumber; j++){
    oneFiltration[j] = (1-alpha)*inputData->at(0)[j] + alpha*inputData->at(1)[j];
  }

  ph = new PersistentHomology();
  ph->setupTriangulation(triangulation_);
  ph->setInputDataPointer(oneFiltration);
  ph->execute<dataType>(true);

  cout << " Got here" << endl;
  
  return 0;
}
