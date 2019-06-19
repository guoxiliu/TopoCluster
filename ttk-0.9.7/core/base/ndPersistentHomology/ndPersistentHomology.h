/// \ingroup base
/// \class ttk::ndPersistentHomology 
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK %ndPersistentHomology processing package.
///
/// %ndPersistentHomology is a TTK processing package that takes a scalar field on the input 
/// and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa ttkndPersistentHomology.cpp %for a usage example.

#pragma once

// base code includes
#include                  <Triangulation.h>
#include                  <Wrapper.h>

#include                  <PersistentHomology.h>

namespace ttk{
  
  class ndPersistentHomology : public Debug{

    public:
        
      ndPersistentHomology();
      
      ~ndPersistentHomology();

      template <class dataType>
        int execute(double point, double angle);
    
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
// #include                  <ndPersistentHomology.cpp>

// template functions
template <class dataType> int ttk::ndPersistentHomology::execute(double point, double angle){

  Timer t;
  
  // check the consistency of the variables -- to adapt
#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation_)
    return -1;
  if(!inputData_)
    return -2;
#endif

  cout << point << " " << angle << endl;

  vector<dataType*>* inputData = (vector<dataType*> *) inputData_;
  SimplexId vertexNumber = triangulation_->getNumberOfVertices();
  
  dataType min1,max1,min2,max2;
  min1 = max1 = inputData->at(0)[0];
  min2 = max2 = inputData->at(1)[0];

  for(int j=0; j<vertexNumber; j++){
    dataType val = inputData->at(0)[j];
    if(val > max1) max1 = val;
    if(val < min1) min1 = val;

    val = inputData->at(1)[j]; 
    if(val > max2) max2 = val;
    if(val < min2) min2 = val;
  }

  vector<dataType> b= {0,0};
  b[0] = ((max1+min1)*point);
  b[1] = ((max2+min2)*point);

  vector<double> m = {cos(angle*3.14/180),sin(angle*3.14/180)};

  dataType* oneFiltration = (dataType*) outputData;

  for(int j=0; j<vertexNumber; j++){
    dataType maxVal = (inputData->at(0)[j]-b[0])/m[0];
    for(int i=1; i<inputData_->size(); i++){
      if(maxVal < (inputData->at(i)[j]-b[i])/m[i])
        maxVal = (inputData->at(i)[j]-b[i])/m[i];
    }

    oneFiltration[j] = m[0]*maxVal;
    for(int i=1; i<inputData_->size(); i++){
      if(oneFiltration[j] > m[i]*maxVal)
        oneFiltration[j] = m[i]*maxVal;  
    }
  }

  ph = new PersistentHomology();
  ph->setupTriangulation(triangulation_);
  ph->setInputDataPointer(oneFiltration);
  ph->execute<dataType>(true);

  cout << " Got here" << endl;
  
  return 0;
}
