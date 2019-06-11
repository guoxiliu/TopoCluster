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
#include <utility>

using namespace std;



namespace ttk{
  
  class PersistentHomology : public Debug{

    public:
        
      PersistentHomology();
      
      ~PersistentHomology();

      template <class dataType>
        int execute(bool onVertices);

      void readHomology(vector<float>&, vector<char>&);

      template <class dataType>
      void readPersistencePairs(vector<float>&, vector<char>&, vector<dataType>&);
      
      void computeBarycenter(Simplex& simplex, vector<float>& coords);

      void indexToSimpl(int index, Simplex& simpl);
      void simplexToVertices(Simplex& simplex, vector<SimplexId>& vertices);

      template <class dataType>
      dataType getFiltration(Simplex& simplex); //return the value for the original function values

      int getFiltration(Simplex& simplex,const vector<int>& filtration); //return the value for the discrete filtration
    
      inline int setInputDataPointer(const void *data){
        inputData_ = (const void*)data;
        return 0;
      }

      inline int setupTriangulation(Triangulation *triangulation){
        triangulation_ = triangulation;
        dimensionality_ = triangulation_->getCellVertexNumber(0) - 1;
        
        if(triangulation_){
            triangulation_->preprocessEdges();

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
    
      const void                  *inputData_;
      Triangulation         *triangulation_;
      int                   dimensionality_;
      
      vector<vector<int> >        matrix;
      map<int,int>                allpairs;
      vector<int>                 homology;
  };
}

#include                  <persistenthomology_template.h>

// template functions
