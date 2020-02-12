/// \ingroup base
/// \class ttk::PreprocessStellar 
/// \author Guoxi Liu <guoxil@g.clemson.edu>
/// \date Nov. 2019.


#pragma once

// base code includes
#include <Triangulation.h>
#include <Wrapper.h>
#include <Octree.h>


namespace ttk{
  
  class PreprocessStellar : public Debug{

    public:
        
      PreprocessStellar();
      
      ~PreprocessStellar();

      int execute(const int &argument) const;
    
      inline int setVerticesPointer(void *data){
        vertices = data;
        return 0;
      }

      inline int setNodesPointer(void *data){
        nodes = data;
        return 0;
      }

      inline int setCellsPointer(void *data){
        cells = data;
        return 0;
      }
     
      /// Setup a (valid) triangulation object for this TTK base object.
      inline int setupTriangulation(Triangulation *triangulation){
        triangulation_ = triangulation;
        return 0;
      }
    
    protected:
    
      void                  *vertices, *nodes, *cells;
      Triangulation         *triangulation_;
  };
}
