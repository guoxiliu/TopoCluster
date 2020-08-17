/// \ingroup base
/// \class ttk::FormangGradient 
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK %formangGradient processing package.
///
/// %FormangGradient is a TTK processing package that takes a scalar field on the input 
/// and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa ttkFormangGradient.cpp %for a usage example.

#pragma once

// base code includes
#include  <Triangulation.h>
#include  <Wrapper.h>
#include <unordered_set>

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/functional/hash.hpp>
#include <boost/optional/optional.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/dynamic_bitset.hpp>

using namespace std;
using namespace ttk;

typedef pair<short int,SimplexId> Simplex;
typedef set<Simplex, boost::function<bool(const Simplex &, const Simplex &)>> SimplexesSet;

namespace ttk{
  
  class FormanGradient : public Debug{

    public:
        
      FormanGradient();
      
      ~FormanGradient();


      void computeGradient(vector<int>* indexing);

      inline int setupTriangulation(Triangulation *triangulation){
        triangulation_ = triangulation;
        dimensionality_ = triangulation_->getCellVertexNumber(0) - 1;

        if(triangulation_){
            triangulation_->preprocessVertexStars();
            triangulation_->preprocessVertexEdges();
            triangulation_->preprocessEdges();
            triangulation_->preprocessEdgeStars();

            if (dimensionality_ >= 2) {
              triangulation_->preprocessTriangles();
              triangulation_->preprocessEdgeTriangles();
              triangulation_->preprocessTriangleEdges();
            }

            if (dimensionality_ == 3) {
              triangulation_->preprocessCellTriangles();
              triangulation_->preprocessTriangleStars();

              triangulation_->preprocessVertexTriangles();
            }

        }
        
        return 0;
      }

      SimplexId extractBoundary(const Simplex& simpl,int dimension_b, int index_b);
      SimplexId extractCoboundary(const Simplex& simpl, int dimension_cb, int index_cb);
      SimplexId sizeCofacets(const Simplex& simpl, int dimension_cb);
      SimplexId indexInsideCoface(const Simplex& simpl, const Simplex& coface);

      void simplexToVertices(Simplex simplex, vector<SimplexId>& vertices);

    protected:
      int homotopyExpansion(int vertex);

      int getIndex(Simplex simplex);

      void setPair(const Simplex& tail, const Simplex& head);
      bool getPair(const Simplex& simpl, Simplex& paired);
      bool getPairLower(const Simplex& simpl, Simplex& paired);
      bool getPairHigher(const Simplex& simpl, Simplex& paired);
      bool isCritical(const Simplex& simpl);

      int numPairableFaces(const Simplex& simplex, const map<SimplexId, int>& critical, SimplexId &chosen);

      bool cmpSimplexesLevelSets(const Simplex& s1, const Simplex& s2);
      void getVertexIndexingOfVertices(Simplex simplex, vector<SimplexId>& vec);
      
      void computeBarycenter(Simplex& simplex, vector<float>& coords);
    
    protected:
    
      vector<SimplexId>           *indexing_;
      Triangulation         *triangulation_;
      SimplexId                   dimensionality_;

      vector<boost::dynamic_bitset<>* > *gradient_;
      vector<set<SimplexId> > criticalSimplices;
  };
}
