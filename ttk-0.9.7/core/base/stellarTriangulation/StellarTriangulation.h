/// \ingroup base
/// \class ttk::StellarTriangulation 
/// \author Guoxi Liu <guoxil@g.clemson.edu>
/// \date Nov. 2019.
///
/// \brief TTK %stellarTriangulation processing package.
///
/// %StellarTriangulation is a TTK processing package that takes a scalar field on the input 
/// and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa ttkStellarTriangulation.cpp %for a usage example.

#pragma once

// base code includes
#include                  <AbstractTriangulation.h>

namespace ttk{
  
  class StellarTriangulation : public AbstractTriangulation{

    public:

      StellarTriangulation();

      ~StellarTriangulation();


      inline int setInputCells(const SimplexId &cellNumber,
                            const LongSimplexId *cellArray){

        if(cellNumber_)
            clear();

        cellNumber_ = cellNumber;
        cellArray_ = cellArray;

        return 0;
      }

      inline int setInputPoints(const SimplexId &pointNumber, const void *pointSet, const int *indexArray, const bool &doublePrecision = false){

        if(vertexNumber_)
            clear();

        vertexNumber_ = pointNumber;
        pointSet_ = pointSet;
        indices = indexArray;
        doublePrecision_ = doublePrecision;
        return 0;
      }

      int getCellEdge(const SimplexId &cellId, 
        const int &localEdgeId, SimplexId &edgeId) const{
        return 0;
      }
        
      SimplexId getCellEdgeNumber(const SimplexId &cellId) const{
        return 0;
      }
      
      const std::vector<std::vector<SimplexId> > *getCellEdges(){
        std::vector<std::vector<SimplexId>> *dummyPointer= new std::vector<std::vector<SimplexId>>();
        return dummyPointer;
      }
      
      int getCellNeighbor(const SimplexId &cellId,
        const int &localNeighborId, SimplexId &neighborId) const{
        return 0;
      }
        
      SimplexId getCellNeighborNumber(const SimplexId &cellId) const{
        return 0;
      }
      
      const std::vector<std::vector<SimplexId> > *getCellNeighbors(){
        std::vector<std::vector<SimplexId>> *dummyPointer= new std::vector<std::vector<SimplexId>>();
        return dummyPointer;
      }
      
      int getCellTriangle(const SimplexId &cellId, 
        const int &localTriangleId, SimplexId &triangleId) const{
        return 0;
      }
        
      SimplexId getCellTriangleNumber(const SimplexId &cellId) const{
        return 0;
      }
        
      const std::vector<std::vector<SimplexId> > *getCellTriangles(){
        std::vector<std::vector<SimplexId>> *dummyPointer= new std::vector<std::vector<SimplexId>>();
        return dummyPointer;
      }
      
      int getCellVertex(const SimplexId &cellId,
        const int &localVertexId, SimplexId &vertexId) const{
        return 0;
      }
    
      SimplexId getCellVertexNumber(const SimplexId &cellId) const{
        return 0;
      }
        
      int getDimensionality() const{
        return 0;
      }
      
      const std::vector<std::pair<SimplexId, SimplexId> > *getEdges(){
        std::vector<std::pair<SimplexId, SimplexId>> *dummyPointer= new std::vector<std::pair<SimplexId, SimplexId>>();
        return dummyPointer;
      }
        
      int getEdgeLink(const SimplexId &edgeId, 
        const int &localLinkId, SimplexId &linkId) const{
        return 0;
      }
        
      SimplexId getEdgeLinkNumber(const SimplexId &edgeId) const{
        return 0;
      }
      
      const std::vector<std::vector<SimplexId> > *getEdgeLinks(){
        std::vector<std::vector<SimplexId>> *dummyPointer= new std::vector<std::vector<SimplexId>>();
        return dummyPointer;
      }
      
      int getEdgeStar(const SimplexId &edgeId, 
        const int &localStarId, SimplexId &starId) const{
        return 0;
      }
        
      SimplexId getEdgeStarNumber(const SimplexId &edgeId) const{
        return 0;
      }
      
      const std::vector<std::vector<SimplexId> > *getEdgeStars(){
        std::vector<std::vector<SimplexId>> *dummyPointer= new std::vector<std::vector<SimplexId>>();
        return dummyPointer;
      }
     
      int getEdgeTriangle(const SimplexId &edgeId, 
        const int &localTriangleId, SimplexId &triangleId) const{
        return 0;
      }
        
      SimplexId getEdgeTriangleNumber(const SimplexId &edgeId) const{
        return 0;
      }
        
      const std::vector<std::vector<SimplexId> > *getEdgeTriangles(){
        std::vector<std::vector<SimplexId>> *dummyPointer= new std::vector<std::vector<SimplexId>>();
        return dummyPointer;
      }
      
      int getEdgeVertex(const SimplexId &edgeId, 
        const int &localVertexId, SimplexId &vertexId) const{
        return 0;
      }
      
      SimplexId getNumberOfCells() const{
        return 0;
      }
      
      SimplexId getNumberOfEdges() const{
        return 0;
      }
      
      SimplexId getNumberOfTriangles() const{
        return 0;
      }
      
      SimplexId getNumberOfVertices() const{
        return 0;
      }
      
      const std::vector<std::vector<SimplexId> > *getTriangles(){
        std::vector<std::vector<SimplexId>> *dummyPointer= new std::vector<std::vector<SimplexId>>();
        return dummyPointer;
      }
      
      int getTriangleEdge(const SimplexId &triangleId,
        const int &localEdgeId, SimplexId &edgeId) const{
        return 0;
      }
      
      SimplexId getTriangleEdgeNumber(const SimplexId &triangleId) const{
        return 0;
      }
      
      const std::vector<std::vector<SimplexId> > *getTriangleEdges(){
        std::vector<std::vector<SimplexId>> *dummyPointer= new std::vector<std::vector<SimplexId>>();
        return dummyPointer;
      }
      
      int getTriangleLink(const SimplexId &triangleId, 
        const int &localLinkId, SimplexId &linkId) const{
          return 0;
        }
        
      SimplexId getTriangleLinkNumber(const SimplexId &triangleId) const{
        return 0;
      }
      
      const std::vector<std::vector<SimplexId> > *getTriangleLinks(){
        std::vector<std::vector<SimplexId>> *dummyPointer= new std::vector<std::vector<SimplexId>>();
        return dummyPointer;
      }
      
      int getTriangleStar(const SimplexId &triangleId,
        const int &localStarId, SimplexId &starId) const{
        return 0;
      }
        
      SimplexId getTriangleStarNumber(const SimplexId &triangleId) const{
        return 0;
      }
      
      const std::vector<std::vector<SimplexId> > *getTriangleStars(){
        std::vector<std::vector<SimplexId>> *dummyPointer= new std::vector<std::vector<SimplexId>>();
        return dummyPointer;
      }
      
      int getTriangleVertex(const SimplexId &triangleId,
        const int &localVertexId, SimplexId &vertexId) const{
        return 0;
      }
      
      int getVertexEdge(const SimplexId &vertexId, 
        const int &localEdgeId, SimplexId &edgeId) const{
        return 0;
      }
        
      SimplexId getVertexEdgeNumber(const SimplexId &vertexId) const{
        return 0;
      }
      
      const std::vector<std::vector<SimplexId> > *getVertexEdges(){
        std::vector<std::vector<SimplexId>> *dummyPointer= new std::vector<std::vector<SimplexId>>();
        return dummyPointer;
      }
      
      int getVertexLink(const SimplexId &vertexId, 
        const int &localLinkId, SimplexId &linkId) const{
        return 0;
      }
        
      SimplexId getVertexLinkNumber(const SimplexId &vertexId) const{
        return 0;
      }
      
      const std::vector<std::vector<SimplexId> > *getVertexLinks(){
        std::vector<std::vector<SimplexId>> *dummyPointer= new std::vector<std::vector<SimplexId>>();
        return dummyPointer;
      }
      
      int getVertexNeighbor(const SimplexId &vertexId, 
        const int &localNeighborId, SimplexId &neighborId) const{
        return 0;
      }
        
      SimplexId getVertexNeighborNumber(const SimplexId &vertexId) const{
        return 0;
      }
      
      const std::vector<std::vector<SimplexId> > *getVertexNeighbors(){
        std::vector<std::vector<SimplexId>> *dummyPointer= new std::vector<std::vector<SimplexId>>();
        return dummyPointer;
      }
      
      int getVertexPoint(const SimplexId &vertexId,
        float &x, float &y, float &z) const{
        return 0;
      }
        
      int getVertexStar(const SimplexId &vertexId, const int &localStarId,
        SimplexId &starId) const{
        return 0;
      }
        
      SimplexId getVertexStarNumber(const SimplexId &vertexId) const{
        return 0;
      }
        
      const std::vector<std::vector<SimplexId> > *getVertexStars(){
        std::vector<std::vector<SimplexId>> *dummyPointer= new std::vector<std::vector<SimplexId>>();
        return dummyPointer;
      }
      
      int getVertexTriangle(const SimplexId &vertexId, 
        const int &localTriangleId, SimplexId &triangleId) const{
        return 0;
      }
        
      SimplexId getVertexTriangleNumber(const SimplexId &vertexId) const{
        return 0;
      }
        
      const std::vector<std::vector<SimplexId> > *getVertexTriangles(){
        std::vector<std::vector<SimplexId>> *dummyPointer= new std::vector<std::vector<SimplexId>>();
        return dummyPointer;
      }
        
      inline bool hasPreprocessedBoundaryEdges() const{
        return hasPreprocessedBoundaryEdges_;
      }
      
      inline bool hasPreprocessedBoundaryTriangles() const{
        return hasPreprocessedBoundaryTriangles_;
      }
       
      inline bool hasPreprocessedBoundaryVertices() const{
        return hasPreprocessedBoundaryVertices_;
      }
     
      inline bool hasPreprocessedCellEdges() const{
        return hasPreprocessedCellEdges_;
      }
      
      inline bool hasPreprocessedCellNeighbors() const{
        return hasPreprocessedCellNeighbors_;
      }
      
      inline bool hasPreprocessedCellTriangles() const{
        return hasPreprocessedCellTriangles_;
      }
       
      inline bool hasPreprocessedEdgeLinks() const{
        return hasPreprocessedEdgeLinks_;
      }
       
      inline bool hasPreprocessedEdgeStars() const{
        return hasPreprocessedEdgeStars_;
      }
      
      inline bool hasPreprocessedEdgeTriangles() const{
        return hasPreprocessedEdgeTriangles_;
      }
        
      inline bool hasPreprocessedEdges() const{
        return hasPreprocessedEdges_;
      }
      
      inline bool hasPreprocessedTriangles() const{
        return hasPreprocessedTriangles_;
      }
      
      inline bool hasPreprocessedTriangleEdges() const{
        return hasPreprocessedTriangleEdges_;
      }
      
      inline bool hasPreprocessedTriangleLinks() const{
        return hasPreprocessedTriangleLinks_;
      }
      
      inline bool hasPreprocessedTriangleStars() const{
        return hasPreprocessedTriangleStars_;
      }
        
      inline bool hasPreprocessedVertexEdges() const{
        return hasPreprocessedVertexEdges_;
      }
      
      inline bool hasPreprocessedVertexLinks() const{
        return hasPreprocessedVertexLinks_;
      }
      
      inline bool hasPreprocessedVertexNeighbors() const{
        return hasPreprocessedVertexNeighbors_;
      }
      
      inline bool hasPreprocessedVertexStars() const{
        return hasPreprocessedVertexStars_;
      }
      
      inline bool hasPreprocessedVertexTriangles() const{
        return hasPreprocessedVertexTriangles_;
      }
      
      bool isEdgeOnBoundary(const SimplexId &edgeId) const{
        return false;
      }
        
      bool isEmpty() const{
        return false;
      }
      
      bool isTriangleOnBoundary(const SimplexId &triangleId) const{
        return false;
      }
      
      bool isVertexOnBoundary(const SimplexId &vertexId) const{
        return false;
      }

      int preprocessBoundaryEdges(){
        preprocessEdges();
        hasPreprocessedBoundaryEdges_ = true;
        return 0;
      }
      
      int preprocessBoundaryTriangles(){
        preprocessTriangles();
        hasPreprocessedBoundaryTriangles_ = true;
        return 0;
      }
      
      int preprocessBoundaryVertices(){
        hasPreprocessedBoundaryVertices_ = true;
        return 0;
      }
      
      int preprocessCellEdges(){
        preprocessEdges();
        hasPreprocessedCellEdges_ = true;
        return 0;
      }
      
      int preprocessCellNeighbors(){
        hasPreprocessedCellNeighbors_ = true;
        return 0;
      }
      
      int preprocessCellTriangles(){
        preprocessTriangles();
        hasPreprocessedCellTriangles_ = true;
        return 0;
      }
      
      int preprocessEdges(){
        hasPreprocessedEdges_ = true;
        return 0;
      }
      
      int preprocessEdgeLinks(){
        preprocessEdges();
        hasPreprocessedEdgeLinks_ = true;
        return 0;
      }
      
      int preprocessEdgeStars(){
        preprocessEdges();
        hasPreprocessedEdgeStars_ = true;
        return 0;
      }
      
      int preprocessEdgeTriangles(){
        preprocessEdges();
        preprocessTriangles();
        hasPreprocessedEdgeTriangles_ = true;
        return 0;
      }
      
      int preprocessTriangles(){
        hasPreprocessedTriangles_ = true;
        return 0;
      }
      
      int preprocessTriangleEdges(){
        preprocessEdges();
        preprocessTriangles();
        hasPreprocessedTriangleEdges_ = true;
        return 0;
      }
      
      int preprocessTriangleLinks(){
        preprocessTriangles();
        hasPreprocessedTriangleLinks_ = true;
        return 0;
      }
      
      int preprocessTriangleStars(){
        preprocessTriangles();
        hasPreprocessedTriangleStars_ = true;
        return 0;
      }
      
      int preprocessVertexEdges() { 
        preprocessEdges();
        hasPreprocessedVertexEdges_ = true;
        return 0;
      }
      
      int preprocessVertexLinks(){
        hasPreprocessedVertexLinks_ = true;
        return 0;
      }
      
      int preprocessVertexNeighbors(){
        hasPreprocessedVertexNeighbors_ = true;
        return 0;
      }
      
      int preprocessVertexStars(){
        hasPreprocessedVertexStars_ = true;
        return 0;
      }
      
      int preprocessVertexTriangles(){
        preprocessTriangles();
        hasPreprocessedVertexTriangles_ = true;
        return 0;
      }

    protected:
      
      // empty wrapping to VTK for now
      bool needsToAbort(){ return false;};
      
      template <class itemType>
        size_t tableFootprint(const std::vector<itemType> &table,
          const std::string tableName = "", 
          std::stringstream *msg = NULL) const{
        
        if((table.size())&&(tableName.length())&&(msg)){
          (*msg) << "[AbstractTriangulation] " << tableName << ": "
            << table.size()*sizeof(itemType) << " bytes" << std::endl;
        }
            
        return table.size()*sizeof(itemType);
      }
      
      template <class itemType>
        size_t tableTableFootprint(
          const std::vector<std::vector<itemType> > &table,
          const std::string tableName = "", 
          std::stringstream *msg = NULL) const;
      
      int updateProgress(const float &progress) {return 0;};
      
      bool                hasPreprocessedBoundaryEdges_,
                          hasPreprocessedBoundaryTriangles_,
                          hasPreprocessedBoundaryVertices_,
                          hasPreprocessedCellEdges_,
                          hasPreprocessedCellNeighbors_,
                          hasPreprocessedCellTriangles_,
                          hasPreprocessedEdges_,
                          hasPreprocessedEdgeLinks_,
                          hasPreprocessedEdgeStars_,
                          hasPreprocessedEdgeTriangles_,
                          hasPreprocessedTriangles_,
                          hasPreprocessedTriangleEdges_,
                          hasPreprocessedTriangleLinks_,
                          hasPreprocessedTriangleStars_,
                          hasPreprocessedVertexEdges_,
                          hasPreprocessedVertexLinks_,
                          hasPreprocessedVertexNeighbors_,
                          hasPreprocessedVertexStars_,
                          hasPreprocessedVertexTriangles_;
      
      std::vector<bool>   boundaryEdges_,
                          boundaryTriangles_,
                          boundaryVertices_;
      std::vector<std::vector<SimplexId> > 
                          cellEdgeList_;
      std::vector<std::vector<SimplexId> >
                          cellNeighborList_;
      std::vector<std::vector<SimplexId> > 
                          cellTriangleList_;
      std::vector<std::vector<SimplexId> >
                          edgeLinkList_;
      std::vector<std::pair<SimplexId, SimplexId> >
                          edgeList_;
      std::vector<std::vector<SimplexId> >
                          edgeStarList_;
      std::vector<std::vector<SimplexId> >
                          edgeTriangleList_;
      std::vector<std::vector<SimplexId> >
                          triangleList_;
      std::vector<std::vector<SimplexId> > 
                          triangleEdgeList_;
      std::vector<std::vector<SimplexId> >
                          triangleLinkList_;
      std::vector<std::vector<SimplexId> >
                          triangleStarList_;
      std::vector<std::vector<SimplexId> > 
                          vertexEdgeList_;
      std::vector<std::vector<SimplexId> >
                          vertexLinkList_;
      std::vector<std::vector<SimplexId> > 
                          vertexNeighborList_;
      std::vector<std::vector<SimplexId> >
                          vertexStarList_;
      std::vector<std::vector<SimplexId> >
                          vertexTriangleList_;


  protected:

      int clear();

      bool                doublePrecision_;
      SimplexId           cellNumber_, vertexNumber_;
      const void          *pointSet_;
      const int           *indices;
      const LongSimplexId *cellArray_;

  };
}
