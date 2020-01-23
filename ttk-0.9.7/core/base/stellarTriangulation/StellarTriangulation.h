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
#include                  <algorithm>
#include                  <AbstractTriangulation.h>

#define VERTEX_ID 0
#define EDGE_ID 1
#define TRIANGLE_ID 2
#define CELL_ID 3

using namespace std;
namespace ttk{
  
  class StellarTriangulation : public AbstractTriangulation{

    public:

      StellarTriangulation();

      ~StellarTriangulation();


      int setInputCells(const SimplexId &cellNumber,
                            const LongSimplexId *cellArray){

        if(cellNumber_)
            clear();

        cellNumber_ = cellNumber;
        cellArray_ = cellArray;

        if((!cellArray_)||(!cellNumber_))
          return -1;
        
        // initialize the array of cell intervals
        cellIntervals_.resize(nodeNumber_+1);
        externalCells_.resize(nodeNumber_+1);
        cellIntervals_[0] = -1;

        SimplexId nodeNum = 1, cid = 0;
        vector<SimplexId> cell;
        for(; cid < cellNumber; cid++){
          SimplexId startPos = (cellArray[0]+1)*cid+1;
          cell = vector<SimplexId>(cellArray+startPos, cellArray+startPos+cellArray[0]);

          if(cell[0] > vertexIntervals_[nodeNum]){
            cellIntervals_[nodeNum++] = cid - 1;
          }

          // create external cell list 
          for(size_t i = 0; i < cell.size(); i++){
            if(cell[i] > vertexIntervals_[nodeNum]){
              int nid = findNodeIndex(cell[i], VERTEX_ID);
              if(externalCells_[nid].empty()){
                externalCells_[nid].push_back(cid);
              }else if(externalCells_[nid].back() != cid){
                externalCells_[nid].push_back(cid);
              }
            }
          }
        }
        cellIntervals_.back() = cid-1;
        
        // TEST: print out the cell intervals
        cout << "Cell intervals: \n";
        for(size_t i = 0; i < cellIntervals_.size(); i++){
          cout << i << ": " << cellIntervals_[i] << endl;
        }

        return 0;
      }

      int setInputPoints(const SimplexId &pointNumber, const void *pointSet, 
        const int *indexArray, const bool &doublePrecision = false){

        if(vertexNumber_)
            clear();

        vertexNumber_ = pointNumber;
        pointSet_ = pointSet;
        doublePrecision_ = doublePrecision;

        // initialize the array of vertex intervals
        // add a dummy node at the beginning
        vertexIntervals_.push_back(-1);

        SimplexId vid = 1;
        for(; vid < pointNumber; vid++){
          if(indexArray[vid] != indexArray[vid-1]){
            vertexIntervals_.push_back(vid-1);
          }
        }
        vertexIntervals_.push_back(vid-1);
        nodeNumber_ = vertexIntervals_.size()-1;

        // TEST: print out the vertex intervals
        cout << "Vertex intervals:\n";
        for(size_t i = 0; i < vertexIntervals_.size(); i++){
          cout << i << ": " << vertexIntervals_[i] << endl;
        }

        return 0;
      }

      int getCellEdge(const SimplexId &cellId, 
        const int &localEdgeId, SimplexId &edgeId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if((cellId < 0)||(cellId >= cellNumber_))
            return -1;
          if((localEdgeId < 0))
            return -2;
        #endif

        SimplexId nid = findNodeIndex(cellId, CELL_ID);
        vector<vector<SimplexId>> localCellEdges;
        if(!getCellEdges(nid, localCellEdges)){
          SimplexId localCellId = cellId - cellIntervals_[nid-1] - 1;
          if(localEdgeId >= (SimplexId) localCellEdges[localCellId].size()){
            return -2;
          }
          edgeId = localCellEdges[localCellId][localEdgeId];
          return 0;
        }
        return -3;
      }
        
      inline SimplexId getCellEdgeNumber(const SimplexId &cellId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(vertexNumber_ <= 0)
            return -1;
          if(cellNumber_ <= 0)
            return -2;
          if(!cellArray_)
          return -3;
        #endif

        return cellArray_[0] * (cellArray_[0]-1) / 2;
      }
      
      const vector<vector<SimplexId> > *getCellEdges(){
        cellEdgeList_.reserve(vertexNumber_);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          vector<vector<SimplexId>> localCellEdges;
          getVertexEdges(nid, localCellEdges);
          cellEdgeList_.insert(cellEdgeList_.end(), localCellEdges.begin(), localCellEdges.end());
        }
        return &cellEdgeList_;
      }
      
      int getCellNeighbor(const SimplexId &cellId,
        const int &localNeighborId, SimplexId &neighborId) const{
        return 0;
      }
        
      SimplexId getCellNeighborNumber(const SimplexId &cellId) const{
        return 0;
      }
      
      const vector<vector<SimplexId> > *getCellNeighbors(){
        vector<vector<SimplexId>> *dummyPointer= new vector<vector<SimplexId>>();
        return dummyPointer;
      }
      
      int getCellTriangle(const SimplexId &cellId, 
        const int &localTriangleId, SimplexId &triangleId) const{
        
        return 0;
      }
        
      SimplexId getCellTriangleNumber(const SimplexId &cellId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(vertexNumber_ <= 0)
            return -1;
          if(cellNumber_ <= 0)
            return -2;
          if(!cellArray_)
            return -3;
        #endif

        return cellArray_[0] * (cellArray_[0]-1) * (cellArray_[0]-2) / 6;
      }
        
      const vector<vector<SimplexId> > *getCellTriangles(){
        vector<vector<SimplexId>> *dummyPointer= new vector<vector<SimplexId>>();
        return dummyPointer;
      }
      
      inline int getCellVertex(const SimplexId &cellId,
        const int &localVertexId, SimplexId &vertexId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if((cellId < 0)||(cellId >= cellNumber_))
            return -1;
          if((localVertexId < 0)
            ||(localVertexId >= cellArray_[0]))
            return -2;
        #endif

        vertexId = cellArray_[(cellArray_[0] + 1)*cellId + localVertexId + 1];
        return 0;
      }
    
      inline SimplexId getCellVertexNumber(const SimplexId &cellId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if((cellId < 0)||(cellId >= cellNumber_))
            return -1;
          if((!cellArray_)||(!cellNumber_))
            return -2;
        #endif
        
        return cellArray_[0];
      }
        
      inline int getDimensionality() const{
        if((cellArray_)&&(cellNumber_)){
          return cellArray_[0] - 1;
        }
        return -1;
      }
      
      const vector<pair<SimplexId, SimplexId> > *getEdges(){
        edgeList_.reserve(edgeIntervals_.back()+1);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          vector<pair<SimplexId, SimplexId>> localInternalEdges, localExternalEdges;
          buildEdgeList(nid, localInternalEdges, localExternalEdges);
          edgeList_.insert(edgeList_.end(), localInternalEdges.begin(), localInternalEdges.end());
        }
        return &edgeList_;
      }
        
      int getEdgeLink(const SimplexId &edgeId, 
        const int &localLinkId, SimplexId &linkId) const{
        return 0;
      }
        
      SimplexId getEdgeLinkNumber(const SimplexId &edgeId) const{
        return 0;
      }
      
      const vector<vector<SimplexId> > *getEdgeLinks(){
        vector<vector<SimplexId>> *dummyPointer= new vector<vector<SimplexId>>();
        return dummyPointer;
      }
      
      int getEdgeStar(const SimplexId &edgeId, 
        const int &localStarId, SimplexId &starId) const{
        return 0;
      }
        
      SimplexId getEdgeStarNumber(const SimplexId &edgeId) const{
        return 0;
      }
      
      const vector<vector<SimplexId> > *getEdgeStars(){
        vector<vector<SimplexId>> *dummyPointer= new vector<vector<SimplexId>>();
        return dummyPointer;
      }
     
      int getEdgeTriangle(const SimplexId &edgeId, 
        const int &localTriangleId, SimplexId &triangleId) const{
        return 0;
      }
        
      SimplexId getEdgeTriangleNumber(const SimplexId &edgeId) const{
        return 0;
      }
        
      const vector<vector<SimplexId> > *getEdgeTriangles(){
        vector<vector<SimplexId>> *dummyPointer= new vector<vector<SimplexId>>();
        return dummyPointer;
      }
      
      int getEdgeVertex(const SimplexId &edgeId, 
        const int &localVertexId, SimplexId &vertexId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if((edgeId < 0)||(edgeId > (SimplexId) edgeIntervals_.back()))
            return -1;
          if((localVertexId != 0)&&(localVertexId != 1))
            return -2;
          if(!hasPreprocessedEdges_)
            return -3;
        #endif

        SimplexId nid = findNodeIndex(edgeId, EDGE_ID);
        vector<pair<SimplexId, SimplexId>> localInternalEdges,localExternalEdges;
        buildEdgeList(nid, localInternalEdges, localExternalEdges);

        SimplexId localEdgeId = edgeId-edgeIntervals_[nid-1]-1;

        if(localVertexId){
          return localInternalEdges[localEdgeId].second;
        }
        return localInternalEdges[localEdgeId].first;
      }
      
      inline SimplexId getNumberOfCells() const { return cellNumber_; }
      
      inline SimplexId getNumberOfEdges() const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(!edgeIntervals_.size())
            return -1;
        #endif

        return edgeIntervals_.back()+1;
      }
      
      SimplexId getNumberOfTriangles() const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(!triangleIntervals_.size())
            return -1;
        #endif

        return triangleIntervals_.back()+1;
      }
      
      inline SimplexId getNumberOfVertices() const { return vertexNumber_; }
      
      const vector<vector<SimplexId> > *getTriangles(){
        vector<vector<SimplexId>> *dummyPointer= new vector<vector<SimplexId>>();
        return dummyPointer;
      }
      
      int getTriangleEdge(const SimplexId &triangleId,
        const int &localEdgeId, SimplexId &edgeId) const{
        return 0;
      }
      
      SimplexId getTriangleEdgeNumber(const SimplexId &triangleId) const{
        #ifndef TTK_ENABLE_KAMIKAZE
          if((triangleId < 0)||(triangleId > triangleIntervals_.back()))
            return -1;
        #endif
        
        return 3;
      }
      
      const vector<vector<SimplexId> > *getTriangleEdges(){
        vector<vector<SimplexId>> *dummyPointer= new vector<vector<SimplexId>>();
        return dummyPointer;
      }
      
      int getTriangleLink(const SimplexId &triangleId, 
        const int &localLinkId, SimplexId &linkId) const{
          return 0;
        }
        
      SimplexId getTriangleLinkNumber(const SimplexId &triangleId) const{
        return 0;
      }
      
      const vector<vector<SimplexId> > *getTriangleLinks(){
        vector<vector<SimplexId>> *dummyPointer= new vector<vector<SimplexId>>();
        return dummyPointer;
      }
      
      int getTriangleStar(const SimplexId &triangleId,
        const int &localStarId, SimplexId &starId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if((triangleId < 0)||(triangleId >= triangleIntervals_.back()))
            return -1;
          if(localStarId < 0)
            return -2;
        #endif

        SimplexId nid = findNodeIndex(triangleId, TRIANGLE_ID);
        vector<vector<SimplexId>> localTriangleStars;
        buildTriangleList(nid, nullptr, nullptr, &localTriangleStars);

        SimplexId localTriangleId = triangleId-triangleIntervals_[nid-1]-1;
        if(localStarId >= (SimplexId) localTriangleStars[localTriangleId].size())
          return -2;
        
        return localTriangleStars[localTriangleId][localStarId];
      }
        
      SimplexId getTriangleStarNumber(const SimplexId &triangleId) const{

         #ifndef TTK_ENABLE_KAMIKAZE
          if((triangleId < 0)||(triangleId >= triangleIntervals_.back()))
            return -1;
        #endif

        SimplexId nid = findNodeIndex(triangleId, TRIANGLE_ID);
        vector<vector<SimplexId>> localTriangleStars;
        buildTriangleList(nid, nullptr, nullptr, &localTriangleStars);
        
        return localTriangleStars[triangleId-triangleIntervals_[nid-1]-1].size();
      }
      
      const vector<vector<SimplexId> > *getTriangleStars(){
        if(!hasPreprocessedTriangleStars_){
          return nullptr;
        }
        triangleStarList_.reserve(triangleIntervals_.back()+1);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          vector<vector<SimplexId>> localTriangleStars;
          buildTriangleList(nid, nullptr, nullptr, &localTriangleStars);
          triangleStarList_.insert(triangleStarList_.end(), localTriangleStars.begin(), localTriangleStars.end());
        }
        return &triangleStarList_;
      }
      
      int getTriangleVertex(const SimplexId &triangleId,
        const int &localVertexId, SimplexId &vertexId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if((triangleId < 0)||(triangleId > triangleIntervals_.back()))
            return -1;
          if((localVertexId < 0) ||(localVertexId > 2))
            return -2;
        #endif

        SimplexId nid = findNodeIndex(triangleId, TRIANGLE_ID);
        vector<vector<SimplexId>> localInternalTriangles;
        buildTriangleList(nid, &localInternalTriangles, nullptr, nullptr);

        return localInternalTriangles[triangleId-triangleIntervals_[nid-1]-1][localVertexId];
      }
      
      int getVertexEdge(const SimplexId &vertexId, 
        const int &localEdgeId, SimplexId &edgeId) const{
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if((vertexId < 0)||(vertexId >= vertexNumber_))
            return -1;
          if(localEdgeId < 0)
            return -2;
        #endif
        
        SimplexId nid = findNodeIndex(vertexId, VERTEX_ID);
        SimplexId localVertexId = vertexId - vertexIntervals_[nid-1] - 1;
        vector<vector<SimplexId>> localVertexEdges;
        if(!getVertexEdges(nid, localVertexEdges)){
          if(localEdgeId >= (SimplexId) localVertexEdges[localVertexId].size())
            return -2;
          edgeId = localVertexEdges[localVertexId][localEdgeId];
          return 0;
        }

        return -3;
      }
      
      SimplexId getVertexEdgeNumber(const SimplexId &vertexId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if((vertexId < 0)||(vertexId >= vertexNumber_))
            return -1;
        #endif

        SimplexId nid = findNodeIndex(vertexId, VERTEX_ID);
        SimplexId localVertexId = vertexId - vertexIntervals_[nid-1] - 1;
        vector<vector<SimplexId>> localVertexEdges;
        if(!getVertexEdges(nid, localVertexEdges)){
          return localVertexEdges[localVertexId].size();
        }
        return -2;
      }
      
      const vector<vector<SimplexId> > *getVertexEdges(){
        if(!hasPreprocessedVertexEdges_)
          return nullptr;

        vertexEdgeList_.reserve(vertexNumber_);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          vector<vector<SimplexId>> localVertexEdges;
          getVertexEdges(nid, localVertexEdges);
          vertexEdgeList_.insert(vertexEdgeList_.end(), localVertexEdges.begin(), localVertexEdges.end());
        }

        return &vertexEdgeList_;
      }
      
      int getVertexLink(const SimplexId &vertexId, 
        const int &localLinkId, SimplexId &linkId) const{
        return 0;
      }
      
      SimplexId getVertexLinkNumber(const SimplexId &vertexId) const{
        return 0;
      }
      
      const vector<vector<SimplexId> > *getVertexLinks(){
        vector<vector<SimplexId>> *dummyPointer= new vector<vector<SimplexId>>();
        return dummyPointer;
      }
      
      int getVertexNeighbor(const SimplexId &vertexId, 
        const int &localNeighborId, SimplexId &neighborId) const{
        #ifndef TTK_ENABLE_KAMIKAZE
          if((vertexId < 0)||(vertexId >= vertexNumber_))
            return -1;
          if(localNeighborId < 0)
            return -2;
        #endif
        
        SimplexId nid = findNodeIndex(vertexId, VERTEX_ID);
        SimplexId localVertexId = vertexId - vertexIntervals_[nid-1] - 1;
        vector<vector<SimplexId>> localVertexNeighbors;
        if(!getVertexNeighbors(nid, localVertexNeighbors)){
          if(localNeighborId < (SimplexId) localVertexNeighbors[localVertexId].size())
            return -2;  
          neighborId = localVertexNeighbors[localVertexId][localNeighborId];
          return 0;
        }
        return -3;
      }
      
      SimplexId getVertexNeighborNumber(const SimplexId &vertexId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if((vertexId < 0)||(vertexId >= vertexNumber_))
            return -1;
        #endif

        SimplexId nid = findNodeIndex(vertexId, VERTEX_ID);
        SimplexId localVertexId = vertexId - vertexIntervals_[nid-1] - 1;
        vector<vector<SimplexId>> localVertexNeighbors;
        if(!getVertexNeighbors(nid, localVertexNeighbors)){
          return localVertexNeighbors[localVertexId].size();
        }
        return -2;
      }
      
      const vector<vector<SimplexId> > *getVertexNeighbors(){
        vertexNeighborList_.reserve(vertexNumber_);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          vector<vector<SimplexId>> localVertexNeighbors;
          getVertexNeighbors(nid, localVertexNeighbors);
          vertexNeighborList_.insert(vertexNeighborList_.end(), localVertexNeighbors.begin(), localVertexNeighbors.end());
        }
        return &vertexNeighborList_;
      }
      
      int getVertexPoint(const SimplexId &vertexId,
        float &x, float &y, float &z) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if((vertexId < 0)||(vertexId >= vertexNumber_))
            return -1;
        #endif
        
        if(doublePrecision_){
          x = ((double *) pointSet_)[3*vertexId];
          y = ((double *) pointSet_)[3*vertexId + 1];
          z = ((double *) pointSet_)[3*vertexId + 2];
        }
        else{
          x = ((float *) pointSet_)[3*vertexId];
          y = ((float *) pointSet_)[3*vertexId + 1];
          z = ((float *) pointSet_)[3*vertexId + 2];
        }
      
        return 0;
      }
      
      int getVertexStar(const SimplexId &vertexId, const int &localStarId,
        SimplexId &starId) const{
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if((vertexId < 0)||(vertexId >= vertexNumber_))
            return -1;
          if(localStarId < 0)
            return -2;
        #endif

        SimplexId nid = findNodeIndex(vertexId, VERTEX_ID);
        vector<vector<SimplexId>> localVertexStars;
        if(!getVertexStars(nid, localVertexStars)){
          SimplexId localVertexId = vertexId - vertexIntervals_[nid-1] - 1;
          if(localStarId >= (SimplexId) localVertexStars[localVertexId].size())
            return -2;
          starId = localVertexStars[localVertexId][localStarId];
          return 0;
        }
        return -3;
      }
      
      SimplexId getVertexStarNumber(const SimplexId &vertexId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if((vertexId < 0)||(vertexId >= vertexNumber_))
            return -1;
        #endif

        SimplexId nid = findNodeIndex(vertexId, VERTEX_ID);
        vector<vector<SimplexId>> localVertexStars;
        if(!getVertexStars(nid, localVertexStars)){
          return localVertexStars[vertexId-vertexIntervals_[nid-1]-1].size();
        }
        return -2;
      }
      
      const vector<vector<SimplexId> > *getVertexStars(){
        vertexStarList_.reserve(vertexNumber_);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          vector<vector<SimplexId>> localVertexStars;
          getVertexStars(nid, localVertexStars);
          for(SimplexId i = 0; i < (SimplexId) localVertexStars.size(); i++){
            vertexStarList_.insert(vertexStarList_.end(), localVertexStars.begin(), localVertexStars.end());
          }
        }
        return &vertexStarList_;
      }
      
      int getVertexTriangle(const SimplexId &vertexId, 
        const int &localTriangleId, SimplexId &triangleId) const{
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if((vertexId < 0)||(vertexId >= (SimplexId) vertexTriangleList_.size()))
            return -1;
          if(localTriangleId < 0)
            return -2;
        #endif

        SimplexId nid = findNodeIndex(vertexId, VERTEX_ID);
        vector<vector<SimplexId>> localVertexTriangles;
        SimplexId localVertexId = vertexId - vertexIntervals_[nid-1] - 1;
        if(!getVertexTriangles(nid, localVertexTriangles)){
          if(localTriangleId >= (SimplexId) localVertexTriangles[localVertexId].size())
            return -2;
          triangleId = localVertexTriangles[localVertexId][localTriangleId];
          return 0;
        }
        return -3;
      }
      
      SimplexId getVertexTriangleNumber(const SimplexId &vertexId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if((vertexId < 0)||(vertexId >= vertexNumber_))
            return -1;
        #endif

        SimplexId nid = findNodeIndex(vertexId, VERTEX_ID);
        vector<vector<SimplexId>> localVertexTriangles;
        if(!getVertexTriangles(nid, localVertexTriangles)){
          return localVertexTriangles[vertexId-vertexIntervals_[nid-1]-1].size();
        }
        return -2;
      }
      
      const vector<vector<SimplexId> > *getVertexTriangles(){
        vertexTriangleList_.reserve(vertexNumber_);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          vector<vector<SimplexId>> localVertexTriangles;
          getVertexTriangles(nid, localVertexTriangles);
          for(SimplexId i = 0; i < (SimplexId) localVertexTriangles.size(); i++){
            vertexTriangleList_.insert(vertexTriangleList_.end(), localVertexTriangles.begin(), localVertexTriangles.end());
          }
        }
        return &vertexTriangleList_;
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

        #ifndef TTK_ENABLE_KAMIKAZE
          if(vertexNumber_ <= 0)
            return -1;
          if(cellNumber_ <= 0)
            return -2;
          if(!cellArray_)
            return -3;
          if(nodeNumber_ <= 0)
            return -4;
        #endif

        if(!edgeIntervals_.size()){
          edgeIntervals_.resize(nodeNumber_+1);
          edgeIntervals_[0] = -1;
          for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
            vector<pair<SimplexId, SimplexId>> localInternalEdges, localExternalEdges;
            SimplexId edgeNum = buildEdgeList(nid, localInternalEdges, localExternalEdges);
            edgeIntervals_[nid] = edgeIntervals_[nid-1] + edgeNum;
          }
        }
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

        #ifndef TTK_ENABLE_KAMIKAZE
          if(vertexNumber_ <= 0)
            return -1;
          if(cellNumber_ <= 0)
            return -2;
          if(!cellArray_)
            return -3;
          if(nodeNumber_ <= 0)
            return -4;
        #endif

        // build triangle interval list
        if(!triangleIntervals_.size()){
          triangleIntervals_.resize(nodeNumber_+1);
          triangleIntervals_[0] = -1;
          for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
            vector<vector<SimplexId>> localInternalTriangles, localExternalTriangles;
            SimplexId triangleNum = buildTriangleList(nid, &localInternalTriangles, &localExternalTriangles, nullptr);
            triangleIntervals_[nid] = triangleIntervals_[nid-1] + triangleNum;
          }
        }

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
        preprocessEdges();
        hasPreprocessedVertexNeighbors_ = true;
        return 0;
      }
      
      int preprocessVertexStars(){

        #ifndef TTK_ENABLE_KAMIKAZE
          if(!cellArray_)
            return -1;
        #endif
        
        vertexStarList_.resize(vertexNumber_);
        for(SimplexId i = 0; i < cellNumber_; i++){
          for(SimplexId j = 0; j < cellArray_[0];j++){
            vertexStarList_[cellArray_[(cellArray_[0] + 1)*i + 1 + j]].push_back(i);
          }
        }

        hasPreprocessedVertexStars_ = true;
        return 0;
      }
      
      int preprocessVertexTriangles(){
        preprocessTriangles();
        hasPreprocessedVertexTriangles_ = true;
        return 0;
      }

      
    protected:

      int clear();

      /**
       * Find the corresponding node index given the id.
       */ 
      SimplexId findNodeIndex(SimplexId id, int idType) const{
          const vector<SimplexId> *intervals = NULL;
          // determine which vector to search
          if(idType == VERTEX_ID){
            intervals = &vertexIntervals_;
          }else if(idType == EDGE_ID){
            intervals = &edgeIntervals_;
          }else if(idType == TRIANGLE_ID){
            intervals = &triangleIntervals_;
          }else if(idType == CELL_ID){
            intervals = &cellIntervals_;
          }else{
            return -1;
          }

          // use binary search to find the first element that is greater than or equal to the given id
          SimplexId low = 0, high = intervals->size()-1;
          while(low <= high){
              SimplexId mid = low + (high-low)/2;
              if(intervals->at(mid) == id){
                  return mid;
              }else if(intervals->at(mid) < id){
                  low = mid + 1;
              }else{
                  high = mid - 1;
              }
          }
          return low;
      }

      /** 
       * Build the whole edge list and return the number of edges in the node.
       */ 
      SimplexId buildEdgeList(SimplexId nodeId, 
        vector<pair<SimplexId, SimplexId>> &internalEdgeList,
        vector<pair<SimplexId, SimplexId>> &externalEdgeList) const{

        SimplexId edgeCount = 0, verticesPerCell = cellArray_[0];
        // loop through all the internal cells first
        internalEdgeList = vector<pair<SimplexId, SimplexId>>(vertexIntervals_[nodeId]-vertexIntervals_[nodeId-1]);
        externalEdgeList = vector<pair<SimplexId, SimplexId>>(vertexIntervals_[nodeId]-vertexIntervals_[nodeId-1]);
        vector<vector<SimplexId>> edgeTable(vertexIntervals_[nodeId]-vertexIntervals_[nodeId-1]);

        // loop through the internal cell list
        for(SimplexId cid = cellIntervals_[nodeId-1]+1; cid <= cellIntervals_[nodeId]; cid++){
          pair<SimplexId, SimplexId> edgeIds;
          SimplexId cellId = (verticesPerCell + 1) * cid;

          // loop through each edge of the cell
          for(SimplexId j = 0; j < verticesPerCell-1; j++){
            edgeIds.first = cellArray_[cellId + j + 1];
            // the edge does not belong to the current node
            if(edgeIds.first > vertexIntervals_[nodeId]){
              break;
            }
            for(SimplexId k = j+1; k < verticesPerCell; k++){
              edgeIds.second = cellArray_[cellId + k + 1];
              
              // search the edge table to see if the edge has already existed
              bool hasFound = false;
              SimplexId localVertexId = edgeIds.first-vertexIntervals_[nodeId-1]-1;
              for(SimplexId l = 0; l < (SimplexId) edgeTable[localVertexId].size(); l++){
                if(edgeIds.second == edgeTable[localVertexId][l]){
                  hasFound = true;
                  break;
                }
              }

              // not found in the edge table - assign new edge id
              if(!hasFound){
                edgeTable[localVertexId].push_back(edgeIds.second);
                internalEdgeList.push_back(edgeIds);
                edgeCount++;
              }
            }
          }
        }

        // loop through the external cell list
        for(SimplexId cid : externalCells_[nodeId]){
          pair<SimplexId, SimplexId> edgeIds;
          SimplexId cellId = (verticesPerCell + 1) * cid;

          // loop through each edge of the cell
          for(SimplexId j = 0; j < verticesPerCell-1; j++){
            for(SimplexId k = j+1; k < verticesPerCell; k++){
              edgeIds.first = cellArray_[cellId + j + 1];
              edgeIds.second = cellArray_[cellId + k + 1];
              
              // the edge is in the current node
              if(edgeIds.first > vertexIntervals_[nodeId-1] && edgeIds.first <= vertexIntervals_[nodeId]){
                bool hasFound = false;
                SimplexId localVertexId = edgeIds.first-vertexIntervals_[nodeId-1]-1;
                for(SimplexId l = 0; l < (SimplexId) edgeTable[localVertexId].size(); l++){
                  if(edgeIds.second == edgeTable[localVertexId][l]){
                    hasFound = true;
                    break;
                  }
                }
                if(!hasFound){
                  edgeTable[localVertexId].push_back(edgeIds.second);
                  internalEdgeList.push_back(edgeIds);
                  edgeCount++;
                }
              }
              // the edge is a bridge between two nodes
              else if(edgeIds.second > vertexIntervals_[nodeId-1] && edgeIds.second <= vertexIntervals_[nodeId]){
                bool hasFound = false;
                SimplexId localVertexId = edgeIds.second-vertexIntervals_[nodeId-1]-1;
                for(SimplexId l = 0; l < (SimplexId) edgeTable[localVertexId].size(); l++){
                  if(edgeIds.first == edgeTable[localVertexId][l]){
                    hasFound = true;
                    break;
                  }
                }
                // not found in the edge table
                if(!hasFound){
                  // no need to increase the edge count cuz it has been counted in the previous node
                  edgeTable[localVertexId].push_back(edgeIds.first);
                  externalEdgeList.push_back(edgeIds);
                }
              }
            }
          }
        }

        return edgeCount;
      }

      /**
       * Build the whole triangle list and return the number of triangles in the node.
       */
      int buildTriangleList(SimplexId nodeId, 
        vector<vector<SimplexId>> *internalTriangleList,
        vector<vector<SimplexId>> *externalTriangleList, 
        vector<vector<SimplexId>> *triangleStars) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 && nodeId > nodeNumber_)
            return -1;
        #endif

        SimplexId triangleCount = 0, verticesPerCell = cellArray_[0];
        SimplexId localVertexNum = vertexIntervals_[nodeId]-vertexIntervals_[nodeId-1];

        if(internalTriangleList){
          // NOTE: 9 is pretty empirical here...
          internalTriangleList->clear();
          internalTriangleList->reserve(9*localVertexNum);
        }
        if(externalTriangleList){
          externalTriangleList->clear();
          externalTriangleList->reserve(9*localVertexNum);
        }
        if(triangleStars){
          triangleStars->clear();
          triangleStars->reserve(9*localVertexNum);
        }
        
        vector<vector<pair<vector<SimplexId>, SimplexId>>> triangleTable(localVertexNum);

        // loop through the internal cell list
        for(SimplexId cid = cellIntervals_[nodeId-1]+1; cid <= cellIntervals_[nodeId]; cid++){
          vector<SimplexId> triangleIds(3);
          SimplexId cellId = (verticesPerCell + 1) * cid;

          // loop through each triangle of the cell
          for(SimplexId j = 0; j < verticesPerCell-2; j++){
            triangleIds[0] = cellArray_[cellId + j + 1];
            // the triangle does not belong to the current node
            if(triangleIds[0] > vertexIntervals_[nodeId]){
              break;
            }
            for(SimplexId k = j+1; k < verticesPerCell-1; k++){
              for(SimplexId l = k+1; l < verticesPerCell-2; l++){
                triangleIds[1] = cellArray_[cellId + k + 1];
                triangleIds[2] = cellArray_[cellId + l + 1];
                
                // search the triangle table to see if the triangle has already existed
                SimplexId triangleId = -1;
                SimplexId localVertexId = triangleIds[0]-vertexIntervals_[nodeId-1]-1;
                for(SimplexId t = 0; t < (SimplexId) triangleTable[localVertexId].size(); t++){
                  if(triangleIds == triangleTable[localVertexId][t].first){
                    triangleId = triangleTable[localVertexId][t].second;
                    break;
                  }
                }

                // not found in the triangle table - assign new triangle id
                if(triangleId == -1){
                  triangleTable[localVertexId].push_back(
                    pair<vector<SimplexId>, SimplexId>(triangleIds, triangleCount));
                  triangleCount++;
                  if(internalTriangleList){
                    internalTriangleList->push_back(triangleIds);
                  }
                  if(triangleStars){
                    triangleStars->push_back(vector<SimplexId>(cid));
                  }
                }
                else{
                  if(triangleStars){
                    (*triangleStars)[triangleId].push_back(cid);
                  }
                }
              }
            }
          }
        }

        // loop through the external cell list
        for(SimplexId cid : externalCells_[nodeId]){
          vector<SimplexId> triangleIds(3);
          SimplexId cellId = (verticesPerCell + 1) * cid;

          // loop through each triangle of the cell
          for(SimplexId j = 0; j < verticesPerCell-2; j++){
            for(SimplexId k = j+1; k < verticesPerCell-1; k++){
              for(SimplexId l = k+1; l < verticesPerCell; l++){
                triangleIds[0] = cellArray_[cellId + j + 1];
                triangleIds[1] = cellArray_[cellId + k + 1];
                triangleIds[2] = cellArray_[cellId + l + 1];

                for(SimplexId v = 0; v < 3; v++){
                  if(triangleIds[v] > vertexIntervals_[nodeId-1] && triangleIds[v] <= vertexIntervals_[nodeId]){
                    SimplexId triangleId = -1;
                    SimplexId localVertexId = triangleIds[v]-vertexIntervals_[nodeId-1]-1;
                    for(SimplexId t = 0; t < (SimplexId) triangleTable[localVertexId].size(); t++){
                      if(triangleIds == triangleTable[localVertexId][t].first){
                        triangleId = triangleTable[localVertexId][t].second;
                        break;
                      }
                    }
                    if(triangleId == -1){
                      triangleTable[localVertexId].push_back(
                        pair<vector<SimplexId>, SimplexId>(triangleIds, triangleCount));
                      if(v == 0){
                        triangleCount++;
                        if(internalTriangleList){
                          internalTriangleList->push_back(triangleIds);
                        }
                        if(triangleStars){
                          triangleStars->push_back(vector<SimplexId>(cid));
                        }
                      }
                      else if(externalTriangleList){
                        externalTriangleList->push_back(triangleIds);
                      }
                    }
                    else{
                      if(v == 0 && triangleStars){
                        (*triangleStars)[triangleId].push_back(cid);
                      }
                    }
                    break;
                  }
                }
              }
            }
          }
        }

        return triangleCount;
      }

      /**
       * Get the cell edges for all cells in a given node. 
       */
      int getCellEdges(const SimplexId &nodeId, vector<vector<SimplexId>> &cellEdges) const{
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 && nodeId > nodeNumber_)
            return -1;
        #endif

        cellEdges = vector<vector<SimplexId>>(cellIntervals_[nodeId]-cellIntervals_[nodeId-1]);
        for(SimplexId i = 0; i < (SimplexId) cellEdges.size(); i++){
          cellEdges[i].reserve(6);
        }
        vector<vector<SimplexId>> vertexEdges;
        getVertexEdges(nodeId, vertexEdges);
        for(SimplexId i = 0; i < (SimplexId) cellEdges.size(); i++){
          SimplexId cellId = (cellArray_[0]+1)*(i+cellIntervals_[nodeId-1]+1);
          for(SimplexId j = 0; j < cellArray_[0]-1; j++){
            for(SimplexId k = j+1; k < cellArray_[0]; k++){
              pair<SimplexId, SimplexId> edgePair(cellArray_[cellId+1+j], cellArray_[cellId+1+k]);
              SimplexId edgeId = getEdgeId(edgePair);
              cellEdges[i].push_back(edgeId);
            }
          }
        }
        return 0;
      }

      /**
       * Get the edge id given a pair of vertex ids.
       * Note: make sure the passed pair is already sorted.
       */
      SimplexId getEdgeId(pair<SimplexId, SimplexId> &edgePair) const{
        SimplexId nid = findNodeIndex(edgePair.first, VERTEX_ID);
        vector<pair<SimplexId, SimplexId>> localInternalEdges, localExternalEdges;
        buildEdgeList(nid, localInternalEdges, localExternalEdges);

        for(SimplexId i = 0; i < (SimplexId) localInternalEdges.size(); i++){
          if(localInternalEdges[i] == edgePair){
            return (i + edgeIntervals_[nid-1] + 1);
          }
        }

        return -1;
      }

      /**
       * Get the edge id given a pair of vertex ids.
       * Note: make sure the passed pair is already sorted.
       */
      SimplexId getTriangleId(vector<SimplexId> &triangle) const{
        SimplexId nid = findNodeIndex(triangle.front(), VERTEX_ID);
        vector<vector<SimplexId>> localInternalTriangles;
        buildTriangleList(nid, &localInternalTriangles, nullptr, nullptr);

        for(SimplexId i = 0; i < (SimplexId) localInternalTriangles.size(); i++){
          if(localInternalTriangles[i] == triangle){
            return (i + triangleIntervals_[nid-1] + 1);
          }
        }

        return -1;
      }
      
      /**
       * Get the vertex edges for all the vertices in a given node.
       */
      int getVertexEdges(const SimplexId &nodeId, vector<vector<SimplexId>> &vertexEdges) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 && nodeId > nodeNumber_)
            return -1;
        #endif

        vertexEdges = vector<vector<SimplexId>>(vertexIntervals_[nodeId]-vertexIntervals_[nodeId-1]);
        vector<pair<SimplexId, SimplexId>> localInternalEdges, localExternalEdges;
        buildEdgeList(nodeId, localInternalEdges, localExternalEdges);

        for(SimplexId i = 0; i < (SimplexId) localInternalEdges.size(); i++){
          vertexEdges[localInternalEdges[i].first-vertexIntervals_[nodeId-1]-1].push_back(edgeIntervals_[nodeId]+i+1);
          if(localInternalEdges[i].second <= vertexIntervals_[nodeId])
            vertexEdges[localInternalEdges[i].second-vertexIntervals_[nodeId-1]-1].push_back(edgeIntervals_[nodeId]+i+1);
        }

        for(SimplexId i = 0; i < (SimplexId) localExternalEdges.size(); i++){
          vertexEdges[localExternalEdges[i].second-vertexIntervals_[nodeId-1]-1].push_back(getEdgeId(localExternalEdges[i]));
        }

        return 0;
      }

      /**
       * Get the vertex neighbors for all the vertices in a given node.
       */
      int getVertexNeighbors(const SimplexId &nodeId, vector<vector<SimplexId>> &vertexNeighbors) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 && nodeId > nodeNumber_)
            return -1;
        #endif

        vertexNeighbors = vector<vector<SimplexId>>(vertexIntervals_[nodeId]-vertexIntervals_[nodeId-1]);
        vector<pair<SimplexId, SimplexId>> localInternalEdges, localExternalEdges;
        buildEdgeList(nodeId, localInternalEdges, localExternalEdges);

        for(SimplexId i = 0; i < (SimplexId) localInternalEdges.size(); i++){
          vertexNeighbors[localInternalEdges[i].first-vertexIntervals_[nodeId-1]-1].push_back(localInternalEdges[i].second);
          if(localInternalEdges[i].second <= vertexIntervals_[nodeId])
            vertexNeighbors[localInternalEdges[i].second-vertexIntervals_[nodeId-1]-1].push_back(i+edgeIntervals_[nodeId]+1);
        }

        for(SimplexId i = 0; i < (SimplexId) localExternalEdges.size(); i++){
          vertexNeighbors[localExternalEdges[i].second-vertexIntervals_[nodeId-1]-1].push_back(localExternalEdges[i].first);
        }

        return 0;
      }

      /** 
       * Get the vertex stars for all the vertices in a given node. 
       * The function is kind of like getVertexCells().
       */ 
      int getVertexStars(const SimplexId &nodeId, vector<vector<SimplexId>> &vertexStars) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 && nodeId > nodeNumber_)
            return -1;
        #endif

        vertexStars = vector<vector<SimplexId>>(vertexIntervals_[nodeId]-vertexIntervals_[nodeId-1]);
        // loop through the internal cell list
        for(SimplexId cid = cellIntervals_[nodeId-1]+1; cid <= cellIntervals_[nodeId]; cid++){
          SimplexId cellId = (cellArray_[0]+1) * cid;
          for(SimplexId j = 0; j < cellArray_[0]; j++){
            // see if it is in the current node
            if(cellArray_[cellId+j+1] > vertexIntervals_[nodeId-1] && cellArray_[cellId+j+1] <= vertexIntervals_[nodeId])
              vertexStars[cellArray_[cellId+j+1]-vertexIntervals_[nodeId-1]-1].push_back(cid);
          }
        }

        // and also external cell list
        for(SimplexId cid : externalCells_[nodeId]){
          SimplexId cellId = (cellArray_[0]+1) * cid;
          for(SimplexId j = 0; j < cellArray_[0]; j++){
            // see if it is in the current node
            if(cellArray_[cellId+j+1] > vertexIntervals_[nodeId-1] && cellArray_[cellId+j+1] <= vertexIntervals_[nodeId])
              vertexStars[cellArray_[cellId+j+1]-vertexIntervals_[nodeId-1]-1].push_back(cid);
          }
        }

        return 0;
      }

      /** 
       * Get the vertex triangles for all the vertices in a given node. 
       */ 
      int getVertexTriangles(const SimplexId &nodeId, vector<vector<SimplexId>> &vertexTriangles) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 && nodeId > nodeNumber_)
            return -1;
        #endif

        vector<vector<SimplexId>> localInternalTriangles, localExternalTriangles;
        buildTriangleList(nodeId, &localInternalTriangles, &localExternalTriangles, nullptr);

        vector<vector<SimplexId>> localVertexTriangles(vertexIntervals_[nodeId]-vertexIntervals_[nodeId-1]);
        SimplexId triangleId = triangleIntervals_[nodeId-1]+1;
        for(SimplexId i = 0; i < (SimplexId) localInternalTriangles.size(); i++){
          for(SimplexId j = 0; j < 3; j++){
            if(localInternalTriangles[i][j] > vertexIntervals_[nodeId-1] && localInternalTriangles[i][j] <= vertexIntervals_[nodeId])
              vertexTriangles[localInternalTriangles[i][j]-vertexIntervals_[nodeId-1]-1].push_back(triangleId);
          }
          triangleId++;
        }
        for(SimplexId i = 0; i < (SimplexId) localExternalTriangles.size(); i++){
          triangleId = getTriangleId(localExternalTriangles[i]);
          for(SimplexId j = 0; j < 3; j++){
            if(localExternalTriangles[i][j] > vertexIntervals_[nodeId-1] && localExternalTriangles[i][j] <= vertexIntervals_[nodeId])
              vertexTriangles[localExternalTriangles[i][j]-vertexIntervals_[nodeId-1]-1].push_back(triangleId);
          }
        }
        
        return 0;
      }

      
      /**
       * Protected class variables.
       */ 
      bool                doublePrecision_;
      SimplexId           cellNumber_, vertexNumber_, nodeNumber_;
      const void          *pointSet_;
      vector<SimplexId>   vertexIntervals_;
      vector<SimplexId>   edgeIntervals_;
      vector<SimplexId>   triangleIntervals_;
      vector<SimplexId>   cellIntervals_;
      const LongSimplexId *cellArray_;
      vector<vector<SimplexId>> externalCells_;
  };
}