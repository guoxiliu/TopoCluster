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

using namespace std;
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

        if((!cellArray_)||(!cellNumber_))
          return -1;
        
        // initialize the array of cell intervals
        cellIntervals_.resize(nodeNumber_);
        externalCells_.resize(nodeNumber_);
        SimplexId nodeNum = 0, cid = 0;
        vector<SimplexId> cell;
        for(; cid < cellNumber; cid++){
          SimplexId startPos = (cellArray[0]+1)*cid+1;
          cell = vector<SimplexId>(cellArray+startPos, cellArray+startPos+cellArray[0]);
          sort(cell.begin(), cell.end());

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

      inline int setInputPoints(const SimplexId &pointNumber, const void *pointSet, 
        const int *indexArray, const bool &doublePrecision = false){

        if(vertexNumber_)
            clear();

        vertexNumber_ = pointNumber;
        pointSet_ = pointSet;
        doublePrecision_ = doublePrecision;

        // initialize the array of vertex intervals
        SimplexId vid = 1;
        for(; vid < pointNumber; vid++){
          if(indexArray[vid] != indexArray[vid-1]){
            vertexIntervals_.push_back(vid-1);
          }
        }
        vertexIntervals_.push_back(vid-1);
        nodeNumber_ = vertexIntervals_.size();

        // TEST: print out the vertex intervals
        cout << "Vertex intervals:\n";
        for(size_t i = 0; i < vertexIntervals_.size(); i++){
          cout << i << ": " << vertexIntervals_[i] << endl;
        }

        return 0;
      }

      int getCellEdge(const SimplexId &cellId, 
        const int &localEdgeId, SimplexId &edgeId) const{
        return 0;
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
        
      const std::vector<std::vector<SimplexId> > *getCellTriangles(){
        std::vector<std::vector<SimplexId>> *dummyPointer= new std::vector<std::vector<SimplexId>>();
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
      
      inline SimplexId getNumberOfCells() const { return cellNumber_; }
      
      inline SimplexId getNumberOfEdges() const{

#ifndef TTK_ENABLE_KAMIKAZE
        if(!edgeIntervals_.size())
          return -1;
#endif
        return edgeIntervals_.back();
      }
      
      SimplexId getNumberOfTriangles() const{
#ifndef TTK_ENABLE_KAMIKAZE
        if(!triangleIntervals_.size())
          return -1;
#endif
        return triangleIntervals_.back();
      }
      
      inline SimplexId getNumberOfVertices() const { return vertexNumber_; }
      
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
        std::vector<std::vector<SimplexId>> *vertexEdges = new std::vector<std::vector<SimplexId>>();

        return vertexEdges;
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

        SimplexId tmpVertexId, verticesPerCell = cellArray_[0];
        vector<vector<SimplexId>> edgeTable(vertexNumber_);
        for(SimplexId i = 0; i < cellNumber_; i++){
          pair<SimplexId, SimplexId> edgeIds;
          SimplexId cellId = (verticesPerCell + 1) * i;
          for(SimplexId j = 0; j < verticesPerCell-1; j++){
            for(SimplexId k = j+1; k < verticesPerCell; k++){
              edgeIds.first = cellArray_[cellId + j + 1];
              edgeIds.second = cellArray_[cellId + k + 1];

              if(edgeIds.first > edgeIds.second){
                tmpVertexId = edgeIds.first;
                edgeIds.first = edgeIds.second;
                edgeIds.second = tmpVertexId;
              }

              bool hasFound = false;
              for(SimplexId l = 0; l < (SimplexId) edgeTable[edgeIds.first].size(); l++){
                if(edgeIds.second == edgeTable[edgeIds.first][l]){
                  hasFound = true;
                  break;
                }
              }
              if(!hasFound){
                edgeTable[edgeIds.first].push_back(edgeIds.second);
              }
            }
          }
        }

        SimplexId edgeCount = 0;
        for(SimplexId i = 0; i < (SimplexId) edgeTable[i].size(); i++){
          edgeCount += edgeTable[i].size();
        }
        edgeList_.resize(edgeCount);

        edgeCount = 0;
        for(SimplexId i = 0; i < (SimplexId) edgeTable.size(); i++){
          for(SimplexId j = 0; j < (SimplexId) edgeTable[i].size(); j++){
            edgeList_[edgeCount].first = i;
            edgeList_[edgeCount].second = edgeTable[i][j];
            edgeCount++;
          }
        }

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
#endif

        if(triangleList_.size()){
          triangleList_.clear();
          triangleList_.reserve(9*vertexNumber_);
        }

        SimplexId triangleNum = 0;
        vector<vector<pair<vector<SimplexId>, SimplexId>>> triangleTable(vertexNumber_);
        for(SimplexId i = 0; i < cellNumber_; i++){
          vector<SimplexId> triangle(3);
          for(int j = 0; j < 4; j++){
            for(int k = 0; k < 3; k++){
              triangle[k] = cellArray_[5*i + 1 + (j + k) % 4];
            }
            sort(triangle.begin(), triangle.end());

            SimplexId triangleId = -1;
            for(SimplexId k = 0; k < (SimplexId) triangleTable[triangle[0]].size(); k++){
              if(triangleTable[triangle[0]][k].first[1] == triangle[1] && triangleTable[triangle[0]][k].first[2] == triangle[2]){
                triangleId = triangleTable[triangle[0]][k].second;
                break;
              }
            }

            if(triangleId == -1){
              triangleId = triangleNum;
              triangleTable[triangle[0]].push_back(make_pair(triangle, triangleNum));
              triangleNum++;
            }

            triangleList_.push_back(triangle);
          }
        }

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

      int clear();

      // Find the corresponding node index given the id.
      int findNodeIndex(int id, int idType) const{
          const vector<SimplexId> *intervals = NULL;
          // determine which vector to search
          if(idType == VERTEX_ID){
              intervals = &vertexIntervals_;
          }else if(idType == EDGE_ID){
              intervals = &edgeIntervals_;
          }else if(idType == TRIANGLE_ID){
              intervals = &cellIntervals_;
          }else{
              return -1;
          }

          // use binary search to find the first element that is greater than or equal to the given id
          int low = 0, high = intervals->size()-1;
          while(low <= high){
              int mid = low + (high-low)/2;
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
      
      bool                doublePrecision_;
      SimplexId           cellNumber_, vertexNumber_, nodeNumber_;
      const void          *pointSet_;
      vector<SimplexId>   vertexIntervals_;
      vector<SimplexId>   edgeIntervals_;
      vector<SimplexId>   triangleIntervals_;
      vector<SimplexId>   cellIntervals_;
      const LongSimplexId *cellArray_;
      vector<vector<SimplexId>> externalEdges_;
      vector<vector<SimplexId>> externalTriangles_;
      vector<vector<SimplexId>> externalCells_;
  };
}