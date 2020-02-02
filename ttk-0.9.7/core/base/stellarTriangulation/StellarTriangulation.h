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
#include                  <set>
#include                  <map>

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

        while(nodeNum <= nodeNumber_){
          SimplexId startPos = (cellArray[0]+1)*cid+1;
          cell = vector<SimplexId>(cellArray+startPos, cellArray+startPos+cellArray[0]);

          if(cell[0] > vertexIntervals_[nodeNum] || cid >= cellNumber){
            cellIntervals_[nodeNum++] = cid - 1;
            continue;
          }
          // create external cell list 
          for(size_t i = 1; i < cell.size(); i++){
            if(cell[i] > vertexIntervals_[nodeNum]){
            SimplexId nid = findNodeIndex(cell[i], VERTEX_ID);
              if(externalCells_[nid].empty()){
                externalCells_[nid].push_back(cid);
              }else if(externalCells_[nid].back() != cid){
                externalCells_[nid].push_back(cid);
              }
            }
          }
          if(cid < cellNumber)
            cid++;
        }
        
        // TEST: print out the cell intervals
        cout << "[StellarTriangulation] Cell intervals: \n";
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
        cout << "[StellarTriangulation] Node number: " << nodeNumber_ << ".\n Vertex intervals:\n";
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
        cellEdgeList_.reserve(cellNumber_);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          vector<vector<SimplexId>> localCellEdges;
          getCellEdges(nid, localCellEdges);
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

        #ifndef TTK_ENABLE_KAMIKAZE
        if((cellId < 0)||(cellId >= cellIntervals_.back()))
          return -1;
        if((localTriangleId < 0)||(localTriangleId >= getCellTriangleNumber(cellId)))
          return -2;
        #endif

        SimplexId nid = findNodeIndex(cellId, CELL_ID);
        vector<vector<SimplexId>> localCellTriangles;
        if(!getCellTriangles(nid, localCellTriangles)){
          triangleId = localCellTriangles[cellId-cellIntervals_[nid-1]-1][localTriangleId];
          return 0;
        }
        return -3;
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
        cellTriangleList_.reserve(edgeIntervals_.back()+1);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          vector<vector<SimplexId>> localCellTriangles;
          getCellTriangles(nid, localCellTriangles);
          cellTriangleList_.insert(cellTriangleList_.end(), localCellTriangles.begin(), localCellTriangles.end());
        }
        return &cellTriangleList_;
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
          vector<pair<SimplexId, SimplexId>> localInternalEdges;
          buildEdgeList(nid, &localInternalEdges, nullptr, nullptr, nullptr);
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
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if((edgeId < 0)||(edgeId > edgeIntervals_.back()))
            return -1;
          if(localStarId < 0)
            return -2;
        #endif

        SimplexId nid = findNodeIndex(edgeId, EDGE_ID);
        vector<vector<SimplexId>> localEdgeStars;
        buildEdgeList(nid, nullptr, nullptr, nullptr, &localEdgeStars);
        SimplexId localEdgeId = edgeId - edgeIntervals_[nid-1] - 1;
        if(localStarId >= (SimplexId) localEdgeStars[localEdgeId].size())
          return -2;
        starId = localEdgeStars[localEdgeId][localStarId];
        return 0;
      }
        
      SimplexId getEdgeStarNumber(const SimplexId &edgeId) const{
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if((edgeId < 0)||(edgeId > edgeIntervals_.back()))
            return -1;
        #endif

        SimplexId nid = findNodeIndex(edgeId, EDGE_ID);
        vector<vector<SimplexId>> localEdgeStars;
        buildEdgeList(nid, nullptr, nullptr, nullptr, &localEdgeStars);

        SimplexId localEdgeId = edgeId - edgeIntervals_[nid-1] - 1;
        return localEdgeStars[localEdgeId].size();
      }
      
      const vector<vector<SimplexId> > *getEdgeStars(){
        edgeStarList_.reserve(edgeIntervals_.back()+1);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          vector<vector<SimplexId>> localEdgeStars;
          buildEdgeList(nid, nullptr, nullptr, nullptr, &localEdgeStars);
          edgeStarList_.insert(edgeStarList_.end(), localEdgeStars.begin(), localEdgeStars.end());
        }
        return &edgeStarList_;
      }
     
      int getEdgeTriangle(const SimplexId &edgeId, 
        const int &localTriangleId, SimplexId &triangleId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
        if((edgeId < 0)||(edgeId > (SimplexId) edgeIntervals_.size()))
          return -1;
        if(localTriangleId < 0)
          return -2;
        #endif

        SimplexId nodeId = findNodeIndex(edgeId, EDGE_ID);
        vector<vector<SimplexId> > localEdgeTriangles;
        getEdgeTriangles(nodeId, localEdgeTriangles);
        SimplexId localEdgeId = edgeId-edgeIntervals_[nodeId-1]-1;
        if(localTriangleId >= (SimplexId) localEdgeTriangles[localEdgeId].size())
          return -2;
        triangleId = localEdgeTriangles[localEdgeId][localTriangleId];

        return 0;
      }
        
      SimplexId getEdgeTriangleNumber(const SimplexId &edgeId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
        if((edgeId < 0)||(edgeId > (SimplexId) edgeIntervals_.size()))
          return -1;
        #endif
        
        if(getDimensionality() == 2){
          return getEdgeStarNumber(edgeId);
        }
        else{
          SimplexId nodeId = findNodeIndex(edgeId,EDGE_ID);
          vector<vector<SimplexId> > localEdgeTriangles;
          getEdgeTriangles(nodeId, localEdgeTriangles);
          return localEdgeTriangles[edgeId-edgeIntervals_[nodeId-1]-1].size();
        }

        return -2;
      }
      
      const vector<vector<SimplexId> > *getEdgeTriangles(){
        edgeTriangleList_.reserve(edgeIntervals_.size()+1);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          vector<vector<SimplexId>> localEdgeTriangles;
          getEdgeTriangles(nid, localEdgeTriangles);
          edgeTriangleList_.insert(edgeTriangleList_.end(), localEdgeTriangles.begin(), localEdgeTriangles.end());
        }
        return &edgeTriangleList_;
      }
      
      int getEdgeVertex(const SimplexId &edgeId, 
        const int &localVertexId, SimplexId &vertexId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if((edgeId < 0)||(edgeId > (SimplexId) edgeIntervals_.back()))
            return -1;
          if((localVertexId != 0)&&(localVertexId != 1))
            return -2;
        #endif

        SimplexId nid = findNodeIndex(edgeId, EDGE_ID);
        vector<pair<SimplexId, SimplexId>> localInternalEdges;
        buildEdgeList(nid, &localInternalEdges, nullptr, nullptr, nullptr);
        SimplexId localEdgeId = edgeId-edgeIntervals_[nid-1]-1;
        if(localVertexId){
          vertexId = localInternalEdges[localEdgeId].second;
        }
        else{
          vertexId = localInternalEdges[localEdgeId].first;
        }
        return 0;
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
                    
        #ifndef TTK_ENABLE_KAMIKAZE
          if((triangleId < 0)||(triangleId > triangleIntervals_.back()))
            return -1;
          if((localEdgeId < 0)||(localEdgeId > 2))
            return -2;
        #endif

        SimplexId nid = findNodeIndex(triangleId, TRIANGLE_ID);
        vector<vector<SimplexId>> localTriangleEdges;
        if(!getTriangleEdges(nid, localTriangleEdges)){
          edgeId = localTriangleEdges[triangleId-triangleIntervals_[nid-1]-1][localEdgeId];
          return 0;
        }
        return -3;
      }
      
      SimplexId getTriangleEdgeNumber(const SimplexId &triangleId) const{
        #ifndef TTK_ENABLE_KAMIKAZE
          if((triangleId < 0)||(triangleId > triangleIntervals_.back()))
            return -1;
        #endif

        SimplexId nid = findNodeIndex(triangleId, TRIANGLE_ID);
        vector<vector<SimplexId>> localTriangleEdges;
        if(!getTriangleEdges(nid, localTriangleEdges)){
          return localTriangleEdges[triangleId-triangleIntervals_[nid-1]-1].size();
        }
        return -2;
      }
      
      const vector<vector<SimplexId> > *getTriangleEdges(){
        if(!hasPreprocessedTriangleEdges_){
          return nullptr;
        }
        triangleEdgeList_.reserve(triangleIntervals_.size()+1);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          vector<vector<SimplexId>> localTriangleEdges;
          getTriangleEdges(nid, localTriangleEdges);
          triangleEdgeList_.insert(triangleEdgeList_.end(), localTriangleEdges.begin(), localTriangleEdges.end());
        }
        return &triangleEdgeList_;
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
          if((triangleId < 0)||(triangleId > triangleIntervals_.back()))
            return -1;
          if(localStarId < 0)
            return -2;
        #endif

        SimplexId nid = findNodeIndex(triangleId, TRIANGLE_ID);
        vector<vector<SimplexId>> localTriangleStars;
        buildTriangleList(nid, nullptr, nullptr, nullptr, &localTriangleStars);

        SimplexId localTriangleId = triangleId-triangleIntervals_[nid-1]-1;
        if(localStarId >= (SimplexId) localTriangleStars[localTriangleId].size())
          return -2;
        
        return localTriangleStars[localTriangleId][localStarId];
      }
        
      SimplexId getTriangleStarNumber(const SimplexId &triangleId) const{

         #ifndef TTK_ENABLE_KAMIKAZE
          if((triangleId < 0)||(triangleId > triangleIntervals_.back()))
            return -1;
        #endif

        SimplexId nid = findNodeIndex(triangleId, TRIANGLE_ID);
        vector<vector<SimplexId>> localTriangleStars;
        buildTriangleList(nid, nullptr, nullptr, nullptr, &localTriangleStars);
        
        return localTriangleStars[triangleId-triangleIntervals_[nid-1]-1].size();
      }
      
      const vector<vector<SimplexId> > *getTriangleStars(){
        if(!hasPreprocessedTriangleStars_){
          return nullptr;
        }
        triangleStarList_.reserve(triangleIntervals_.back()+1);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          vector<vector<SimplexId>> localTriangleStars;
          buildTriangleList(nid, nullptr, nullptr, nullptr, &localTriangleStars);
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
        buildTriangleList(nid, &localInternalTriangles, nullptr, nullptr, nullptr);
        vertexId = localInternalTriangles[triangleId-triangleIntervals_[nid-1]-1][localVertexId];
        return 0;
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
          if((vertexId < 0)||(vertexId >= vertexNumber_))
            return -1;
          if(localTriangleId < 0)
            return -2;
        #endif

        SimplexId nid = findNodeIndex(vertexId, VERTEX_ID);
        vector<vector<SimplexId>> localVertexTriangles;
        SimplexId localVertexId = vertexId - vertexIntervals_[nid-1] - 1;
        if(!getVertexTriangles(nid, localVertexTriangles)){
          if(localTriangleId >= (SimplexId) localVertexTriangles[localVertexId].size()){
            return -2;
          }
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

        if(!hasPreprocessedEdges_){
          edgeIntervals_.resize(nodeNumber_+1);
          edgeIntervals_[0] = -1;
          for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
            // map<pair<SimplexId,SimplexId>,SimplexId> localInternalEdgeMap;
            // vector<pair<SimplexId, SimplexId>> localInternalEdges, localExternalEdges;
            // SimplexId edgeNum = buildEdgeList(nid, &localInternalEdges, &localExternalEdges, &localInternalEdgeMap);
            vector<pair<SimplexId, SimplexId>> localInternalEdges;
            SimplexId edgeNum = buildEdgeList(nid, &localInternalEdges, nullptr, nullptr, nullptr);
            edgeIntervals_[nid] = edgeIntervals_[nid-1] + edgeNum;
          }
          hasPreprocessedEdges_ = true;
        }

        cout << "Edge intervals: " << endl;
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          cout << nid << ": " << edgeIntervals_[nid] << endl;
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
          if(nodeNumber_ <= 0)
            return -4;
        #endif

        // build triangle interval list
        if(!hasPreprocessedTriangles_){
          triangleIntervals_.resize(nodeNumber_+1);
          triangleIntervals_[0] = -1;
          for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
            // vector<vector<SimplexId>> localInternalTriangles, localExternalTriangles;
            // map<vector<SimplexId>, SimplexId> localInternalTriangleMap;
            // SimplexId triangleNum = buildTriangleList(nid, &localInternalTriangles, &localExternalTriangles, &localInternalTriangleMap, nullptr);
            vector<vector<SimplexId>> localInternalTriangles;
            SimplexId triangleNum = buildTriangleList(nid, &localInternalTriangles, nullptr, nullptr, nullptr);
            triangleIntervals_[nid] = triangleIntervals_[nid-1] + triangleNum;
          }
          hasPreprocessedTriangles_ = true;
        }

        cout << "Triangle intervals: " << endl;
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          cout << nid << ": " << triangleIntervals_[nid] << endl;
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
        const vector<SimplexId> *intervals = nullptr;
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

        vector<SimplexId>::const_iterator low = lower_bound(intervals->begin(), intervals->end(), id);
        return (low-intervals->begin());
      }

      /** 
       * Build the whole edge list and return the number of edges in the node.
       */ 
      SimplexId buildEdgeList(SimplexId nodeId, 
        vector<pair<SimplexId, SimplexId>> * const internalEdgeList,
        vector<pair<SimplexId, SimplexId>> * const externalEdgeList,
        map<pair<SimplexId, SimplexId>, SimplexId> * const internalEdgeMap,
        vector<vector<SimplexId>> * const edgeStars) const{
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 || nodeId > nodeNumber_)
            return -1;
        #endif

        SimplexId edgeCount = 0, verticesPerCell = cellArray_[0];
        map<pair<SimplexId, SimplexId>, SimplexId> *edgeMap = new map<pair<SimplexId, SimplexId>, SimplexId>();
        // loop through all the internal cells first
        if(internalEdgeList){
          internalEdgeList->clear();
          internalEdgeList->reserve(6*(vertexIntervals_[nodeId]-vertexIntervals_[nodeId-1]));
        }
        if(externalEdgeList){
          externalEdgeList->clear();
          externalEdgeList->reserve(6*(vertexIntervals_[nodeId]-vertexIntervals_[nodeId-1]));
        }
        if(edgeStars){
          edgeStars->clear();
          edgeStars->reserve(6*(vertexIntervals_[nodeId]-vertexIntervals_[nodeId-1]));
        }
        if(internalEdgeMap){
          edgeMap = internalEdgeMap;
        }

        set<pair<SimplexId, SimplexId>> externalEdgeSet;

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
              
              // not found in the edge map - assign new edge id
              if(edgeMap->find(edgeIds) == edgeMap->end()){
                edgeCount++;
                (*edgeMap)[edgeIds] = edgeCount+edgeIntervals_[nodeId-1];
                if(internalEdgeList)
                  internalEdgeList->push_back(edgeIds);
                if(edgeStars){
                  edgeStars->push_back(vector<SimplexId>{cid});
                }
              }
              else if(edgeStars){
                (*edgeStars)[edgeMap->at(edgeIds)-edgeIntervals_[nodeId-1]-1].push_back(cid);
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
                if(edgeMap->find(edgeIds) == edgeMap->end()){
                  edgeCount++;
                  (*edgeMap)[edgeIds] = edgeCount+edgeIntervals_[nodeId-1];
                  if(internalEdgeList)
                    internalEdgeList->push_back(edgeIds);
                  if(edgeStars)
                    edgeStars->push_back(vector<SimplexId>{cid});
                }
                else if(edgeStars){
                  (*edgeStars)[edgeMap->at(edgeIds)-edgeIntervals_[nodeId-1]-1].push_back(cid);
                }
              }
              // the edge is a bridge between two nodes
              else if(externalEdgeList && edgeIds.second > vertexIntervals_[nodeId-1] && edgeIds.second <= vertexIntervals_[nodeId]){
                externalEdgeSet.insert(edgeIds);
              }
            }
          }
        }

        if(externalEdgeList)
          externalEdgeList->insert(externalEdgeList->end(), externalEdgeSet.begin(), externalEdgeSet.end());
        
        return edgeCount;
      }

      /**
       * Build the whole triangle list and return the number of triangles in the node.
       */
      int buildTriangleList(SimplexId nodeId, 
        vector<vector<SimplexId>> * const internalTriangleList,
        vector<vector<SimplexId>> * const externalTriangleList, 
        map<vector<SimplexId>,SimplexId> * const internalTriangleMap,
        vector<vector<SimplexId>> * const triangleStars) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 || nodeId > nodeNumber_)
            return -1;
        #endif

        SimplexId triangleCount = 0, verticesPerCell = cellArray_[0];
        SimplexId localVertexNum = vertexIntervals_[nodeId]-vertexIntervals_[nodeId-1];
        map<vector<SimplexId>,SimplexId> *triangleMap = new map<vector<SimplexId>,SimplexId>();
        set<vector<SimplexId>> externalTriangleSet;

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
        if(internalTriangleMap){
          triangleMap = internalTriangleMap;
        }


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
              for(SimplexId l = k+1; l < verticesPerCell; l++){
                triangleIds[1] = cellArray_[cellId + k + 1];
                triangleIds[2] = cellArray_[cellId + l + 1];
                
                if(triangleMap->find(triangleIds) == triangleMap->end()){
                  triangleCount++;
                  (*triangleMap)[triangleIds] = triangleCount + triangleIntervals_[nodeId-1];
                  if(internalTriangleList){
                    internalTriangleList->push_back(triangleIds);
                  }
                  if(triangleStars){
                    triangleStars->push_back(vector<SimplexId>{cid});
                  }
                }
                else if(triangleStars){
                  (*triangleStars)[triangleMap->at(triangleIds)-triangleIntervals_[nodeId-1]-1].push_back(cid);
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
                    if(triangleMap->find(triangleIds) == triangleMap->end()){
                      if(v == 0){
                        triangleCount++;
                        (*triangleMap)[triangleIds] = triangleCount + triangleIntervals_[nodeId-1];
                        if(internalTriangleList){
                          internalTriangleList->push_back(triangleIds);
                        }
                        if(triangleStars){
                          triangleStars->push_back(vector<SimplexId>{cid});
                        }
                      }
                      else if(externalTriangleList){
                        externalTriangleSet.insert(triangleIds);
                      }
                    }
                    else{
                      if(v == 0 && triangleStars){
                        (*triangleStars)[triangleMap->at(triangleIds)-triangleIntervals_[nodeId-1]-1].push_back(cid);
                      }
                    }
                  }
                }
              }
            }
          }
        }

        if(externalTriangleList){
          externalTriangleList->insert(externalTriangleList->end(), externalTriangleSet.begin(), externalTriangleSet.end());
        }

        return triangleCount;
      }

      /**
       * Get the cell edges for all cells in a given node. 
       */
      int getCellEdges(const SimplexId &nodeId, vector<vector<SimplexId>> &cellEdges) const{
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 || nodeId > nodeNumber_)
            return -1;
        #endif

        cellEdges.clear();
        cellEdges.resize(cellIntervals_[nodeId]-cellIntervals_[nodeId-1]);
        for(SimplexId i = 0; i < (SimplexId) cellEdges.size(); i++){
          cellEdges[i].reserve(6);
        }
        map<pair<SimplexId,SimplexId>, SimplexId> internalEdgeMap;
        buildEdgeList(nodeId, nullptr, nullptr, &internalEdgeMap, nullptr);

        for(SimplexId i = 0; i < (SimplexId) cellEdges.size(); i++){
          SimplexId cellId = (cellArray_[0]+1)*(i+1+cellIntervals_[nodeId-1]);
          for(SimplexId j = 0; j < cellArray_[0]-1; j++){
            for(SimplexId k = j+1; k < cellArray_[0]; k++){
              pair<SimplexId, SimplexId> edgePair(cellArray_[cellId+1+j], cellArray_[cellId+1+k]);
              if(edgePair.first > vertexIntervals_[nodeId]){
                cellEdges[i].push_back(getEdgeId(edgePair));
              }
              else{
                cellEdges[i].push_back(internalEdgeMap[edgePair]);
              }
            }
          }
        }
        return 0;
      }

      /**
       * Get the cell triangles for all cells in a given node. 
       */
      int getCellTriangles(const SimplexId &nodeId, vector<vector<SimplexId>> &cellTriangles) const{
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 || nodeId > nodeNumber_)
            return -1;
        #endif

        cellTriangles.clear();
        cellTriangles.resize(cellIntervals_[nodeId]-cellIntervals_[nodeId-1]);
        for(SimplexId i = 0; i < (SimplexId) cellTriangles.size(); i++){
          cellTriangles[i].reserve(4);
        }
        for(SimplexId i = 0; i < (SimplexId) cellTriangles.size(); i++){
          SimplexId cellId = (cellArray_[0]+1)*(i+cellIntervals_[nodeId-1]+1);
          vector<SimplexId> triangleVec(3);
          for(SimplexId j = 0; j < cellArray_[0]-2; j++){
            triangleVec[0] = cellArray_[cellId+1+j]; 
            for(SimplexId k = j+1; k < cellArray_[0]-1; k++){
              triangleVec[1] = cellArray_[cellId+1+k];
              for(SimplexId l = k+1; l < cellArray_[0]; l++){
                triangleVec[2] = cellArray_[cellId+1+l];
                cellTriangles[i].push_back(getTriangleId(triangleVec));
              }
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
        map<pair<SimplexId, SimplexId>, SimplexId> localEdgeMap; 
        buildEdgeList(nid, nullptr, nullptr, &localEdgeMap, nullptr);

        if(localEdgeMap.find(edgePair) == localEdgeMap.end()){
          return -1;
        }

        return localEdgeMap[edgePair];
      }

      /**
       * Get the edge triangles for all the edges in a given node.
       */
      int getEdgeTriangles(const SimplexId& nodeId, vector<vector<SimplexId> >& edgeTriangles) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 || nodeId > nodeNumber_)
              return -1;
        #endif

        edgeTriangles.clear();
        edgeTriangles = vector<vector<SimplexId> >(edgeIntervals_[nodeId] - edgeIntervals_[nodeId-1]);

        vector<vector<SimplexId>> internalTriangles, externalTriangles;
        map<pair<SimplexId, SimplexId>, SimplexId> internalEdgeMap;
        map<vector<SimplexId>, SimplexId> internalTriangleMap;

        buildEdgeList(nodeId, nullptr, nullptr, &internalEdgeMap, nullptr);
        buildTriangleList(nodeId, &internalTriangles, &externalTriangles, &internalTriangleMap, nullptr);

        // for internal triangles
        pair<SimplexId, SimplexId> edgeIds;
        for(SimplexId tid = 0; tid < (SimplexId) internalTriangles.size(); tid++){
          for(SimplexId j = 0; j < 2; j++){
            if(internalTriangles[tid][j] > vertexIntervals_[nodeId]){
              break;
            }
            for(SimplexId k = j+1; k < 3; k++){
              edgeIds = pair<SimplexId,SimplexId>(internalTriangles[tid][j], internalTriangles[tid][k]);
              edgeTriangles[internalEdgeMap[edgeIds]-edgeIntervals_[nodeId-1]-1].push_back(
                tid + triangleIntervals_[nodeId-1] + 1);
            }
          }
        }

        // for external triangles
        for(vector<SimplexId> triangle : externalTriangles){
          // loop through each edge of the cell
          SimplexId triangleId = getTriangleId(triangle);
          for(SimplexId j = 0; j < 2; j++){
            for(SimplexId k = j+1; k < 3; k++){
              edgeIds.first = triangle[j];
              edgeIds.second = triangle[k];
              
              // the edge is in the current node
              if(edgeIds.first > vertexIntervals_[nodeId-1] && edgeIds.first <= vertexIntervals_[nodeId]){
                edgeTriangles[internalEdgeMap[edgeIds]-edgeIntervals_[nodeId-1]-1].push_back(triangleId);
              }
            }
          }
        }

        return 0;
      }

      /**
       * Get the edge id given a pair of vertex ids.
       * Note: make sure the passed pair is already sorted.
       */
      SimplexId getTriangleId(vector<SimplexId> &triangle) const{
        SimplexId nid = findNodeIndex(triangle.front(), VERTEX_ID);
        map<vector<SimplexId>,SimplexId> localTriangleMap;
        buildTriangleList(nid, nullptr, nullptr, &localTriangleMap, nullptr);

        if(localTriangleMap.find(triangle) == localTriangleMap.end())
          return -1;

        return localTriangleMap[triangle];
      }

      /**
       * Get the triangle edges for all the triangles in a given node.
       */
      int getTriangleEdges(const SimplexId &nodeId, vector<vector<SimplexId>> &triangleEdges) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 || nodeId > nodeNumber_)
            return -1;
        #endif

        triangleEdges.clear();
        triangleEdges.resize(triangleIntervals_[nodeId]-triangleIntervals_[nodeId-1], vector<SimplexId>(3));
        vector<vector<SimplexId>> internalTriangles;
        map<pair<SimplexId,SimplexId>, SimplexId> internalEdgeMap;
        buildEdgeList(nodeId, nullptr, nullptr, &internalEdgeMap, nullptr);
        buildTriangleList(nodeId, &internalTriangles, nullptr, nullptr, nullptr);

        for(SimplexId tid = 0; tid < (SimplexId) internalTriangles.size(); tid++){
          // since the first vertex of the triangle is in the node ...
          pair<SimplexId,SimplexId> edgePair(internalTriangles[tid][0], internalTriangles[tid][1]);
          triangleEdges[tid][0] = internalEdgeMap[edgePair];
          edgePair.second = internalTriangles[tid][2];
          triangleEdges[tid][1] = internalEdgeMap[edgePair];
          edgePair.first = internalTriangles[tid][1];
          if(edgePair.first > vertexIntervals_[nodeId-1] && edgePair.second <= vertexIntervals_[nodeId]){
            triangleEdges[tid][2] = internalEdgeMap[edgePair];
          }else{
            triangleEdges[tid][2] = getEdgeId(edgePair);
          }
        }

        return 0;
      }

      /**
       * Get the vertex edges for all the vertices in a given node.
       */
      int getVertexEdges(const SimplexId &nodeId, vector<vector<SimplexId>> &vertexEdges) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 || nodeId > nodeNumber_)
            return -1;
        #endif

        vertexEdges.clear();
        vertexEdges.resize(vertexIntervals_[nodeId]-vertexIntervals_[nodeId-1]);
        vector<pair<SimplexId, SimplexId>> localInternalEdges, localExternalEdges;
        buildEdgeList(nodeId, &localInternalEdges, &localExternalEdges, nullptr, nullptr);

        for(SimplexId i = 0; i < (SimplexId) localInternalEdges.size(); i++){
          vertexEdges[localInternalEdges[i].first-vertexIntervals_[nodeId-1]-1].push_back(edgeIntervals_[nodeId-1]+i+1);
          // the second vertex id of the edge must be greater than the first one
          if(localInternalEdges[i].second <= vertexIntervals_[nodeId])
            vertexEdges[localInternalEdges[i].second-vertexIntervals_[nodeId-1]-1].push_back(edgeIntervals_[nodeId-1]+i+1);
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
          if(nodeId <= 0 || nodeId > nodeNumber_)
            return -1;
        #endif

        vertexNeighbors.clear();
        vertexNeighbors.resize(vertexIntervals_[nodeId]-vertexIntervals_[nodeId-1]);
        vector<pair<SimplexId, SimplexId>> localInternalEdges, localExternalEdges;
        buildEdgeList(nodeId, &localInternalEdges, &localExternalEdges, nullptr, nullptr);

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
       * The function is similar as getVertexCells().
       */ 
      int getVertexStars(const SimplexId &nodeId, vector<vector<SimplexId>> &vertexStars) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 || nodeId > nodeNumber_)
            return -1;
        #endif

        vertexStars.clear();
        vertexStars.resize(vertexIntervals_[nodeId]-vertexIntervals_[nodeId-1]);
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
          if(nodeId <= 0 || nodeId > nodeNumber_)
            return -1;
        #endif

        vertexTriangles.clear();
        vertexTriangles.resize(vertexIntervals_[nodeId] - vertexIntervals_[nodeId-1]);

        vector<vector<SimplexId>> localInternalTriangles, localExternalTriangles;
        buildTriangleList(nodeId, &localInternalTriangles, &localExternalTriangles, nullptr, nullptr);

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
