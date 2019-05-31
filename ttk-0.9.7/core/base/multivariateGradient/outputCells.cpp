#include                  <MultivariateGradient.h>

using namespace std;
using namespace ttk;

void MultivariateGradient::readOutputGradientPairs(std::vector<float>& pairedPoints){

  int tot = triangulation_->getNumberOfTriangles() + triangulation_->getNumberOfEdges() + triangulation_->getNumberOfVertices();
  pairedPoints = std::vector<float>((tot-numberCriticalSimplices)*3, 0);

  int stepPairs=0;

  for(int i=0;i<=dimensionality_;i++){
    for(int j=0; j<(*gradient)[i].size(); j++){
      Simplex simplex = Simplex(i,j);

      if(!isCritical(simplex)){
        Simplex other;
        getPair(simplex, other);

        if(other.first == simplex.first+1){

          std::vector<float> barycenter;
          computeBarycenter(simplex, barycenter);

          pairedPoints[stepPairs++] = barycenter[0];
          pairedPoints[stepPairs++] = barycenter[1];
          pairedPoints[stepPairs++] = barycenter[2];

          computeBarycenter(other, barycenter);
          pairedPoints[stepPairs++] = barycenter[0];
          pairedPoints[stepPairs++] = barycenter[1];
          pairedPoints[stepPairs++] = barycenter[2];

        }
      }
    }
  }
}

void MultivariateGradient::readOutputCriticalPoints(int index,
                                                          std::vector<float>& points,
                                                          std::vector<char>& dimensions,
                                                          std::vector<char>& type,
                                                          std::vector<int>& labels){


  if(index != -1){

    int cluster = cclusters->getSortedClusterIndex(index);

    if(cclusters->isDead(cluster) || cluster == -1)
        return;


    int size = cclusters->getClusterSize(cluster);

    points = std::vector<float>(size*3, 0);
    dimensions = std::vector<char>(size, 0);
    type = std::vector<char>(size, 0);
    labels = std::vector<int>(size, 0);

    int step=0;

    vector<Simplex> simplices;
    cclusters->getSimplices(cluster, simplices);

    for(auto simplex : simplices){
      std::vector<float> barycenter;
      computeBarycenter(simplex, barycenter);

      points[3*step] = barycenter[0];
      points[3*step+1] = barycenter[1];
      points[3*step+2] = barycenter[2];

      dimensions[step] = simplex.first;
      labels[step] = cluster;
      type[step++] = cclusters->getType(cluster);
    }

  }
  else{

    points = std::vector<float>(numberCriticalSimplices*3, 0);
    dimensions = std::vector<char>(numberCriticalSimplices, 0);
    type = std::vector<char>(numberCriticalSimplices, 0);
    labels = std::vector<int>(numberCriticalSimplices, 0);

    int step=0;
    for(int i=0; i<cclusters->size(); i++){
      if(cclusters->isDead(i))
        continue;

      vector<Simplex> simplices;
      cclusters->getSimplices(i, simplices);

      for(auto simplex : simplices){
        std::vector<float> barycenter;
        computeBarycenter(simplex, barycenter);

        points[3*step] = barycenter[0];
        points[3*step+1] = barycenter[1];
        points[3*step+2] = barycenter[2];

        dimensions[step] = simplex.first;
        labels[step] = i;
        type[step++] = cclusters->getType(i);
      }
    }
  }
}


void MultivariateGradient::getOutputClusterGraph(list<float>& points, list<pair<int,int> >& cells,
                                                        std::vector<int>& csize,
                                                        std::vector<char>& type,
                                                        std::vector<int>& labels){

    unordered_map<int,int> reindex;
    int iter=0;
    for(int i=0; i<cclusters->size(); i++){
        if(cclusters->isDead(i))
          continue;

        reindex[i]=iter++;

        Simplex simpl;
        cclusters->getClusterBestDraw(i,simpl);
        vector<float> barycenter;
        computeBarycenter(simpl,barycenter);

        for(auto val : barycenter){
          points.push_back(val);
        }

        type.push_back(cclusters->getType(i));
        csize.push_back(cclusters->getClusterSize(i));
        labels.push_back(i);

        unordered_map<int,int> arcs;
        cclusters->getConnections(i,arcs,true);
        for(auto val : arcs){
          if(val.first < i){
            cells.push_back(pair<int,int>(reindex[i],reindex[val.first]));
          }
        }

        cclusters->getConnections(i,arcs,false);
        for(auto val : arcs){
          if(val.first < i){
            cells.push_back(pair<int,int>(reindex[i],reindex[val.first]));
          }

        }
    }

}

void MultivariateGradient::fullSetOfClusters(int cluster, unordered_set<int>& all_clusters){

    queue<int> clusters;
    clusters.push(cluster);

    all_clusters.insert(cluster);

    while(!clusters.empty()){
      int lcl = clusters.front();
      clusters.pop();

      unordered_set<int> local_clusters=unordered_set<int>();
      cclusters->getMergedClusters(lcl, local_clusters);

      for(auto ret_cl : local_clusters){
        if(all_clusters.find(ret_cl) == all_clusters.end()){
          clusters.push(ret_cl);
          all_clusters.insert(ret_cl);
        }
      }
    }

    cout << "Done" << endl;
}

void MultivariateGradient::readOutputDescendingCell(int clIndex, list<float>& point_coords, list<vector<int> >& cells, vector<int>& nCells){

  int mainCluster = cclusters->getSortedClusterIndex(clIndex);
  cout << "New Index Cluster for Descending Region " << mainCluster << " " << cclusters->size() << endl;

  if(cclusters->isDead(mainCluster) || mainCluster == -1)
      return;

  unordered_set<int> found_edges;
  unordered_set<int> found_triangles;
  unordered_set<int> found_tetrahedra;

  unordered_set<int> all_clusters=unordered_set<int>();
  fullSetOfClusters(mainCluster,all_clusters);
  all_clusters.insert(mainCluster);

  cout << "Clusters collected " << endl;

  for(auto cluster : all_clusters){
    vector<Simplex> cSimplices;
    cclusters->getSimplices(cluster, cSimplices);

    //for each critical simplex we extract the descending geometry
    for(auto s : cSimplices){
        if(s.first > 0){
            CMap critical_reached;
            if(s.first == 1){
                if(cclusters->isDead(cluster))
                  found_edges.insert(s.second);
                extractOutgoingPath(s, critical_reached, true, found_edges);
            }
            else if(s.first == 2){
                if(cclusters->isDead(cluster))
                  found_triangles.insert(s.second);
                extractOutgoingPath(s, critical_reached, true, found_triangles);
            }
            else{
                if(cclusters->isDead(cluster))
                  found_tetrahedra.insert(s.second);
                extractOutgoingPath(s, critical_reached, true, found_tetrahedra);
            }
        }
    }
  }

  cout << "Region collected " << endl;

  //create set of points
  unordered_map<int,int> points;
  int nVertices=0;
  int nTetra=0;
  int nTriangles=0;
  int nEdges=0;

  unordered_set<int>* ptrSimplSet;
  for(int i=1; i<4; i++){

    if(i == 1) ptrSimplSet = &found_edges;
    else if(i == 2) ptrSimplSet = &found_triangles;
    else ptrSimplSet = &found_tetrahedra;

    for(unordered_set<int>::iterator it=ptrSimplSet->begin(); it != ptrSimplSet->end(); it++){
      vector<SimplexId> vertices;
      simplexToVertices(Simplex(i,*it), vertices);

      for(auto v : vertices){
          if(points.find(v) == points.end()){
              points[v]=nVertices++;
          }
      }

      //if true we need this simplex because it is not completely on the boundary of someone else
      vector<SimplexId> newExplicitCell(vertices.size(),0);
      for(int vi = 0; vi<vertices.size(); vi++){
          newExplicitCell[vi] = points[vertices[vi]];
      }

      cells.push_back(newExplicitCell);

      nCells[i] = nCells[i] + 1;
    }
  }

  vector<int> true_points(points.size());
  for(auto i : points){
    true_points[i.second] = i.first;
  }

  //prepare final number of simplices
  nCells[0] = nEdges;
  nCells[1] = nTriangles;
  nCells[2] = nTetra;

  //prepare point coordinates
  for(auto i : true_points){
      vector<float> vCoords(3);
      triangulation_->getVertexPoint(i, vCoords[0], vCoords[1], vCoords[2]);
      point_coords.push_back(vCoords[0]);
      point_coords.push_back(vCoords[1]);
      point_coords.push_back(vCoords[2]);
  }

  cout << "I should have " << points.size() << " points and " << nEdges << " " << nTriangles << " " << nTetra << endl;
}

void MultivariateGradient::readOutputAscendingCell(int clIndex, list<float>& point_coords, list<vector<int> >& cells, vector<int>& nCells){

  int mainCluster = cclusters->getSortedClusterIndex(clIndex);
  cout << "New Index Cluster for Ascending Region " << mainCluster << endl;

  if(cclusters->isDead(mainCluster) || mainCluster == -1)
      return;

  unordered_set<int> all_clusters=unordered_set<int>();
  fullSetOfClusters(mainCluster,all_clusters);
  all_clusters.insert(mainCluster);

  unordered_set<int> found_edges;
  unordered_set<int> found_triangles;
  unordered_set<int> found_tetrahedra;

  for(auto cluster : all_clusters){
    vector<Simplex> cSimplices;
    cclusters->getSimplices(cluster, cSimplices);

    //for each critical simplex we extract the descending geometry
    for(auto s : cSimplices){
      if(s.first < dimensionality_){
          unordered_set<int> critical_reached;
          if(s.first == 0)
              extractIngoingPath(s, critical_reached, true, found_edges);
          else if(s.first == 1){
            if(cclusters->isDead(cluster))
              found_edges.insert(s.second);

            extractIngoingPath(s, critical_reached, true, found_triangles);
          }
          else{
            if(cclusters->isDead(cluster))
              found_triangles.insert(s.second);
            extractIngoingPath(s, critical_reached, true, found_tetrahedra);
          }
      }
      else{
        if(cclusters->isDead(cluster))
          found_tetrahedra.insert(s.second);
      }
    }
  }


  //create set of points
  unordered_map<int,int> points;
  int nVertices=0;
  int nTetra=0;
  int nTriangles=0;
  int nEdges=0;

  unordered_set<int>* ptrSimplSet;
  for(int i=1; i<4; i++){

    if(i == 1) ptrSimplSet = &found_edges;
    else if(i == 2) ptrSimplSet = &found_triangles;
    else ptrSimplSet = &found_tetrahedra;

    for(unordered_set<int>::iterator it=ptrSimplSet->begin(); it != ptrSimplSet->end(); it++){
      vector<SimplexId> vertices;
      simplexToVertices(Simplex(i,*it), vertices);

      for(auto v : vertices){
          if(points.find(v) == points.end()){
              points[v]=nVertices++;
          }
      }

      //if true we need this simplex because it is not completely on the boundary of someone else
      vector<SimplexId> newExplicitCell(vertices.size(),0);
      for(int vi = 0; vi<vertices.size(); vi++){
          newExplicitCell[vi] = points[vertices[vi]];
      }

      cells.push_back(newExplicitCell);
      nCells[i] = nCells[i] + 1;
    }
  }

  vector<int> true_points(points.size());
  for(auto i : points){
    true_points[i.second] = i.first;
  }

  //prepare final number of simplices
  nCells[0] = nEdges;
  nCells[1] = nTriangles;
  nCells[2] = nTetra;

  //prepare point coordinates
  for(auto i : true_points){
      vector<float> vCoords(3);
      triangulation_->getVertexPoint(i, vCoords[0], vCoords[1], vCoords[2]);
      point_coords.push_back(vCoords[0]);
      point_coords.push_back(vCoords[1]);
      point_coords.push_back(vCoords[2]);
  }

}
