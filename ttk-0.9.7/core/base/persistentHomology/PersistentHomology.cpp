#include                  <PersistentHomology.h>

#include <unordered_set>

using namespace std;
using namespace ttk;

PersistentHomology::PersistentHomology(){

  inputData_ = NULL;
  triangulation_ = NULL;
}

PersistentHomology::~PersistentHomology(){
  
}




void PersistentHomology::computeBarycenter(Simplex& simplex, vector<float>& coords){

  std::vector<SimplexId> vertices;
  simplexToVertices(simplex, vertices);

  coords = vector<float>(3,0);

  for(auto vertex : vertices){
      vector<float> vCoords(3);
      triangulation_->getVertexPoint(vertex, vCoords[0], vCoords[1], vCoords[2]);

      for(int i=0;i<3;i++)
          coords[i] += vCoords[i];
  }

  for(int i=0;i<3;i++)
      coords[i] = coords[i]/static_cast<float>(vertices.size());

}

int PersistentHomology::getFiltration(Simplex& simplex, const vector<int>& filtration){

  std::vector<SimplexId> vertices;
  simplexToVertices(simplex, vertices);

  //from the vertices select the field value of the simplex
  int fieldValue = 0;

  for(SimplexId v=0; v<vertices.size(); v++){
    fieldValue = fieldValue > filtration[vertices[v]] ? fieldValue : filtration[vertices[v]];
  }

  return fieldValue;
}

void PersistentHomology::indexToSimpl(int i, Simplex& simpl){
      if(i < triangulation_->getNumberOfVertices()){
        simpl = Simplex(0,i);   
      }
      else if(i < triangulation_->getNumberOfVertices() + triangulation_->getNumberOfEdges()){
        simpl = Simplex(1,i-triangulation_->getNumberOfVertices());
      }
      else if(i < triangulation_->getNumberOfVertices() + triangulation_->getNumberOfEdges() + triangulation_->getNumberOfTriangles()){
        simpl = Simplex(2,i-(triangulation_->getNumberOfVertices() + triangulation_->getNumberOfEdges()));
      }
      else{
        simpl = Simplex(3,i-(triangulation_->getNumberOfVertices() + triangulation_->getNumberOfEdges() + triangulation_->getNumberOfTriangles()));
      }
}

void PersistentHomology::simplexToVertices(Simplex& simplex, vector<SimplexId>& vertices){

  vertices = vector<SimplexId>(simplex.first+1);
  for(int i=0; i<simplex.first+1;i++){
    SimplexId vertex;

    switch (simplex.first) {
      case 0:
          vertex = simplex.second;
          break;
      case 1:
          triangulation_->getEdgeVertex(simplex.second,i,vertex);
          break;
      case 2:
          triangulation_->getTriangleVertex(simplex.second,i,vertex);
          break;
      case 3:
          triangulation_->getCellVertex(simplex.second,i,vertex);
          break;
      default:
        cout << "Error, accepted values between 1 and 3 not " << simplex.first << endl;
        break; 
    }

    vertices[i]=vertex;
  }
}

void PersistentHomology::readHomology(vector<float>& points, vector<char>& dims){

  vector<float> coords;
  for(int i=0; i<homology.size(); i++){
      Simplex simpl;
      indexToSimpl(homology[i],simpl);
      computeBarycenter(simpl, coords);
      points.push_back(coords[0]);
      points.push_back(coords[1]);
      points.push_back(coords[2]);
      dims.push_back(simpl.first);
    
  } 
}