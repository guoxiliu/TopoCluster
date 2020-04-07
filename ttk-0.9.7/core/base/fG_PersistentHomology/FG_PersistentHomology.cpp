#include                  <FG_PersistentHomology.h>

using namespace std;
using namespace ttk;

FG_PersistentHomology::FG_PersistentHomology(){

    matrix =       vector<vector<int> >();
    allpairs =     map<int,int>();
    homology =     vector<int>();
}

FG_PersistentHomology::~FG_PersistentHomology(){

}

void FG_PersistentHomology::getBoundarySimplices(Simplex simpl, vector<Simplex>& boundary){

    boundary.clear();

    int dim = simpl.first;
    queue<SimplexId> visitCell;
    visitCell.push(simpl.second);

    set<Simplex> bd_set;

    while(!visitCell.empty()){

      SimplexId simpl = visitCell.front();
      visitCell.pop();

      for(int i=0; i<dim+1; i++){
          SimplexId face = extractBoundary(Simplex(dim,simpl),dim-1,i);

          Simplex pair;
          bool isPaired = getPair(Simplex(dim-1,face), pair);

          if(isPaired && pair.first == dim && pair.second != simpl ){

              visitCell.push(pair.second);
          }
          else if(!isPaired){
            if(bd_set.find(Simplex(dim-1,face)) == bd_set.end())
              bd_set.insert(Simplex(dim-1,face));
            else
              bd_set.erase(Simplex(dim-1,face));
          }
      }
  }

  boundary = vector<Simplex>(bd_set.begin(), bd_set.end());
}

void FG_PersistentHomology::readHomology(vector<float>& points, vector<char>& dims){

  vector<float> coords;
  for(int i=0; i<homology.size(); i++){
      Simplex simpl = index_to_simplex[homology[i]];
      computeBarycenter(simpl, coords);
      points.push_back(coords[0]);
      points.push_back(coords[1]);
      points.push_back(coords[2]);
      dims.push_back(simpl.first);

  }
}


void FG_PersistentHomology::collectSeeds(int dim, SimplexId face, vector<SimplexId>& cofaces){

  queue<SimplexId> visitCell;
  visitCell.push(face);

  set<int> visited;
  visited.insert(face);

  while(!visitCell.empty()){

      SimplexId simplex = visitCell.front();
      visitCell.pop();

      for(int i=0; i<sizeCofacets(Simplex(dim,simplex), dim+1); i++){
          SimplexId coface = extractCoboundary(Simplex(dim,simplex),dim+1,i);

          Simplex pair;
          Simplex head = Simplex(dim+1,coface);
          bool isPaired = getPair(head, pair);

          if(isPaired && pair.first == dim && pair.second != simplex && visited.find(pair.second) == visited.end()){
              visitCell.push(pair.second);
              visited.insert(pair.second);
          }
          else if(!isPaired){
              cofaces.push_back(coface);
          }
      }
  }
}

void FG_PersistentHomology::readGradientVectors(vector<float> &points, vector<float> &vects) {

    for(int i=0; i<triangulation_->getNumberOfVertices(); i++){
        Simplex simpl(0,i);
        Simplex next;
        bool isPaired = getPair(simpl, next);

        if(isPaired){
            vector<float> coords;
            computeBarycenter(simpl, coords);
            points.push_back(coords[0]);
            points.push_back(coords[1]);
            points.push_back(coords[2]);

            vector<float> coords2;
            computeBarycenter(next, coords2);
            vects.push_back(coords2[0] - coords[0]);
            vects.push_back(coords2[1] - coords[1]);
            vects.push_back(coords2[2] - coords[2]);
        }
    }

    for(int i=0; i<triangulation_->getNumberOfEdges(); i++){
        Simplex simpl(1,i);
        Simplex next;
        bool isPaired = getPair(simpl, next);

        if(isPaired && next.first == 2){
            vector<float> coords;
            computeBarycenter(simpl, coords);
            points.push_back(coords[0]);
            points.push_back(coords[1]);
            points.push_back(coords[2]);

            vector<float> coords2;
            computeBarycenter(next, coords2);
            vects.push_back(coords2[0] - coords[0]);
            vects.push_back(coords2[1] - coords[1]);
            vects.push_back(coords2[2] - coords[2]);
        }
    }

    if(triangulation_->getDimensionality() == 3){
        for(int i=0; i<triangulation_->getNumberOfTriangles(); i++){
            Simplex simpl(2,i);
            Simplex next;
            bool isPaired = getPair(simpl, next);

            if(isPaired && next.first == 3){
                vector<float> coords;
                computeBarycenter(simpl, coords);
                points.push_back(coords[0]);
                points.push_back(coords[1]);
                points.push_back(coords[2]);

                vector<float> coords2;
                computeBarycenter(next, coords2);
                vects.push_back(coords2[0] - coords[0]);
                vects.push_back(coords2[1] - coords[1]);
                vects.push_back(coords2[2] - coords[2]);
            }
        }
    }
}
