#include                  <GradientForPH.h>

using namespace std;
using namespace ttk;

void GradientForPH::readOutputGradientPairs(std::vector<float>& pairedPoints){

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

void GradientForPH::readOutputCriticalPoints(int index, std::vector<float>& points, std::vector<char>& dimensions){

    //TODO IMPLEMENT HERE
}
