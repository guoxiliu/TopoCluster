#define PERSISTENTHOMOLOGY_TEMPLATE_H

#include "PersistentHomology.h"


using namespace std;
using namespace ttk;

template <class dataType> int ttk::PersistentHomology::execute(bool onVertices){
 
  // check the consistency of the variables -- to adapt
#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation_)
    return -1;
  if(!inputData_)
    return -2;
#endif

  Timer t;
   
  vector<int> filtration;
  BoundaryMatrix<dataType> bdmatrix;

  if(onVertices){


    //init filtration integers
    filtration = vector<int>(triangulation_->getNumberOfVertices(), 0);
    
    cout << "Preparing filtration" << endl;
    //create discrete filtration for vertices
    t.reStart();
    map<dataType, vector<int> > organize_values;
    for(int i=0; i<triangulation_->getNumberOfVertices(); i++){
      
      Simplex simplex = Simplex(0,i);
      dataType f = getFiltration<dataType>(simplex);

      if(organize_values.find(f) == organize_values.end())
        organize_values[f] = vector<int>();
      organize_values[f].push_back(i);
    }

    int index=0;
    for(auto lvlset : organize_values){
      for(auto v : lvlset.second){
        filtration[v]=index++;
      }
    }

    cout << t.getElapsedTime() << endl;
  t.reStart();
  cout << "Feeding matrix " << endl;
   int global_index=0;   
   for(int i=0; i<triangulation_->getNumberOfVertices(); i++){
        Simplex simplex = Simplex(0,i);
        int f = getFiltration(simplex, filtration);
        bdmatrix.addValue(global_index++,-1,f);
    }

    for(int i=0; i<triangulation_->getNumberOfEdges(); i++){
        Simplex simplex = Simplex(1,i);
        int f = getFiltration(simplex, filtration);
        
        for(int j=0; j<2; j++){
            SimplexId v;
            triangulation_->getEdgeVertex(i,j,v);
            bdmatrix.addValue(global_index,v,f);
        }

        global_index++;
    }

    for(int i=0; i<triangulation_->getNumberOfTriangles(); i++){
        Simplex simplex = Simplex(2,i);
        int f = getFiltration(simplex, filtration);
        
        for(int j=0; j<3; j++){
            SimplexId s;
            triangulation_->getTriangleEdge(i,j,s);
            bdmatrix.addValue(global_index,s+triangulation_->getNumberOfVertices(),f);
        }

        global_index++;
    }

    if(dimensionality_ > 2)
      for(int i=0; i<triangulation_->getNumberOfCells(); i++){
          Simplex simplex = Simplex(3,i);
          int f = getFiltration(simplex, filtration);
          
          for(int j=0; j<4; j++){
              SimplexId s;
              triangulation_->getCellTriangle(i,j,s);
              bdmatrix.addValue(global_index,s+triangulation_->getNumberOfVertices()+triangulation_->getNumberOfEdges(),f);
          }

          global_index++;
      }

    cout << t.getElapsedTime() << endl;
  }
  else{
    //TODO
    cout << "Not yet implemented the case where filtration values are defined on cells" << endl;
  }

  cout << "Sorting" << endl;
  t.reStart();
  bdmatrix.sort();
  cout << t.getElapsedTime() << endl;

  cout << "Reducing" << endl;
  t.reStart();
  bdmatrix.reduce();
  cout << t.getElapsedTime() << endl;

  bdmatrix.getMatrix(matrix);
  bdmatrix.getPairs(allpairs);
  bdmatrix.getHomology(homology);

  return 0;
}

template <class dataType>
dataType ttk::PersistentHomology::getFiltration(Simplex& simplex){

  dataType *inputData = (dataType *) inputData_;

  vector<SimplexId> vertices = vector<SimplexId>();
  simplexToVertices(simplex, vertices);
  
  //from the vertices select the field value of the simplex
  dataType fieldValue = inputData[vertices[0]];
  for(SimplexId v=1; v<vertices.size(); v++){
    fieldValue = fieldValue > inputData[vertices[v]] ? fieldValue : inputData[vertices[v]];
  }

  return fieldValue;
}

template <class dataType>
void PersistentHomology::readPersistencePairs(vector<float>& points, vector<char>& dims, vector<dataType>& filtration){

  vector<float> coords;

  for(auto ppair :allpairs){

    Simplex simpl1,simpl2;
    indexToSimpl(ppair.first, simpl1);
    indexToSimpl(ppair.second,simpl2);

    computeBarycenter(simpl1, coords);
    points.push_back(coords[0]);
    points.push_back(coords[1]);
    points.push_back(coords[2]);
    dims.push_back(simpl1.first);

    computeBarycenter(simpl2, coords);
    points.push_back(coords[0]);
    points.push_back(coords[1]);
    points.push_back(coords[2]);
    dims.push_back(simpl2.first); 

    dataType f1 = getFiltration<dataType>(simpl1);
    dataType f2 = getFiltration<dataType>(simpl2);

    filtration.push_back(abs(f1-f2));

  } 

}




//#endif // PERSISTENTHOMOLOGY_TEMPLATE_H