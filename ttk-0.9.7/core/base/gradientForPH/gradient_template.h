//#ifndef MultivariateGradient_TEMPLATE_H
#define GRADIENT_TEMPLATE_H

#include "GradientForPH.h"

using namespace std;
using namespace ttk;


template <typename dataType>
int GradientForPH::computeGradient(){

    Timer t;

    // check the consistency of the variables -- to adapt
    #ifndef TTK_ENABLE_KAMIKAZE
    if(!triangulation_)
    return -1;
    if(!inputData_)
    return -2;
    #endif

    //prepare gradient
    gradient = new vector<vector<vector<SimplexId> > >();

    (*gradient).push_back(vector<vector<SimplexId> >(triangulation_->getNumberOfVertices(), vector<SimplexId>(1,-1)));
    (*gradient).push_back(vector<vector<SimplexId> >(triangulation_->getNumberOfEdges(), vector<SimplexId>(2,-1)));

    if(dimensionality_ == 2){
        (*gradient).push_back(vector<vector<SimplexId> >(triangulation_->getNumberOfTriangles(), vector<SimplexId>(1,-1)));
    }
    else{
        (*gradient).push_back(vector<vector<SimplexId> >(triangulation_->getNumberOfTriangles(), vector<SimplexId>(2,-1)));
        (*gradient).push_back(vector<vector<SimplexId> >(triangulation_->getNumberOfCells(), vector<SimplexId>(1,-1)));
    }

    //setup indexing
    computeIndexing<dataType>();

    cout << "Indexing computed" << endl;

    //compute indexing on the data
    SimplexId vertexNumber = triangulation_->getNumberOfVertices();

    #ifdef TTK_ENABLE_OPENMP
    #pragma omp parallel for 
    #endif
    for(SimplexId i=0; i<vertexNumber; i++){
        //   cout << "Working on vertex " << i << endl;
        vector<vector<SimplexId> > lower_star;

        //extract lower star for each vertex based on first component
        computeLowerStar(i, lower_star);
        //   cout << "Computed lower star" << endl;

        //compute homotopy expansion
        int local = assignGradient<dataType>(lower_star);
        //   cout << "Paired gradients" << endl;
        //

        #ifdef TTK_ENABLE_OPENMP
        #pragma omp critical
        #endif
        numberCriticalSimplices += local;
    }

    cout << "Gradient computed" << endl;


//#ifdef TTK_ENABLE_OPENMP
//#pragma omp parallel for num_threads(threadNumber_)
//#endif
//  for(SimplexId i = 0; i < vertexNumber; i++){
//    cout <<  triangulation_->getVertexStarNumber(i) << " " << triangulation_->getVertexTriangleNumber(i) << " - ";
//    for(SimplexId j = 0; j < triangulation_->getVertexEdgeNumber(i); j++){
//        SimplexId edge;
//        triangulation_->getVertexEdge(i,j,edge);
//        cout << edge << " ";
//    }
//    cout << endl;
//  }

  return 0;
}

template <typename dataType>
void GradientForPH::computeIndexing() {

    dataType* inputData = (dataType*) inputData_;
    SimplexId vertexNumber = triangulation_->getNumberOfVertices();

    vertexIndexing = new std::vector<SimplexId>(vertexNumber,-1);

    map<dataType, vector<int> > sorted;
    for(int i=0; i<vertexNumber; i++){
        dataType val = inputData[i];
        if(sorted.find(val) == sorted.end()){
            sorted[val]=vector<int>();
        }
        sorted[val].push_back(i);
    } 
    int index=0;
    for(auto levelset : sorted){
        for(auto v : levelset.second){
            (*vertexIndexing)[v]=index++;
        }
    }

}

template <typename dataType>
dataType ttk::GradientForPH::getFieldValue(Simplex simplex){

    dataType* inputData = (dataType*) inputData_;

    //from the simplex extract its boundary vertices
    vector<SimplexId> vertices;
    simplexToVertices(simplex,vertices);

    //from the vertices select the field value of the simplex
    dataType fieldValue = inputData[vertices[0]];

    for(SimplexId v=1; v<vertices.size(); v++){
            fieldValue = fieldValue > inputData[vertices[v]] ? fieldValue : inputData[vertices[v]];
    }

    return fieldValue;
}


template <typename dataType>
int GradientForPH::assignGradient(SimplexesVec& vec){

    auto lexycographic = bind(&GradientForPH::cmpSimplexesLevelSets, this, _1, _2);
    vector<SimplexesSet > sset(2, SimplexesSet(lexycographic));

    //count the number of simplices in the lower star
    int toClassify = 0;
    int local_critical=0;
    for(int d=0; d<vec.size(); d++){
        toClassify += vec[d].size();
    }

    //populate the set of vertices
    for(auto v : vec[0]){
        sset[1].insert(Simplex(0,v));
    }

    //set the working dimension. It means we will start searching for edges to be paired with vertices
    short int workingDim = 1;

    //continue until you have classified all the simplices as either paired or critical
    while(toClassify != 0){

        if(workingDim <= dimensionality_){

            //populate d and d-1 simplices sets (d == workingDim)
            sset[0] = sset[1];
            sset[1].clear();

            for(auto v : vec[workingDim])
                sset[1].insert(Simplex(workingDim,v));


            //start searching for pairings
            while(sset[0].size() != 0){
                //here pair simplices

                list<Simplex> toRemove = list<Simplex>();
                for(Simplex s : sset[1]){


                    Simplex chosen;
                    int numPairable = numPairableFaces<dataType>(s, sset[0], chosen);

                    if(numPairable == 1){

                        setPair(chosen,s);
                        sset[0].erase(chosen);
                        toRemove.push_back(s);
                        toClassify -= 2;
                    }

                }

                if(!toRemove.empty()){
                    for(auto s : toRemove){
                        sset[1].erase(s);
                    }
                }
                else{
                    Simplex critical=*sset[0].begin();
                    sset[0].erase(critical);
                    toClassify--;
                    local_critical++;
                    
                }

            }

            workingDim++;

        }
        else{
            //here we should classify the remaining simplices as critical
            Simplex critical=*sset[0].begin();
            sset[0].erase(critical);
            toClassify--;
            local_critical++;
        }

    }

    return local_critical;
}

template <typename dataType>
int GradientForPH::numPairableFaces(Simplex simplex, SimplexesSet& pairable, Simplex &chosen){

    dataType field = getFieldValue<dataType>(simplex);
    int numBoundary = simplex.first+1;
    uint num_pairable=0;

    SimplexId id;
    for(int i=0; i<numBoundary; i++){
        switch(simplex.first){
        case 1:
            triangulation_->getEdgeVertex(simplex.second, i, id);
            break;

        case 2:
            triangulation_->getTriangleEdge(simplex.second, i, id);
            break;

        case 3:
            triangulation_->getCellTriangle(simplex.second, i, id);
            break;
        }

        Simplex candidate = Simplex(simplex.first-1, id);

        //see if contained
        if(pairable.find(candidate) != pairable.end())
        {

            //see if same field value
            dataType field2 = getFieldValue<dataType>(candidate);

            if(field == field2){
                chosen = candidate;
                num_pairable++;
            }
        }
    }

    return num_pairable;
}


//#endif // MultivariateGradient_TEMPLATE_H
