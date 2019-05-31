//#ifndef MultivariateGradient_TEMPLATE_H
#define MULTIVARIATEGRADIENT_TEMPLATE_H

#include "MultivariateGradient.h"

using namespace std;
using namespace ttk;


template <typename dataType>
int ttk::MultivariateGradient::computeGradient(){

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
    #pragma omp parallel for num_threads(threadNumber)
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



  {
    std::stringstream msg;
    msg << "[MultivariateGradient] Data-set (" << vertexNumber
      << " points) processed in "
      << t.getElapsedTime() << " s. (" << threadNumber
      << " thread(s))."
      << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
void ttk::MultivariateGradient::computeIndexing() {

    std::vector<dataType*>* inputData = (std::vector<dataType*>*) inputData_;
    SimplexId vertexNumber = triangulation_->getNumberOfVertices();

    vertexIndexing = new std::vector<SimplexId>(vertexNumber,-1);

    //build up a vector to sort
    vector<vector<dataType> > verticesFieldValues(vertexNumber);
    for(int i=0; i<vertexNumber; i++){
        vector<dataType> values(inputData->size()+1);

        for(int f=0; f<inputData->size(); f++){
            values[f] = inputData->at(f)[i];
        }
        values[inputData->size()] = static_cast<dataType>(i);
        verticesFieldValues[i]=values;
    }

    //create the indexing based on the field values
    sort(verticesFieldValues.begin(), verticesFieldValues.end(),
         [](const std::vector<dataType>& a, const std::vector<dataType>& b) {
            for(int i=0; i<a.size(); i++){

                if(a[i]==b[i])
                    continue;

                return a[i] < b[i];
            }

            return false;
         });

    //save the indexing function
    for(int i=0;i<vertexNumber;i++){
        int index = static_cast<int>(verticesFieldValues[i][inputData->size()]);
        (*vertexIndexing)[index]=i;
    }

}

template <typename dataType>
std::vector<dataType> ttk::MultivariateGradient::getFieldValue(Simplex simplex){

    std::vector<dataType*>* inputData = (std::vector<dataType*>*) inputData_;

    //from the simplex extract its boundary vertices
    vector<SimplexId> vertices;
    simplexToVertices(simplex,vertices);

    //from the vertices select the field value of the simplex
    vector<dataType> fieldValue(inputData->size());

    for(int f=0; f<inputData->size(); f++){
        fieldValue[f] = inputData->at(f)[vertices[0]];
    }

    for(SimplexId v=1; v<vertices.size(); v++){
        for(int f=0; f<inputData->size(); f++){
            fieldValue[f] = fieldValue[f] > inputData->at(f)[vertices[v]] ? fieldValue[f] : inputData->at(f)[vertices[v]];
        }
    }

    return fieldValue;
}


template <typename dataType>
int MultivariateGradient::assignGradient(SimplexesVec& vec){

    auto lexycographic = bind(&MultivariateGradient::cmpSimplexesLevelSets, this, _1, _2);
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

                    // vector<SimplexId> vertices;
                    // simplexToVertices(s,vertices);

                    // for(auto v : vertices){
                    //     vector<dataType> values = getFieldValue<dataType>(Simplex(0,v));
                    //     cout << values[0] << " (" << (*vertexIndexing)[v] << ") ";
                    // }

                    // cout << endl;

                    Simplex chosen;
                    int numPairable = numPairableFaces<dataType>(s, sset[0], chosen);

                    if(numPairable == 1){

                        // vector<SimplexId> vertices;
                        // simplexToVertices(chosen, vertices);
                        // for(auto v : vertices){
                        //     vector<dataType> values = getFieldValue<dataType>(Simplex(0,v));
                        //     cout << values[0] << " (" << (*vertexIndexing)[v] << ") ";
                        // }
                        // cout << " selected" << endl;

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
int MultivariateGradient::numPairableFaces(Simplex simplex, SimplexesSet& pairable, Simplex &chosen){

    vector<dataType> field = getFieldValue<dataType>(simplex);
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
            vector<dataType> field2 = getFieldValue<dataType>(candidate);

            if(field == field2){
                chosen = candidate;
                num_pairable++;
            }
        }
    }

    return num_pairable;
}

template <typename dataType>
void MultivariateGradient::readOutputCriticalCells(std::list<float>& points,
                                                          std::list<char>& dimensions){


    //print gradient
    for(int i=0;i<=dimensionality_;i++){

        for(int j=0; j<(*gradient)[i].size(); j++){
            Simplex simplex = Simplex(i,j);
            if(isCritical(simplex)){

                //is a critical cell according to our definition?

                vector<SimplexId> vertices;
                simplexToVertices(simplex, vertices);
                bool samelvlset = false;

                if(simplex.first != 0){

                    for(auto v : vertices){
                        vector<dataType> fieldV = getFieldValue<dataType>(Simplex(0,v));
                        vector<dataType> fieldS = getFieldValue<dataType>(simplex);

                        if(getFieldValue<dataType>(Simplex(0,v)) == getFieldValue<dataType>(simplex)){
                            samelvlset=true;
                            break;
                        }
                    }
                }

                if(!samelvlset){
                    //computeBarycenter
                    std::vector<float> barycenter;
                    computeBarycenter(simplex, barycenter);

                    for(auto b : barycenter)
                        points.push_back(b);

                    //dimensions
                    dimensions.push_back(i);
                }
            }
        }
    }
}

//#endif // MultivariateGradient_TEMPLATE_H
