#define FG_SEGMENTATION_TEMPLATE_H

#include "FG_Segmentation.h"

using namespace std;
using namespace ttk;

template <class dataType> void ttk::FG_Segmentation::computeIndexing(){


    dataType* field = (dataType*)inputData_;

    int vertices = triangulation_->getNumberOfVertices();

    vector<pair<dataType,SimplexId> > thepairs(vertices);
    for(int i=0; i<vertices; i++){
        thepairs[i] = pair<dataType,SimplexId>(field[i],i);
    }

    sort(thepairs.begin(), thepairs.end());

    indexing_ = new vector<SimplexId>(vertices);
    for(SimplexId i=0; i<vertices; i++){
        (*indexing_)[thepairs[i].second]=i;
    }

    computeGradient(indexing_);
}