/*********************************************************************************************
* Definition of the Layer Scheduler (non-compactt mapping)
*********************************************************************************************/

#ifndef LAYERSCHEDULER_H_
#define LAYERSCHEDULER_H_

#include <vector>
#include "HierarchyRoot.h"
#include "HierarchyObject.h"

using namespace std;

class LayerScheduler {
//private:
public:
    LayerScheduler();
    virtual ~LayerScheduler() {}

    /* Functions */
    void Initialize(int _layerIdx, const vector<double> _layerStructure,
                    int _hlevelMappingUnit, vector<vector<int>> _idxOffsetVector,
                    const HierarchyObject* hTop);
    void CheckLinearArray(const HierarchyObject* hObject);

    vector<double> HRootScheduling(const HierarchyRoot* hRoot, double idxRow, double idxCol,
                                double weightMatrixRow, double weightMatrixCol);
    // scheduler type00 - base (tile-wise mapping)
    vector<vector<double>> HObjectScheduling_00(const HierarchyObject* hObject, double idxRow, double idxCol,
                                double wRow, double wCol, double wInC, double wOutC,
                                bool offsetObject, bool lastObject);
    // scheduler type01 - compact
    vector<vector<double>> HObjectScheduling_01(const HierarchyObject* hObject, double idxRow, double idxCol,
                                double numRowObjectAvailable, double numColObjectAvailable,
                                double wRow, double wCol, double wInC, double wOutC);

    /* Scheduling Parameter (Architecture/Network info) */
    int numBitInput, numColMuxed, numCellPerSynapse;
    int lengthInfoRead;
    int hlevelRowSystolic, hlevelColSystolic;
    int numHObjectRowSAExt, numHObjectColSAExt; // number of hObject lie in the row/col dim for Systolic Array Extension
                                                // NOTE: assume that the extension can be finished in the parent hObject

    /* Layer Properties */
    int layerIdx;
    vector<double> layerStructure;
    double inW, inH, inC, kW, kH, outC, fanIn, fanOut, outW, outH, numConv;
    int hlevelMappingUnit; // hLevel used for mapping unit
    int hlevelTop;
    vector<vector<int>> idxOffsetVector; // hlevel, idxOffsetRow, idxOffsetCol
    vector<vector<int>> nextIdxOffsetVector; // hlevel, idxOffsetRow, idxOffsetCol

    /* Scheduling Progress Checker */
    bool doingAct, doingMaxPool;
    bool doneAct, doneMaxPool;
 
    /* Scheduling Result (for each object)*/
    int numSubObject, numUsedSubObjectRow, numUsedSubObjectCol; // number of object
    int numUsedSubObject_top;


}; /* class LayerScheduler */

#endif /* LAYERSCHEDULER_H_ */
