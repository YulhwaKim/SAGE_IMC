/*********************************************************************************************
* Definition of the Hierarchy Root
*********************************************************************************************/

#ifndef HIERARCHYROOT_H_
#define HIERARCHYROOT_H_

#include <vector>
#include "InputParameter.h"
#include "Technology.h"
#include "MemCell.h"
#include "CIMArray.h"

using namespace std;

class HierarchyRoot {
//private:
public:
    HierarchyRoot(InputParameter& _inputParameter, Technology& _tech, MemCell& _cell);
    virtual ~HierarchyRoot() {}
    InputParameter& inputParameter;
    Technology& tech;
    MemCell& cell;
    int hlevel = 0; // hlevel is 0

    /* Components */
    CIMArray * cimArray;

    /* Functions */
    void Initialize();
    void CalculateArea();
    void CalculateLatency(vector<double> infoReadCIM);
    void CalculateLatency(bool CalculateFreq, double weightMatrixRow, double weightMatrixCol, 
                        double numBitInput, double numCellPerSynapse, double numRead, double *clkPeriod);
    void CalculatePower(vector<double> infoReadCIM);
    void CalculatePower(double weightMatrixRow, double weightMatrixCol, 
                        double numBitInput, double numCellPerSynapse, double numRead);
    void CalculateLeakage();
    void PrintProperty();
    void GetColumnResistance();

    /* properties */
    bool initialized;   /* Initialization flag */

    int numRow, numCol;
    double height, width, area, usedArea;
    double maxConductance, minConductance;
    double numOutBit;

    double clkFreq;

    double inputActiveRatio; 
    vector<double> weightLevelRatioVector; /* ratio of each level (low -> high)*/
    double columnRes;

    /* Performance */
    vector<double> areaVector; /* Vector (total, array, ADC, accum, dummy(buffer, ic, other)) */
    vector<double> latencyVector; /* Vector (total, array (include ADC), accum, dummy(buffer, ic, other)) */
    vector<double> readDynamicEnergyVector; /* Vector (total, array, ADC, accum, dummy(buffer, ic, other)) */
    double leakage;
    

}; /* class HierarchyRoot */

#endif /* HIERARCHYROOT_H_ */
