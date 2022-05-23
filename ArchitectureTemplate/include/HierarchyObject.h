/*********************************************************************************************
* Definition of the Hierarchy Object
*********************************************************************************************/

#ifndef HIERARCHYOBJECT_H_
#define HIERARCHYOBJECT_H_

#include <vector>
#include "InputParameter.h"
#include "Technology.h"
#include "MemCell.h"
#include "DigitalElements.h"
#include "InterConnect.h"
#include "BufferUnit.h"
#include "HierarchyRoot.h"

using namespace std;

class HierarchyObject {
//private:
public:
    HierarchyObject(InputParameter& _inputParameter, Technology& _tech, MemCell& _cell, 
                    const int _hlevel, const HierarchyRoot* _rootObject, const HierarchyObject* _subObject,
                    const vector<double> _designObject);
    virtual ~HierarchyObject() {}
    InputParameter& inputParameter;
    Technology& tech;
    MemCell& cell;
    const int hlevel; // level of the hierarchy (hlevel 0: HierarchyRoot - CIMArray)
    const HierarchyRoot *rootObject;
    const HierarchyObject *subObject;
    const vector<double> designObject;
    vector<double> designDE;
    vector<double> designIC;
    vector<double> designBU;
    
    const int lenDesignDE = 8;
    const int lenDesignIC = 5;
    const int lenDesignBU = 7;

    /* Components */
    DigitalElements *digitalElements;
    InterConnect *interConnect;
    BufferUnit *bufferUnit;

    /* Functions */
    void Initialize(double _clkFreq);
    void CalculateArea();
    void CalculateLatency(const vector<double> infoRead, const vector<double> subLatencyVector);
    void CalculatePower(const vector<double> infoRead, const vector<double> subReadDynamicEnergyVector);
    void CalculateICLatency(const vector<double> infoRead);
    void CalculateICPower(const vector<double> infoRead);
    void CalculateLeakage();
    //void PrintProperty(const char* str);

    const int lenInfoReadDE = 4;
    const int lenInfoReadIC = 6;
    const int lenInfoReadBU = 6;


    /* Properties */
    bool initialized;
    bool inputBuffer; // true: additional buffer for input , false: input&output share single buffer
   
    int numRow, numCol, numRowSubObject, numColSubObject; // number of row, col for this object / subobject
    int numInC, numOutC, numInCSubObject, numOutCSubObject; // number of channels that this object / subobject can handle
    int numSubObject, numSubObjectRow, numSubObjectCol; // number of subObject lie in the row/col dim
    double subWidth, subHeight;
    double height, width, area;
    double clkFreq;
    int numOutBit, numOutBitSubObject;

    /* Performance */
    vector<double> areaVector; /* Vector (total, array, ADC, accum, buffer, ic, other) */
    vector<double> latencyVector; /* Vector (total, array, accum, buffer, ic, other) */
    vector<double> readDynamicEnergyVector; /* Vector (total, array, ADC, accum, buffer, ic, other) */
    double leakage;

    vector<double> areaVector2; /* Vector (total, subObject, accum, buffer, ic, other) */
    vector<double> latencyVector2; /* Vector (total, subObject, accum, buffer, ic, other) */
    vector<double> readDynamicEnergyVector2; /* Vector (total, subObject, accum, buffer, ic, other) */
    

}; /* class HierarchyObject */

#endif /* HIERARCHYOBJECT_H_ */
