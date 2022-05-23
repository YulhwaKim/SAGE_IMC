/*********************************************************************************************
* Definition of the Digital Elements
*********************************************************************************************/

#ifndef DIGITALELEMENTS_H_
#define DIGITALELEMENTS_H_

#include <vector>
#include "InputParameter.h"
#include "Technology.h"
#include "AdderTree.h"
#include "BitShifter.h"
#include "MaxPooling.h"

using namespace std;

class DigitalElements {
//private:
public:
    DigitalElements(const InputParameter& _inputParameter, const Technology& _tech, const vector<double> _designDE);
    virtual ~DigitalElements() {}
    const InputParameter& inputParameter;
    const Technology& tech;
    const vector<double> designDE; // AdderTree - numUnit, numAdderBit, numAdd / reLU - numUnit, numBit / max - numUnit, numBit, window

    /* Components */
    AdderTree *adderTree;
    BitShifter *reLu;
    MaxPooling *maxPooling;

    /* Functions */
    void Initialize(bool _fixedDataFlow, double _clkFreq);
    void CalculateArea(double newWidth);
    /* infoReadDE: AdderTree - numRead, numUnitAdd / reLu - numRead / max -numRead */
    void CalculateLatency(const vector<double> infoReadDE); 
    void CalculatePower(const vector<double> infoReadDE);
    void CalculateLeakage();
    void PrintProperty(const char* str);

    /* properties */
    bool initialized;   /* Initialization flag */
    bool fixedDataFlow;

    bool placeAdderTree; /* placement of digital circuits */
    bool placeReLu;
    bool placeMaxPooling;

    double clkFreq;

    /* Performance Vector (total, addertree, relu, maxpool) */
    vector<double> areaVector;
    vector<double> latencyVector;
    vector<double> readDynamicEnergyVector;
    double leakage;

}; /* class DigitalElements */

#endif /* DIGITALELEMENTS_H_ */
