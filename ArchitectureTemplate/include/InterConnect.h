/*********************************************************************************************
* Definition of InterConnect
*********************************************************************************************/

#ifndef INTERCONNECT_H_
#define INTERCONNECT_H_

#include <vector>
#include "InputParameter.h"
#include "Technology.h"
#include "Bus.h"
#include "LinearArray.h"
#include "Mesh.h"
#include "HBus.h"

using namespace std;

class InterConnect {
//private
public:
    InterConnect(const InputParameter& _inputParameter, const Technology& _tech, const vector<double> _designIC);
    virtual ~InterConnect() {}
    const InputParameter& inputParameter;
    const Technology& tech;
    const vector<double> designIC; // delaytolerance, outType, outBusWidth, inType, inBusWidth (delay, 2(2D Mesh), flit, port, -)

    /* Components */
    Bus             *outBus;
    Bus             *inBus;
    LinearArray     *inLinear;
    LinearArray     *outLinear;
    Mesh            *mesh;
    HBus            *outHBus;
    HBus            *inHBus;

    /* Functions */
    void Initialize(int numRow, int numCol, double _unitHeight, double _unitWidth, int inBUSize, double _clkFreq);
    void CalculateArea();
    /*inforReadIC: numOutRead, numInRead, x_init, y_init, x_end, y_end*/
    void CalculateLatency(const vector<double> infoReadIC); 
    void CalculatePower(const vector<double> inforReadIC);
    void CalculateLeakage();
    void PrintProperty(const char* str);

    /* properties */
    bool initialized;   /* Initialization flag */
    int inType, outType; /* 0: bus, 1:LinearArray, 2: 2D mesh NoC, 3: hierarchical bus*/
    double delaytolerance, outBusWidth, inBusWidth;
    int numPort, flitSize; /* used for mesh NoC only */
    double unitHeight, unitWidth;
    double clkFreq;
    BusMode inBusMode;

    /* Performance */
    double area;
    double readLatency;
    double leakage;
    double readDynamicEnergy;
    double overlapLatency, inputCounter, prevInputCounter; // for input latency overlap

}; /* class InterConnect */

#endif /* INTERCONNECT_H_ */
