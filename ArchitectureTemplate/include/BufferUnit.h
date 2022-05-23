/*********************************************************************************************
* Definition of Buffer Unit
*********************************************************************************************/

#ifndef BUFFERUNIT_H_
#define BUFFERUNIT_H_

#include <vector>
#include "InputParameter.h"
#include "Technology.h"
#include "MemCell.h"
#include "DFF.h"
#include "Buffer.h"

class BufferUnit {
//private
public:
    BufferUnit(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell, const vector<double> _designBU);
    virtual ~BufferUnit() {}
    const InputParameter& inputParameter;
    const Technology& tech;
    const MemCell& cell;
    const vector<double> designBU; // buType, outBUSize, outBUCoreBW, numOutBUCore, inputBUSize, inputBUCoreBW, numInBUCore
                                /* buType: 0 - DFF, 1 - register file, 2 - SRAM */

    /* Components */
    DFF *outDff;
    DFF *inDff;
    Buffer *outBuffer;
    Buffer *inBuffer;

    /* Functions */
    void Initialize(double _unitWireRes, double _clkFreq);
    void CalculateArea(double newHeight, double newWidth);
    /* infoReadBU: numOutRead, numOutWrite, outParallelism, numInRead, numInWrite, inParallelism */
    void CalculateLatency(const vector<double> infoReadBU);
    void CalculatePower(const vector<double> inforReadBU);
    void CalculateLeakage();
    void PrintProperty(const char* str);

    /* properties */
    bool initialized;   /* Initialization flag */
    int buType; /* buType: 0 - DFF, 1 - register file, 2 - SRAM */
    int outBUSize, inBUSize, outBUCoreBW, inBUCoreBW, numOutBUCore, numInBUCore;
    double clkFreq;

    /* Performance */
    double area;
    double latency;
    double leakage;
    double dynamicEnergy;

}; /* class BufferUnit */

#endif /* BUFFEREUNIT_H_ */
