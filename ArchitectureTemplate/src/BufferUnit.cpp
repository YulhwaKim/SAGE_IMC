#include <iostream>
#include <cmath>
#include "constant.h"
#include "formula.h"
#include "BufferUnit.h"


BufferUnit::BufferUnit(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell, const vector<double> _designBU):
inputParameter(_inputParameter), tech(_tech), cell(_cell), designBU(_designBU) {

    buType = (int)designBU[0];
    outBUSize = (int)designBU[1];
    inBUSize = (int)designBU[4];
    outBUCoreBW = (int)designBU[2];
    inBUCoreBW = (int)designBU[5];
    numOutBUCore = (int)designBU[3];
    numInBUCore = (int)designBU[6];

    if ( buType == 0 ) {
        outDff = new DFF(inputParameter, tech);
        if ( inBUSize > 0 ) {
            inDff = new DFF(inputParameter, tech);
        }
    } else {
        outBuffer = new Buffer(inputParameter, tech, cell);
        if ( inBUSize > 0 ) {
            inBuffer = new Buffer(inputParameter, tech, cell);
        }
    }
    // not initialized
    initialized = false;

}

/* Initialize digital module */
void BufferUnit::Initialize(double _unitWireRes, double _clkFreq) {

    // set clock frequency
    clkFreq = _clkFreq;
 
    // initialize buffer
    if ( buType == 0 ) {
        outDff->Initialize(outBUSize, clkFreq);
        if ( inBUSize > 0 ) {
            inDff->Initialize(inBUSize, clkFreq);
        }
    } else {
        bool SRAM = (buType==2)? true : false;
        outBuffer->Initialize(outBUSize, (int)designBU[2]/*interface_width*/, (int)designBU[3]/*num_interface*/, 
                        _unitWireRes, clkFreq, SRAM);
        if ( inBUSize > 0 ) {
            inBuffer->Initialize(inBUSize, (int)designBU[5]/*interface_width*/, (int)designBU[6]/*num_interface*/, 
                            _unitWireRes, clkFreq, SRAM);
        }
        //// DICE debugging
        //printf("out Buffer Init\n");
        //printf("%-20s %10d\n", "outBUSize", (int)outBUSize);
        //printf("%-20s %10d\n", "interface_width", (int)designBU[2]);
        //printf("%-20s %10d\n", "num_interface", (int)designBU[3]);
    }
    // set initialized flag
    initialized = true;

}

/* Calculate area of digital module */
void BufferUnit::CalculateArea(double newHeight, double newWidth) {
   
    area = 0;

    if ( buType == 0 ) {
        outDff->CalculateArea(NULL, newWidth, NONE);
        area = outDff->area;
        if ( inBUSize > 0 ) {
            inDff->CalculateArea(newHeight, NULL, NONE);
            area += inDff->area;
        }
    } else {
        outBuffer->CalculateArea(NULL, newWidth, NONE);
        area = outBuffer->area;
        if ( inBUSize > 0 ) {
            inBuffer->CalculateArea(newHeight, NULL, NONE);
            area += inBuffer->area;
        }
        //// DICE debugging
        //printf("out Buffer Area\n");
        //printf("%-20s %10.4e\n", "newWidth", newWidth);
        //printf("%-20s %10.4e\n", "area", outBuffer->area);
    }

}

/* Calculate Latency of digital module */
void BufferUnit::CalculateLatency(const vector<double> infoReadBU) {

    latency = 0;

    if ( buType == 0 ) {
        outDff->CalculateLatency(0, infoReadBU[0]/*numOutRead*/);
        latency = outDff->readLatency;
        if ( inBUSize > 0 ) {
            inDff->CalculateLatency(0, infoReadBU[3]/*numInRead*/);
            latency += inDff->readLatency;
        }
    } else {
        outBuffer->CalculateLatency(outBuffer->interface_width, 
                                ceil(infoReadBU[0])/*numOutRead*/,  // ceil for access counting
                                outBuffer->interface_width, 
                                ceil(infoReadBU[1])/*numOutWrite*/);
        latency =  ( outBuffer->readLatency + outBuffer->writeLatency ) / infoReadBU[2];
        if ( inBUSize > 0 ) {
            inBuffer->CalculateLatency(inBuffer->interface_width, 
                                    ceil(infoReadBU[3])/*numInRead*/,  // ceil for access counting
                                    inBuffer->interface_width, 
                                    ceil(infoReadBU[4])/*numInWrite*/);
            latency +=  ( inBuffer->readLatency + inBuffer->writeLatency ) / infoReadBU[5];
        }
        latency = ceil(latency); // ceil for cycle counting 
        ////// DICE debugging 
        //printf("Buffer Latency breakdown (Cycle)\n");
        //printf("%-20s %15.1f\n", "outBuffer", outBuffer->readLatency + outBuffer->writeLatency);
        //printf("%-20s %15.1f\n", "outBuffer read", outBuffer->readLatency);
        //printf("%-20s %15.1f\n", "outBuffer write", outBuffer->writeLatency);
        //printf("%-20s %15d\n", "outBuffer interface_width", outBuffer->interface_width);
        //if ( inBUSize > 0 ) {
        //    printf("%-20s %15.1f\n", "inBuffer", inBuffer->readLatency + inBuffer->writeLatency);
        //    printf("%-20s %15d\n", "inBuffer interface_width", inBuffer->interface_width);
        //}
        //printf("%-20s %15.1f\n", "infoReadBU[2]", infoReadBU[2]);
        //printf("%-20s %15.1f\n", "infoReadBU[5]", infoReadBU[5]);
        //printf("%-20s %15.1f\n", "latency", latency);
    }

}

/* Calculate readDynamicEnergy & leakage of digital module */
void BufferUnit::CalculatePower(const vector<double> infoReadBU) {

    leakage = 0;
    dynamicEnergy = 0;

    if ( buType == 0 ) {
        outDff->CalculatePower(infoReadBU[0]/*numOutRead*/, outDff->numDff, false);
        leakage = outDff->leakage;
        dynamicEnergy = outDff->readDynamicEnergy;
        if ( inBUSize > 0 ) {
            inDff->CalculatePower(infoReadBU[3]/*numInRead*/, inDff->numDff, false);
            leakage += inDff->leakage;
            dynamicEnergy += inDff->readDynamicEnergy;
        }
    } else {
        outBuffer->CalculatePower(outBuffer->interface_width, 
                                infoReadBU[0] /* numOutRead */, 
                                outBuffer->interface_width,
                                infoReadBU[1] /* numOutWrite */);
        leakage = outBuffer->leakage;
        dynamicEnergy = outBuffer->readDynamicEnergy + outBuffer->writeDynamicEnergy;
        if ( inBUSize > 0 ) {
            inBuffer->CalculatePower(inBuffer->interface_width, 
                                    infoReadBU[3] /* numInRead */, 
                                    inBuffer->interface_width,
                                    infoReadBU[4] /* numInWrite */);
            leakage += inBuffer->leakage;
            dynamicEnergy += ( inBuffer->readDynamicEnergy + inBuffer->writeDynamicEnergy );
        }
        //// DICE debugging 
        //printf("Buffer Energy breakdown\n");
        //printf("%-20s %15.4e pJ\n", "outBuffer read energy", outBuffer->readDynamicEnergy*1e12);
        //printf("%-20s %15.4e pJ\n", "outBuffer write energy", outBuffer->writeDynamicEnergy*1e12);
        ////printf("%-20s %15d\n", "outBuffer interface_width", outBuffer->interface_width);
        //printf("%-20s %15.1f\n", "outBuffer numOutRead", infoReadBU[0]);
        //printf("%-20s %15.1f\n", "outBuffer numOutWrite", infoReadBU[1]);
        //if ( inBUSize > 0 ) {
        //    printf("%-20s %15.4e pJ\n", "inBuffer read energy", inBuffer->readDynamicEnergy*1e12);
        //    printf("%-20s %15.4e pJ\n", "inBuffer write energy", inBuffer->writeDynamicEnergy*1e12);
        //    //printf("%-20s %15d\n", "inBuffer interface_width", inBuffer->interface_width);
        //    printf("%-20s %15.1f\n", "inBuffer numOutRead", infoReadBU[3]);
        //    printf("%-20s %15.1f\n", "inBuffer numOutWrite", infoReadBU[4]);
        //}
    }

}

/* Calculate leakage of digital module */
void BufferUnit::CalculateLeakage() {

    leakage = 0;

    if ( buType == 0 ) {
        outDff->CalculatePower(0, 0, false);
        leakage = outDff->leakage;
        if ( inBUSize > 0 ) {
            inDff->CalculatePower(0, 0, false);
            leakage += inDff->leakage;
        }
    } else {
        outBuffer->CalculatePower(0, 0, 0, 0);
        leakage = outBuffer->leakage;
        if ( inBUSize > 0 ) {
            inBuffer->CalculatePower(0, 0, 0, 0);
            leakage += inBuffer->leakage;
        }
    }

}

void BufferUnit::PrintProperty(const char* str) {
    if ( buType == 0) {
        outDff->PrintProperty(str);
        if ( inBUSize > 0 ) {
            inDff->PrintProperty(str);
        }
    } else {
        outBuffer->PrintProperty(str);
        if ( inBUSize > 0 ) {
            inBuffer->PrintProperty(str);
        }
    }
}





