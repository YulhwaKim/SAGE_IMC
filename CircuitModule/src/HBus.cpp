/*******************************************************************************
* Copyright (c) 2015-2017
* School of Electrical, Computer and Energy Engineering, Arizona State University
* PI: Prof. Shimeng Yu
* All rights reserved.
* 
* This source code is part of NeuroSim - a device-circuit-algorithm framework to benchmark 
* neuro-inspired architectures with synaptic devices(e.g., SRAM and emerging non-volatile memory). 
* Copyright of the model is maintained by the developers, and the model is distributed under 
* the terms of the Creative Commons Attribution-NonCommercial 4.0 International Public License 
* http://creativecommons.org/licenses/by-nc/4.0/legalcode.
* The source code is free and you can redistribute and/or modify it
* by providing that the following conditions are met:
* 
*  1) Redistributions of source code must retain the above copyright notice,
*     this list of conditions and the following disclaimer.
* 
*  2) Redistributions in binary form must reproduce the above copyright notice,
*     this list of conditions and the following disclaimer in the documentation
*     and/or other materials provided with the distribution.
* 
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
* SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
* 
* Developer list: 
*   Pai-Yu Chen	    Email: pchen72 at asu dot edu 
*                    
*   Xiaochen Peng   Email: xpeng15 at asu dot edu
********************************************************************************/

#include <cmath>
#include <iostream>
#include "constant.h"
#include "typedef.h"
#include "formula.h"
#include "HBus.h"
#include "Param.h"

using namespace std;

extern Param *param;

HBus::HBus(const InputParameter& _inputParameter, const Technology& _tech): inputParameter(_inputParameter), tech(_tech), FunctionUnit() {
	initialized = false;
}

void HBus::Initialize(BusMode _mode, bool _inputWire, int _numRow, int _numCol, double _delaytolerance, double _busWidth, 
                    double _unitHeight, double _unitWidth, double _foldedratio, double _clkFreq){
	if (initialized)
		cout << "[HBus] Warning: Already initialized!" << endl;

    mode = _mode;
    inputWire = _inputWire;
	numRow = _numRow;
	numCol = _numCol;     // num of Row and Col in tile/pe level

    unitHeight = _unitHeight;
    unitWidth = _unitWidth;

	delaytolerance = _delaytolerance;
	busWidth = _busWidth;

    foldedratio = _foldedratio;
	clkFreq = _clkFreq;
    switchingRatio = 0.25;

	unitLengthWireResistance = param->unitLengthWireResistance;
	unitLengthWireCap = 0.2e-15/1e-6;;   // 0.2 fF/mm
	
	// define min INV resistance and capacitance to calculate repeater size
	widthMinInvN = MIN_NMOS_SIZE * tech.featureSize;
	widthMinInvP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	CalculateGateArea(INV, 1, widthMinInvN, widthMinInvP, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &hMinInv, &wMinInv);
	CalculateGateCapacitance(INV, 1, widthMinInvN, widthMinInvP, hMinInv, tech, &capMinInvInput, &capMinInvOutput);
	double resOnRep = CalculateOnResistance(widthMinInvN, NMOS, inputParameter.temperature, tech) + CalculateOnResistance(widthMinInvP, PMOS, inputParameter.temperature, tech);


	// optimal repeater design to achieve highest speed
	repeaterSize = floor((double)sqrt( (double) resOnRep*unitLengthWireCap/capMinInvInput/unitLengthWireResistance));
	minDist = sqrt(2*resOnRep*(capMinInvOutput+capMinInvInput)/(unitLengthWireResistance*unitLengthWireCap));
	CalculateGateArea(INV, 1, MIN_NMOS_SIZE * tech.featureSize * repeaterSize, tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize * repeaterSize, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &hRep, &wRep);
	CalculateGateCapacitance(INV, 1, MIN_NMOS_SIZE * tech.featureSize * repeaterSize, tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize * repeaterSize, hRep, tech, &capRepInput, &capRepOutput);
	resOnRep = CalculateOnResistance(MIN_NMOS_SIZE * tech.featureSize * repeaterSize, NMOS, inputParameter.temperature, tech) + CalculateOnResistance(tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize * repeaterSize, PMOS, inputParameter.temperature, tech);
	double minUnitLengthDelay = 0.7*(resOnRep*(capRepInput+capRepOutput+unitLengthWireCap*minDist)+0.5*unitLengthWireResistance*minDist*unitLengthWireCap*minDist+unitLengthWireResistance*minDist*capRepInput)/minDist;
	double maxUnitLengthEnergy = (capRepInput+capRepOutput+unitLengthWireCap*minDist)*tech.vdd*tech.vdd/minDist;
	
	if (delaytolerance) {   // tradeoff: increase delay to decrease energy
		double delay = 0;
		double energy = 100;
		while(delay<minUnitLengthDelay*(1+delaytolerance) && (repeaterSize >= 1)) {
			repeaterSize -= 1;
			minDist *= 0.9;
			CalculateGateArea(INV, 1, MIN_NMOS_SIZE * tech.featureSize * repeaterSize, tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize * repeaterSize, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &hRep, &wRep);
			CalculateGateCapacitance(INV, 1, MIN_NMOS_SIZE * tech.featureSize * repeaterSize, tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize * repeaterSize, hRep, tech, &capRepInput, &capRepOutput);
			resOnRep = CalculateOnResistance(MIN_NMOS_SIZE * tech.featureSize * repeaterSize, NMOS, inputParameter.temperature, tech) + CalculateOnResistance(tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize * repeaterSize, PMOS, inputParameter.temperature, tech);
			delay = 0.7*(resOnRep*(capRepInput+capRepOutput+unitLengthWireCap*minDist)+0.5*unitLengthWireResistance*minDist*unitLengthWireCap*minDist+unitLengthWireResistance*minDist*capRepInput)/minDist;
			energy = (capRepInput+capRepOutput+unitLengthWireCap*minDist)*tech.vdd*tech.vdd/minDist;
		}
	}
	

	widthInvN = repeaterSize * MIN_NMOS_SIZE * tech.featureSize;
	widthInvP = repeaterSize * tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	
	// INV
	CalculateGateArea(INV, 1, widthInvN, widthInvP, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &hInv, &wInv);
    CalculateGateCapacitance(INV, 1, widthInvN, widthInvP, hInv, tech, &capInvInput, &capInvOutput);   

    // mode checking
    if (mode == HORIZONTAL) {
        hopLength = unitWidth;
        //wireLength = unitWidth * (numCol - 1);
        wireLength = unitWidth * (numCol - 0.5);
        //numBus = ceil(numRow / 2.0);
        numBus = numRow;
        //printf("%-20s %15.4e %15.e %5d\n", "HORIZONTAL", hopLength, wireLength, numBus);
    } else {
        hopLength = unitHeight;
        //wireLength = unitHeight * (numRow - 1);
        wireLength = unitHeight * (numRow - 0.5);
        //numBus = ceil(numCol / 2.0);
        numBus = numCol;
        //printf("%-20s %15.4e %15.e %5d\n", "VERTICAL", hopLength, wireLength, numBus);
    }

    numRepeater = busWidth * ceil(wireLength / minDist);

    // wireWidth
    wireWidth = busWidth * param->wireWidth * 1e-9 * 2; // 2x for spacing

    // input wire (required when the location of bus is not matched with the location of buffer)
    if ( inputWire ) {
        if ( mode == HORIZONTAL ) {
            inputHopLength = unitHeight;
            inputWireLength = unitHeight * (numRow - 1);
        } else {
            inputHopLength = unitWidth;
            inputWireLength = unitWidth * (numCol - 1);
        }
        numRepeaterInputWire = busWidth * ceil(inputWireLength/minDist);
    }
	
	initialized = true;
}

void HBus::CalculateArea() {
	if (!initialized) {
		cout << "[HBus] Error: Require initialization first!" << endl;
	} else {
        // area of wire ( repeater area + wire area )
		area = hInv * wInv * numRepeater + wireWidth * wireLength;
        // repeat for the number of bus
        area *= numBus;

        // area of inputWire
        if ( inputWire ) {
            area += hInv * wInv * numRepeaterInputWire + inputWireLength * wireWidth;
        }
	}
}

void HBus::CalculateLatency(int numHopsRow, int numHopsCol, double numRead){ // numHops = max idxRow/Col
	if (!initialized) {
		cout << "[HBus] Error: Require initialization first!" << endl;
	} else {

		readLatency = 0;
		
		double resOnRep = CalculateOnResistance(widthInvN, NMOS, inputParameter.temperature, tech) 
                        + CalculateOnResistance(widthInvP, PMOS, inputParameter.temperature, tech);
		unitLatencyRep = 0.7 * ( resOnRep * (capInvInput + capInvOutput + unitLengthWireCap * minDist)
                       + 0.5 * unitLengthWireResistance * minDist * unitLengthWireCap * minDist
                       + unitLengthWireResistance * minDist * capInvInput ) / minDist;
		unitLatencyWire = 0.7 * unitLengthWireResistance * minDist * unitLengthWireCap * minDist / minDist;
    
        // get numHops & idxBus
        int numHops, idxBus;
        if (mode == HORIZONTAL) {
            numHops = numHopsCol;
            //idxBus = 1 + floor( numHopsRow / 2 ) * 2;
            idxBus = numHopsRow;
        } else {
            numHops = numHopsRow;
            //idxBus = 1 + floor( numHopsCol / 2 ) * 2;
            idxBus = numHopsCol;
        }

        // latency of wire
		if (numRepeater > 0) {
			readLatency = hopLength * (numHops + 0.5) * unitLatencyRep;
		} else {
			readLatency = hopLength * (numHops + 0.5) * unitLatencyWire;
		}

        // latency of input wire
        if ( inputWire ) {
            if ( numRepeaterInputWire > 0 ) {
                readLatency += inputHopLength * idxBus * unitLatencyRep;
            } else {
                readLatency += inputHopLength * idxBus * unitLatencyWire;
            }
        }
		
		if (param->synchronous) {
			readLatency = ceil(readLatency*clkFreq);
		}
		readLatency *= numRead; 	
	}
}

void HBus::CalculatePower(int numHopsRow, int numHopsCol, double numBitAccess, double numRead) {
	if (!initialized) {
		cout << "[HBus] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;

        // leakage of wire ( repeater )
		repeaterLeakage = CalculateGateLeakage(INV, 1, widthInvN, widthInvP, inputParameter.temperature, tech) * tech.vdd;
        leakage += repeaterLeakage * numRepeater;
        leakage *= numBus;

        // leakage of inputWire repeater
        if ( inputWire ) {
            leakage += repeaterLeakage * numRepeaterInputWire;
        }

		unitLengthEnergyRep = (capInvInput + capInvOutput + unitLengthWireCap * minDist) 
                                * tech.vdd * tech.vdd / minDist * 0.25;
		unitLengthEnergyWire = (unitLengthWireCap * minDist) * tech.vdd * tech.vdd / minDist * 0.25;

        // get numHops & idxBus
        int numHops, idxBus;
        if (mode == HORIZONTAL) {
            numHops = numHopsCol;
            //idxBus = 1 + floor( numHopsRow / 2 ) * 2;
            idxBus = numHopsRow;
        } else {
            numHops = numHopsRow;
            //idxBus = 1 + floor( numHopsCol / 2 ) * 2;
            idxBus = numHopsCol;
        }

        // dynamicE of wire
		if (numRepeater > 0) {
			readDynamicEnergy = hopLength * (numHops + 0.5) * unitLengthEnergyRep;
		} else {
			readDynamicEnergy = hopLength * (numHops + 0.5) * unitLengthEnergyWire;
		}

        // dynamicE of inputWire
        if ( inputWire ) {
		    if (numRepeater > 0) {
		    	readDynamicEnergy += inputHopLength * idxBus * unitLengthEnergyRep;
		    } else {
		    	readDynamicEnergy += inputHopLength * idxBus * unitLengthEnergyWire;
		    }
        }

        readDynamicEnergy *= switchingRatio;
        readDynamicEnergy *= numBitAccess;
		readDynamicEnergy *= numRead;
	}
}

void HBus::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}

