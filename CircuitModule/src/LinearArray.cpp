#include <cmath>
#include <iostream>
#include "constant.h"
#include "typedef.h"
#include "formula.h"
#include "LinearArray.h"
#include "Param.h"

using namespace std;

extern Param *param;

LinearArray::LinearArray(const InputParameter& _inputParameter, const Technology& _tech): inputParameter(_inputParameter), tech(_tech), dff(_inputParameter, _tech), FunctionUnit() {
	initialized = false;
}

void LinearArray::Initialize(BusMode _mode, bool _inputWire, int _numRow, int _numCol, double _delaytolerance,
                            double _busWidth, double _unitHeight, double _unitWidth, double _foldedratio, double _clkFreq){
	if (initialized)
		cout << "[LinearArray] Warning: Already initialized!" << endl;

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
    switchingRatio = 0.5;

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
        wireLength = unitWidth;
        numNode = numCol;
        numBranch = numRow;
    } else {
        wireLength = unitHeight;
        numNode = numRow;
        numBranch = numCol;
    }

    numRepeater = busWidth * ceil(wireLength/minDist);

    // wireWidth
    wireWidth = busWidth * param->wireWidth * 1e-9 * 2; // 2x for spacing


    // initialize dff for each connection
    dff.Initialize((int)busWidth, clkFreq);

    // input wire (required when the location of bus is not matched with the location of buffer)
    if ( inputWire ) {
        if ( mode == HORIZONTAL ) {
            inputWireLength = unitHeight * (numRow - 1);
        } else { 
            inputWireLength = unitWidth * (numCol - 1);
        }
        numRepeaterInputWire = busWidth * ceil(inputWireLength/minDist);
    }
	
	initialized = true;
}

void LinearArray::CalculateArea() {
	if (!initialized) {
		cout << "[LinearArray] Error: Require initialization first!" << endl;
	} else {
        // area of wire ( repeater area + wire area )
        area = hInv * wInv * numRepeater + wireWidth * wireLength;

        //// area of DFF -- assume that the linear array use dff of each submodule for buffering & work as systolic array
        //if (mode == HORIZONTAL) {
        //    dff.CalculateArea(unitHeight, NULL, NONE);
        //} else {
        //    dff.CalculateArea(NULL, unitWidth, NONE);
        //}
        //area += dff.area;

        // reapeat for the number of connection
        area *= ( numNode - 1 ) * numBranch;

        // area of input wire
        if ( inputWire ) {
            area += hInv * wInv * numRepeaterInputWire + inputWireLength * wireWidth;
        }
    }
}

void LinearArray::CalculateLatency(int numActiveRow, int numActiveCol, double numRead){
	if (!initialized) {
		cout << "[LinearArray] Error: Require initialization first!" << endl;
	} else {
		readLatency = 0;
		
		double resOnRep = CalculateOnResistance(widthInvN, NMOS, inputParameter.temperature, tech) 
                        + CalculateOnResistance(widthInvP, PMOS, inputParameter.temperature, tech);
		unitLatencyRep = 0.7 * ( resOnRep * (capInvInput + capInvOutput + unitLengthWireCap * minDist)
                       + 0.5 * unitLengthWireResistance * minDist * unitLengthWireCap * minDist
                       + unitLengthWireResistance * minDist * capInvInput ) / minDist;
		unitLatencyWire = 0.7 * unitLengthWireResistance * minDist * unitLengthWireCap * minDist / minDist;

        /* NOTE: Systolic array does not care about the case that the data is not fully filled.
        * For layer-by-layer operation, this corner case does not influence the latency info that much. 
        * However, this corner case should be considered when multiple-layer is operated on the chip simulaneously, 
        * so that we cannot guarantee that the data is fed every cycle. 
        */
        // latency of wire
		if (numRepeater > 0) {
			readLatency = wireLength * unitLatencyRep;
		} else {
			readLatency = wireLength * unitLatencyWire;
		}

        // latency of input wire
        if ( inputWire ) {
            if ( numRepeaterInputWire > 0 ) {
			    readLatency += inputWireLength * unitLatencyRep * numActiveRow;
            } else {
			    readLatency += inputWireLength * unitLatencyWire * numActiveRow;
            }
        }
    
		if (param->synchronous) {
			readLatency = ceil(readLatency * clkFreq);
            //printf("Sys unit latency: %15.4e \n", readLatency);
		} else {
            //// latency of DFF
            //dff.CalculateLatency(1e20, 1);
            //readLatency += dff.readLatency;
        }

		readLatency *= numRead; 	
	}
}

void LinearArray::CalculatePower(int numActiveRow, int numActiveCol, double numBitAccess, double numRead) {
	if (!initialized) {
		cout << "[LinearArray] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;

        //// leakage & dynamicE of dff
        //dff.CalculatePower(1, (int)numBitAccess, param->validated);
        //leakage = dff.leakage;
        //readDynamicEnergy = dff.readDynamicEnergy;

        // leakage of wire ( repeater )
		repeaterLeakage = CalculateGateLeakage(INV, 1, widthInvN, widthInvP, inputParameter.temperature, tech) * tech.vdd;
		leakage += repeaterLeakage * numRepeater;
        leakage *= ( numNode - 1) * numBranch;

        // leakage of inputWire repeater
        if ( inputWire ) {
            leakage += repeaterLeakage * numRepeaterInputWire;
        }
        

		unitLengthEnergyRep = (capInvInput + capInvOutput + unitLengthWireCap * minDist) 
                                * tech.vdd * tech.vdd / minDist * 0.25;
		unitLengthEnergyWire = (unitLengthWireCap * minDist) * tech.vdd * tech.vdd / minDist * 0.25;
	
        // dynamicE of wire
		if (numRepeater > 0) {
			readDynamicEnergy = wireLength * unitLengthEnergyRep * numBitAccess;
		} else {
			readDynamicEnergy = wireLength * unitLengthEnergyWire * numBitAccess;
		}

        // get how many hops do each data cross
        int numActiveNode, numActiveBranch;
        if (mode == HORIZONTAL) {
            numActiveNode = numActiveCol;
            numActiveBranch = numActiveRow;
        } else {
            numActiveNode = numActiveRow;
            numActiveBranch = numActiveCol;
        }
        readDynamicEnergy *= ( numActiveNode - 1 ) * numActiveBranch;

        // dynamicE of inputWire
        if ( inputWire ) {
		    if (numRepeater > 0) {
		    	readDynamicEnergy += inputWireLength * unitLengthEnergyRep * numActiveRow * numBitAccess;
		    } else {
		    	readDynamicEnergy += inputWireLength * unitLengthEnergyWire * numActiveRow * numBitAccess;
		    }
        }
        readDynamicEnergy *= switchingRatio;
		readDynamicEnergy *= numRead;
    }		
}

void LinearArray::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}



