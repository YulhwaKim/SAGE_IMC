#include <cmath>
#include <iostream>
#include "constant.h"
#include "typedef.h"
#include "formula.h"
#include "Mesh.h"
#include "Param.h"

using namespace std;

extern Param *param;

Mesh::Mesh(const InputParameter& _inputParameter, const Technology& _tech): inputParameter(_inputParameter), tech(_tech), mux(_inputParameter, _tech), dff(_inputParameter, _tech), FunctionUnit() {
	initialized = false;
}

void Mesh::Initialize(int _numPort, int _flitSize, int _numRow, int _numCol, double _delaytolerance,
                    double _unitHeight, double _unitWidth, double _foldedratio, double _clkFreq) {
	if (initialized)
		cout << "[Mesh] Warning: Already initialized!" << endl;

    numPort = _numPort; // numPort: 5 (1 PE, 4 connection) or numPort: 8 (4 PE, 4 connection)
    if ( !(numPort==5 || numPort==8) ) {
        cerr << "[Mesh] Mesh only support numPort 5 or 8!" << endl;
        exit(-1);
    }
    flitSize = _flitSize;
    numRow = _numRow;
    numCol = _numCol;

    unitHeight = _unitHeight;
    unitWidth = _unitWidth;

    delaytolerance = _delaytolerance;
    foldedratio = _foldedratio;
    switchingRatio = 0.25;

    clkFreq = _clkFreq;

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


    // Calculate the number of Router & wireLength
    if ( numPort == 5 ) {
        // numRouter
        numRouterRow = numRow; 
        numRouterCol = numCol;
        numRouter = numRouterRow * numRouterCol;
        // wireLength bw router
        wireLengthH = unitWidth;
        wireLengthV = unitHeight;
    } else if ( numPort == 8 ) {
        // numRouter
        numRouterRow = ceil(numRow / 2.0);
        numRouterCol = ceil(numCol / 2.0);
        numRouter = numRouterRow * numRouterCol;
        // wireLength bw router
        wireLengthH = 2 * unitWidth;
        wireLengthV = 2 * unitHeight;
    }

    // number of repeater
    numRepeaterH = flitSize * ceil(wireLengthH/minDist);
    numRepeaterV = flitSize * ceil(wireLengthV/minDist);

    // wireWidth
    wireWidthH = flitSize * param->wireWidth * 1e-9 * 2; // 2x for spacing
    wireWidthV = flitSize * param->wireWidth * 1e-9 * 2; // 2x for spacing

    // DFF
    numDff = 4; 
    dff.Initialize(flitSize, clkFreq); // NOTE:BUG: change to use DFF or not based on the buffer usage of subObject
    // mux for crossbar switch
    numMux = numPort * flitSize;
    mux.Initialize(numPort - 1);
	
	initialized = true;
}

void Mesh::CalculateArea() {
	if (!initialized) {
		cout << "[Mesh] Error: Require initialization first!" << endl;
	} else {
        area = 0;

        // area of wire ( repeater area + wire area )
        if (numRouterCol > 1) {        
            double repeaterAreaH = hInv * wInv * numRepeaterH;
            double wireAreaH = wireWidthH * wireLengthH;
            area += ( repeaterAreaH + wireAreaH ) * ( numRouterCol - 1 ) * numRouterRow * 2; // 2 for Rx & Tx wire
        }
        if (numRouterRow > 1) {
            double repeaterAreaV = hInv * wInv * numRepeaterV;
            double wireAreaV = wireWidthV * wireLengthV;
            area += ( repeaterAreaV + wireAreaV ) * ( numRouterRow - 1 ) * numRouterCol * 2; // 2 for Rx & Tx wire
        }

        // area of DFF
        dff.CalculateArea(NULL, NULL, NONE);
        area += dff.area * numDff * numRouter;

        // area of mux
        mux.CalculateArea(NULL, NULL, NONE);
        area += mux.area * numMux * numRouter;
    }
}

void Mesh::CalculateLatency(int numHopsRow, int numHopsCol, double numRead){
	if (!initialized) {
		cout << "[Mesh] Error: Require initialization first!" << endl;
	} else {
		readLatency = 0;
		
		double resOnRep = CalculateOnResistance(widthInvN, NMOS, inputParameter.temperature, tech) 
                        + CalculateOnResistance(widthInvP, PMOS, inputParameter.temperature, tech);
		unitLatencyRep = 0.7 * ( resOnRep * (capInvInput + capInvOutput + unitLengthWireCap * minDist)
                       + 0.5 * unitLengthWireResistance * minDist * unitLengthWireCap * minDist
                       + unitLengthWireResistance * minDist * capInvInput ) / minDist;
		unitLatencyWire = 0.7 * unitLengthWireResistance * minDist * unitLengthWireCap * minDist / minDist;

        /* NOTE: Router does not consider bw blocking*/
        // latency of wire
		if (numRepeaterH > 0) {
			readLatency += wireLengthH * unitLatencyRep * numHopsCol;
		} else {
			readLatency += wireLengthH * unitLatencyWire * numHopsCol;
		}

		if (numRepeaterV > 0) {
			readLatency += wireLengthV * unitLatencyRep * numHopsRow;
		} else {
			readLatency += wireLengthV * unitLatencyWire * numHopsRow;
		}
    
		if (param->synchronous) {
			readLatency = ceil(readLatency * clkFreq);
		} else {
            // latency of DFF
            dff.CalculateLatency(1e20, 1);
            readLatency += dff.readLatency * (numHopsRow + numHopsCol);
        }

		readLatency *= numRead; 	
	}
}

void Mesh::CalculatePower(int numHopsRow, int numHopsCol, double numBitAccess, double numRead) {
	if (!initialized) {
		cout << "[Mesh] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;

        // leakage & dynamicE of dff
        dff.CalculatePower(1, (int)numBitAccess, param->validated);
        leakage = dff.leakage * numDff * numRouter;
        readDynamicEnergy = dff.readDynamicEnergy * (numHopsRow + numHopsCol);

        // leakage & dynamicE of mux
        mux.CalculatePower(numRead);
        leakage += mux.leakage * numMux * numRouter;
        readDynamicEnergy += mux.readDynamicEnergy * (numHopsRow + numHopsCol);

        // leakage of wire
		double repeaterLeakage = CalculateGateLeakage(INV, 1, widthInvN, widthInvP, inputParameter.temperature, tech) * tech.vdd;
        if (numRouterCol > 0) {        
            leakage += repeaterLeakage * numRepeaterH * (numRouterCol - 1) * numRouterRow * 2; // 2 for Rx & Tx wire
        }
        if (numRouterRow > 0) {
            leakage += repeaterLeakage * numRepeaterV * (numRouterRow - 1) * numRouterCol * 2; // 2 for Rx & Tx wire
        }
	
        // dynamicE of wire
		unitLengthEnergyRep = (capInvInput + capInvOutput + unitLengthWireCap * minDist)
                            * tech.vdd * tech.vdd / minDist * 0.25;
		unitLengthEnergyWire = (unitLengthWireCap * minDist) * tech.vdd * tech.vdd / minDist * 0.25;

		if (numRepeaterH > 0) {
			readDynamicEnergy += wireLengthH * unitLengthEnergyRep * numBitAccess * numHopsCol;
		} else {
			readDynamicEnergy += wireLengthH * unitLengthEnergyWire * numBitAccess * numHopsCol;
		}
		if (numRepeaterV > 0) {
			readDynamicEnergy += wireLengthV * unitLengthEnergyRep * numBitAccess * numHopsRow;
		} else {
			readDynamicEnergy += wireLengthV * unitLengthEnergyWire * numBitAccess * numHopsRow;
		}
        
        readDynamicEnergy *= switchingRatio; // BUG: what is the right switching ratio for mux???
		readDynamicEnergy *= numRead;
    }		
}

void Mesh::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}



