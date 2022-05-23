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
#include "formula.h"
#include "Param.h"
#include "DigitalMux.h"

using namespace std;

extern Param *param;

DigitalMux::DigitalMux(const InputParameter& _inputParameter, const Technology& _tech): inputParameter(_inputParameter), tech(_tech), FunctionUnit() {
	initialized = false;
}

void DigitalMux::Initialize(int _numSelection){ // composed of 2-to-1 Mux
	if (initialized)
		cout << "[DigitalMux] Warning: Already initialized!" << endl;

	numSelection = _numSelection;	/* Typically numMuxed */
    numStage  = (int)ceil(log2(numSelection));
    num2to1Mux = numSelection - 1;

	// DigitalMux
	widthTgN = MIN_NMOS_SIZE * tech.featureSize;
	widthTgP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
	EnlargeSize(&widthTgN, &widthTgP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech);

    // numInv, Nor, NAND required for each 2to1Mux
	numInv2to1Mux = 3;
    numNand2to1Mux = 2;
    numNor2to1Mux = 1;

    // numInv for control signal
    numInvControl = numStage;

    numInv = numInv2to1Mux * num2to1Mux + numInvControl;
    numNand = numNand2to1Mux * num2to1Mux;
    numNor = numNor2to1Mux * num2to1Mux;

	initialized = true;
}

void DigitalMux::CalculateArea(double _newHeight, double _newWidth, AreaModify _option){
	if (!initialized) {
		cout << "[DigitalMux] Error: Require initialization first!" << endl;
	} else {
		double hInv, wInv, hNand, wNand, hNor, wNor;
		area = 0;
		height = 0;
		width = 0;

		// gate area
		CalculateGateArea(INV, 1, widthTgN, widthTgP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hInv, &wInv);
		CalculateGateArea(NAND, 2, 2*widthTgN, widthTgP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hNand, &wNand);
		CalculateGateArea(NOR, 2, widthTgN, 2*widthTgP, tech.featureSize*MAX_TRANSISTOR_HEIGHT, tech, &hNor, &wNor);

        // overall area
		if (_newHeight && _option==NONE) { // Inv/Nand/Nor in multiple rows given the total height
			// Calculate the number of Inv/Nand/Nor per Col
			if ( (_newHeight < hInv) || (_newHeight < hNand) || (_newHeight < hNor) ) {
				cout << "[DigitalMux] Error: DigitalMux height is even larger than the assigned height !" << endl;
			} else {
                // INV
                int numColInv = 0; // Number of cols of INV
				int numInvPerCol = (int)ceil(_newHeight/hInv);
				if (numInvPerCol > numInv) {
					numInvPerCol = numInv;
				}
				numColInv = (int)ceil((double)numInv / numInvPerCol);

                // NAND
                int numColNand = 0; // Number of colss of NAND
				int numNandPerCol = (int)ceil(_newHeight/hNand);
				if (numNandPerCol > numNand) {
					numNandPerCol = numNand;
				}
				numColNand = (int)ceil((double)numNand / numNandPerCol);

                // NOR
                int numColNor = 0; // Number of Cols of NOR
				int numNorPerCol = (int)ceil(_newHeight/hNor);
				if (numNorPerCol > numNor) {
					numNorPerCol = numNor;
				}
				numColNor = (int)ceil((double)numNor / numNorPerCol);

				height = _newHeight;
				width = wInv * numColInv + wNand * numColNand + wNor * numColNor;
			}
        } else if (_newWidth && _option==NONE) { // Tg in multiple rows given the total width
			// Calculate the number of Tg per row
			if ( (_newWidth < wInv) || (_newWidth < wNand) || (_newWidth < wNor) ) {
				cout << "[DigitalMux] Error: DigitalMux width is even larger than the assigned width !" << endl;
			} else {
                // INV
                int numRowInv = 0; // Number of rows of INV
				int numInvPerRow = (int)ceil(_newWidth/wInv);
				if (numInvPerRow > numInv) {
					numInvPerRow = numInv;
				}
				numRowInv = (int)ceil((double)numInv / numInvPerRow);

                // NAND
                int numRowNand = 0; // Number of rows of NAND
				int numNandPerRow = (int)ceil(_newWidth/wNand);
				if (numNandPerRow > numNand) {
					numNandPerRow = numNand;
				}
				numRowNand = (int)ceil((double)numNand / numNandPerRow);

                // NOR
                int numRowNor = 0; // Number of rows of NOR
				int numNorPerRow = (int)ceil(_newWidth/wNor);
				if (numNorPerRow > numNor) {
					numNorPerRow = numNor;
				}
				numRowNor = (int)ceil((double)numNor / numNorPerRow);

				width = _newWidth;
				height = hInv * numRowInv + hNand * numRowNand + hNor * numRowNor;
			}
		} else {    // Assume one row of Tg by default
			width = ( wInv * numInv2to1Mux + wNand * numNand2to1Mux + wNor * numNor2to1Mux ) + wInv; // wInv for control
            height = MAX(hNand * num2to1Mux, hNor * num2to1Mux); 
		}

		area = height * width;

		// Modify layout
		newHeight = _newHeight;
		newWidth = _newWidth;
		switch (_option) {
			case MAGIC:
				MagicLayout();
				break;
			case OVERRIDE:
				OverrideLayout();
				break;
			default:    // NONE
				break;
		}


        // Cap
		CalculateGateCapacitance(INV, 1, widthTgN, widthTgP, hInv, tech, &capInvInput, &capInvOutput);
		CalculateGateCapacitance(NAND, 2, 2*widthTgN, widthTgP, hNand, tech, &capNandInput, &capNandOutput);
		CalculateGateCapacitance(NOR, 2, widthTgN, 2*widthTgP, hNor, tech, &capNorInput, &capNorOutput);
	}
}

void DigitalMux::CalculateLatency(double _rampInput, double numRead) { 
	if (!initialized) {
		cout << "[DigitalMux] Error: Require initialization first!" << endl;
	} else {
        if (param->synchronous) {
            readLatency = 0;        // combinationa -> no latency
        } else {
            readLatency = 0; // latency for non-sync case is not calculated here.
        } 
	}
}

void DigitalMux::CalculatePower(double numRead) {
	if (!initialized) {
		cout << "[DigitalMux] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;

		// Leakage power
		leakage += CalculateGateLeakage(INV, 1, widthTgN, widthTgP, inputParameter.temperature, tech) * tech.vdd * numInv;
		leakage += CalculateGateLeakage(NAND, 2, 2*widthTgN, widthTgP, inputParameter.temperature, tech) * tech.vdd * numNand;
		leakage += CalculateGateLeakage(NOR, 2, widthTgN, 2*widthTgP, inputParameter.temperature, tech) * tech.vdd * numNor;

		// Read dynamic energy (assuming only one addr is selected)
        // 2to1Mux
        readDynamicEnergy += ( (capInvInput + capInvOutput + capNandInput + capNandOutput) * 2 
                                + capNorInput + capNorOutput + capInvInput + capInvOutput ) * tech.vdd * tech.vdd * numStage;
        // inverters for control
		readDynamicEnergy += (capInvInput + capInvOutput) * tech.vdd * tech.vdd * numInvControl;
		
		readDynamicEnergy *= numRead;
		
		if(param->validated){
			readDynamicEnergy *= param->epsilon; 	// switching activity of control circuits, epsilon = 0.05 by default
		}
	}
}

void DigitalMux::PrintProperty(const char* str) {
	FunctionUnit::PrintProperty(str);
}

void DigitalMux::SaveOutput(const char* str) {
	FunctionUnit::SaveOutput(str);
}

