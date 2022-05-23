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
#include <vector>
#include "constant.h"
#include "formula.h"
#include "CIMArray.h"
#include "Param.h"


using namespace std;

extern Param *param;

CIMArray::CIMArray(InputParameter& _inputParameter, Technology& _tech, MemCell& _cell):
						inputParameter(_inputParameter), tech(_tech), cell(_cell),
						wllevelshifter(_inputParameter, _tech, _cell),
						sllevelshifter(_inputParameter, _tech, _cell),
						bllevelshifter(_inputParameter, _tech, _cell),
						wlDecoder(_inputParameter, _tech, _cell),
						wlDecoderOutput(_inputParameter, _tech, _cell),
						wlNewDecoderDriver(_inputParameter, _tech, _cell),
						wlNewSwitchMatrix(_inputParameter, _tech, _cell),
						rowCurrentSenseAmp(_inputParameter, _tech, _cell),
						mux(_inputParameter, _tech, _cell),
						muxDecoder(_inputParameter, _tech, _cell),
						slSwitchMatrix(_inputParameter, _tech, _cell),
						blSwitchMatrix(_inputParameter, _tech, _cell),
						wlSwitchMatrix(_inputParameter, _tech, _cell),
						deMux(_inputParameter, _tech, _cell),
						readCircuit(_inputParameter, _tech, _cell),
						precharger(_inputParameter, _tech, _cell),
						senseAmp(_inputParameter, _tech, _cell),
						wlDecoderDriver(_inputParameter, _tech, _cell),
						sramWriteDriver(_inputParameter, _tech, _cell),
						adder(_inputParameter, _tech),
						dff(_inputParameter, _tech),
						shiftAddInput(_inputParameter, _tech, _cell),
						shiftAddCell(_inputParameter, _tech, _cell),
						multilevelSenseAmp(_inputParameter, _tech, _cell),
						multilevelSAEncoder(_inputParameter, _tech, _cell),
						sarADC(_inputParameter, _tech, _cell) {
	initialized = false;
	readDynamicEnergyArray = writeDynamicEnergyArray = 0;
} 

void CIMArray::Initialize(int _numRow, int _numCol, double _unitWireRes){  //initialization module
	
	numRow = _numRow;    //import parameters
	numCol = _numCol;
	unitWireRes = _unitWireRes;
	
	double MIN_CELL_HEIGHT = MAX_TRANSISTOR_HEIGHT;  //set real layout cell height
	double MIN_CELL_WIDTH = (MIN_GAP_BET_GATE_POLY + POLY_WIDTH) * 2;  //set real layout cell width

	if (cell.memCellType == Type::RRAM ||  cell.memCellType == Type::FeFET) {  //if array is RRAM
		double cellHeight = cell.heightInFeatureSize; 
		double cellWidth = cell.widthInFeatureSize;  
		if (cell.accessType == CMOS_access) {  // 1T1R
			if (relaxArrayCellWidth) {
				lengthRow = (double)numCol * MAX(cellWidth, MIN_CELL_WIDTH*2) * tech.featureSize;	// Width*2 because generally switch matrix has 2 pass gates per column, even the SL/BL driver has 2 pass gates per column in traditional 1T1R memory
			} else {
				lengthRow = (double)numCol * cellWidth * tech.featureSize;
			}
			if (relaxArrayCellHeight) {
				lengthCol = (double)numRow * MAX(cellHeight, MIN_CELL_HEIGHT) * tech.featureSize;
			} else {
				lengthCol = (double)numRow * cellHeight * tech.featureSize;
			}
		} else {	// Cross-point, if enter anything else except 'CMOS_access'
            /* DICE: REMOVED */
		}
	}      //finish setting array size
	
	capRow1 = lengthRow * 0.2e-15/1e-6;	// BL for 1T1R, WL for Cross-point and SRAM
	capRow2 = lengthRow * 0.2e-15/1e-6;	// WL for 1T1R
	capCol = lengthCol * 0.2e-15/1e-6;
	
	resRow = lengthRow * unitWireRes; 
	resCol = lengthCol * unitWireRes;

    
    //// DICE: debugging
    //printf("CIMarray\n");
    //printf("%-20s %15.4e\n", "lengthRow", lengthRow); // 3.3792e-05
    //printf("%-20s %15.4e\n", "lengthCol", lengthCol); // 1.1264e-05

    numColPerRead = numCol / numColMuxed;
	
	//start to initializing the subarray modules
    if (cell.memCellType == Type::RRAM || cell.memCellType == Type::FeFET) {
		if (cell.accessType == CMOS_access) {	// 1T1R
			cell.resCellAccess = cell.resistanceOn * IR_DROP_TOLERANCE;    //calculate access CMOS resistance
			cell.widthAccessCMOS = CalculateOnResistance(tech.featureSize, NMOS, inputParameter.temperature, tech) 
                                * LINEAR_REGION_RATIO / cell.resCellAccess;   //get access CMOS width
			if (cell.widthAccessCMOS > cell.widthInFeatureSize) {	// Place transistor vertically
				printf("Transistor width of 1T1R=%.2fF is larger than the assigned cell width=%.2fF in layout\n", cell.widthAccessCMOS, cell.widthInFeatureSize);
				exit(-1);
			}

			cell.resMemCellOn = cell.resCellAccess + cell.resistanceOn;        //calculate single memory cell resistance_ON
			cell.resMemCellOff = cell.resCellAccess + cell.resistanceOff;      //calculate single memory cell resistance_OFF
			cell.resMemCellAvg = cell.resCellAccess + cell.resistanceAvg;      //calculate single memory cell resistance_AVG

			capRow2 += CalculateGateCap(cell.widthAccessCMOS * tech.featureSize, tech) * numCol;          //sum up all the gate cap of access CMOS, as the row cap
			capCol += CalculateDrainCap(cell.widthAccessCMOS * tech.featureSize, NMOS, cell.widthInFeatureSize * tech.featureSize, tech) * numRow;	// If capCol is found to be too large, increase cell.widthInFeatureSize to relax the limit
		} else {	// Cross-point
            /* DICE: REMOVED */
		}
		
		if (cell.writeVoltage > 1.5) {
			wllevelshifter.Initialize(numRow, clkFreq);
			bllevelshifter.Initialize(numRow, clkFreq);
			sllevelshifter.Initialize(numCol, clkFreq);
		}
		
		if (conventionalSequential) {  
            /* DICE: REMOVED */
		} else if (conventionalParallel) { 
			// double resTg = cell.resMemCellOn / numRow;
			double resTg = cell.resMemCellAvg / numRow;

			wlSwitchMatrix.Initialize(ROW_MODE, numRow, resRow, true, false, 
                                numWriteCellPerOperationMemory, numWriteCellPerOperationNeuro, 1, clkFreq);
			slSwitchMatrix.Initialize(COL_MODE, numCol, resTg, true, false, // resTg --> parallel connection of memCell
                                numWriteCellPerOperationMemory, numWriteCellPerOperationNeuro, 1, clkFreq);     
			if (numColMuxed>1) {
				mux.Initialize(numColPerRead, numColMuxed, resTg, FPGA);       
				muxDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numColMuxed)), true, false);
			}
			if (SARADC) {
				sarADC.Initialize(numColPerRead, levelOutput, clkFreq);
			} else {
				multilevelSenseAmp.Initialize(numColPerRead, levelOutput, clkFreq, true, currentMode);
				multilevelSAEncoder.Initialize(levelOutput, numColPerRead);
			}

            numOutBit = (int)ceil(log2(levelOutput));

			if (numReadPulse > 1) {
				shiftAddInput.Initialize(numColPerRead, numOutBit, 
                                    clkFreq, spikingMode, numReadPulse, numColMuxed);
                numOutBit = shiftAddInput.numOutBit;
			}
			if (numCellPerSynapse > 1) {
                int tmp_numBit = (numReadPulse > 1)? shiftAddInput.numOutBit : (int)ceil(log2(levelOutput));
                SpikingMode tmp_spkMode = (numReadPulse > 1)? NONSPIKING : spikingMode;
				shiftAddCell.Initialize(numColPerRead, tmp_numBit, 
                                    clkFreq, tmp_spkMode, numCellPerSynapse, 1);
                numOutBit = shiftAddCell.numOutBit;
			}
			
		} else if (BNNsequentialMode || XNORsequentialMode) {       
            /* DICE: REMOVED */
		} else if (BNNparallelMode || XNORparallelMode) {      
            /* DICE: REMOVED */
		}
	} 
	initialized = true;  //finish initialization
}



void CIMArray::CalculateArea() {  //calculate layout area for total design
	if (!initialized) {
		cout << "[Subarray] Error: Require initialization first!" << endl;  //ensure initialization first
	} else {  //if initialized, start to do calculation
		area = 0;
		usedArea = 0;

	    if (cell.memCellType == Type::RRAM || cell.memCellType == Type::FeFET) {
			// Array only
			heightArray = lengthCol;
			widthArray = lengthRow;
			areaArray = heightArray * widthArray;
			
			// Level shifter for write
			if (cell.writeVoltage > 1.5) {
				wllevelshifter.CalculateArea(heightArray, NULL, NONE);
				bllevelshifter.CalculateArea(heightArray, NULL, NONE);
				sllevelshifter.CalculateArea(NULL, widthArray, NONE);				
			}
			
			if (conventionalSequential) { 				
                /* DICE: REMOVED */
			} else if (conventionalParallel) { 
				wlSwitchMatrix.CalculateArea(heightArray, NULL, NONE);
				slSwitchMatrix.CalculateArea(NULL, widthArray, NONE);
				if (numColMuxed > 1) {
					mux.CalculateArea(NULL, widthArray, NONE);
					muxDecoder.CalculateArea(NULL, NULL, NONE);
					double minMuxHeight = MAX(muxDecoder.height, mux.height);
					mux.CalculateArea(minMuxHeight, widthArray, OVERRIDE);
				}
				if (SARADC) {
					sarADC.CalculateUnitArea();
					sarADC.CalculateArea(NULL, widthArray, NONE);
				} else {
					multilevelSenseAmp.CalculateArea(NULL, widthArray, NONE);
					multilevelSAEncoder.CalculateArea(NULL, widthArray, NONE);
				}
				if (numReadPulse > 1) {
					shiftAddInput.CalculateArea(NULL, widthArray, NONE);
				}
				if (numCellPerSynapse > 1) {
					shiftAddCell.CalculateArea(NULL, widthArray, NONE);
				}
				height = ( (cell.writeVoltage > 1.5)? sllevelshifter.height : 0 ) 
                        + slSwitchMatrix.height + heightArray 
                        + ( (numColMuxed > 1)? mux.height : 0 ) 
                        + multilevelSenseAmp.height + multilevelSAEncoder.height  
                        + shiftAddInput.height + shiftAddCell.height
                        + sarADC.height;
				width = MAX( ( (cell.writeVoltage > 1.5)? (wllevelshifter.width + bllevelshifter.width) : 0 ) 
                                + wlNewSwitchMatrix.width + wlSwitchMatrix.width
                            , ( (numColMuxed > 1)? muxDecoder.width : 0 ) )
                        + widthArray;
				usedArea = areaArray 
                        + ( (cell.writeVoltage > 1.5)? (wllevelshifter.area + bllevelshifter.area + sllevelshifter.area) : 0 ) 
                        + wlSwitchMatrix.area + wlNewSwitchMatrix.area + slSwitchMatrix.area 
                        + ( (numColMuxed > 1)? (mux.area + muxDecoder.area) : 0 ) 
                        + multilevelSenseAmp.area  + multilevelSAEncoder.area 
                        + shiftAddInput.area + shiftAddCell.area
                        + sarADC.area;
				
				areaADC = multilevelSenseAmp.area + multilevelSAEncoder.area + sarADC.area;
				areaAccum = shiftAddInput.area + shiftAddCell.area;
				areaOther = ((cell.writeVoltage > 1.5)? (wllevelshifter.area + bllevelshifter.area + sllevelshifter.area) : 0) 
                        + wlNewSwitchMatrix.area + wlSwitchMatrix.area + slSwitchMatrix.area 
                        + ( (numColMuxed > 1)? (mux.area + muxDecoder.area) : 0 );
				
				area = height * width;				
				emptyArea = area - usedArea;
                // assume that the layout can be magically done in the square
                height = sqrt(area);
                width = area/height;
			} else if (BNNsequentialMode || XNORsequentialMode) {    
                /* DICE: REMOVED */
			} else if (BNNparallelMode || XNORparallelMode) {      
                /* DICE: REMOVED */
			}
		} 
	}
}

/* read operation only*/
void CIMArray::CalculateLatency(double columnRes, double weightMatrixRow, double weightMatrixCol, 
                            double numBitInput, bool CalculateclkFreq) {   //calculate latency for different mode 
    if ( weightMatrixCol > (double)numCol ) {
        cerr << "[Error] weightMatrixCol of the CIMarray cannot exceed the numCol!!" << endl;
        exit(-1);
    }
	if (!initialized) {
		cout << "[Subarray] Error: Require initialization first!" << endl;
	} else {
		
		readLatency = 0;
		readLatencyADC = 0;
		readLatencyAccum = 0;
		readLatencyOther = 0;
		writeLatency = 0;

        double numMuxing = ceil(weightMatrixCol / numColPerRead); 
        double numCIM = numMuxing * numBitInput;

	    if (cell.memCellType == Type::RRAM || cell.memCellType == Type::FeFET) {
			if (conventionalSequential) {
                /* DICE: REMOVED */
			} else if (conventionalParallel) {
				double capBL = lengthCol * 0.2e-15/1e-6;
				double colRamp = 0;
				double tau = (capCol)*(cell.resMemCellAvg/(numRow/2));
				colDelay = horowitz(tau, 0, 1e20, &colRamp);
				colDelay = tau * 0.2;  // assume the 15~20% voltage drop is enough for sensing
				if (CalculateclkFreq || !param->synchronous) {				
					wlSwitchMatrix.CalculateLatency(1e20, capRow2, resRow, 1, 0);
					slSwitchMatrix.CalculateLatency(1e20, capCol, resCol, 1, 0);
					if (numColMuxed>1) {
						mux.CalculateLatency(colRamp, 0, 1);
						muxDecoder.CalculateLatency(1e20, mux.capTgGateN*numColPerRead, mux.capTgGateP*numColPerRead, 1, 0);
					}
					if (SARADC) {
						sarADC.CalculateLatency(1);
					} else {
						multilevelSenseAmp.CalculateLatency(columnRes, 1, 1);
						multilevelSAEncoder.CalculateLatency(1e20, 1);
					}				
					if (CalculateclkFreq) {
						readLatency += MAX( wlSwitchMatrix.readLatency + slSwitchMatrix.readLatency , 
                                            ((numColMuxed > 1)==true? (mux.readLatency+muxDecoder.readLatency):0));
						readLatency += colDelay;
						readLatency += multilevelSenseAmp.readLatency;
						readLatency += multilevelSAEncoder.readLatency;
						readLatency += sarADC.readLatency;
						readLatency *= (validated==true? param->beta : 1);	// latency factor of sensing cycle, beta = 1.4 by default									
                        //// DICE: debugging
                        //// slSwitchMatrix is the major source of the latency change
                        //printf("CIMarray CalculateLatency Breakdown\n");
                        //printf("%-20s %15.4e ns\n", "capRow2", capRow2);
                        //printf("%-20s %15.4e ns\n", "resRow", resRow);
                        //printf("%-20s %15.4e ns\n", "capCol", capCol);
                        //printf("%-20s %15.4e ns\n", "resCol", resCol);
                        //printf("%-20s %15.4e ns\n", "wlSwitchMatrix", wlSwitchMatrix.readLatency*1e9);
                        //printf("%-20s %15.4e ns\n", "slSwitchMatrix", slSwitchMatrix.readLatency*1e9); 
                        //printf("%-20s %15.4e ns\n", "mux", mux.readLatency*1e9);
                        //printf("%-20s %15.4e ns\n", "muxDecoder", muxDecoder.readLatency*1e9);
                        //printf("%-20s %15.4e ns\n", "colDelay", colDelay*1e9);
                        //printf("%-20s %15.4e ns\n", "multilevelSenseAmp", multilevelSenseAmp.readLatency*1e9);
                        //printf("%-20s %15.4e ns\n", "multilevelSAEncoder", multilevelSAEncoder.readLatency*1e9);
                        //printf("%-20s %15.4e ns\n", "sarADC", sarADC.readLatency*1e9);
                        //if ( validated ) {
                        //    printf("%-20s %15.4e\n", "beta", param->beta);
                        //}
					}
				}
				if (!CalculateclkFreq) {
					if (numReadPulse > 1) {
						shiftAddInput.CalculateLatency(numCIM);
					}
                    if (numCellPerSynapse > 1) {  
                        double numSA = (numCellPerSynapse > numMuxing)? numCellPerSynapse : numMuxing;
						shiftAddCell.CalculateLatency(numSA);	
                    }
					if (param->synchronous) {
                        //DICE: readLatencyADC only includes the ADC latency
						readLatencyOther = numCIM;
					} else {
                        //DICE: readLatencyADC only includes the ADC latency
						readLatencyADC = ( multilevelSenseAmp.readLatency + multilevelSAEncoder.readLatency 
                                            + sarADC.readLatency ) 
                                        * numCIM * (validated==true? param->beta : 1);
						readLatencyOther = colDelay * numCIM * (validated==true? param->beta : 1);
						readLatencyOther += MAX( wlSwitchMatrix.readLatency
                                                , ( (numColMuxed > 1)==true? (mux.readLatency + muxDecoder.readLatency) : 0 ) ) 
                                            * numCIM * (validated==true? param->beta : 1);
					}
					readLatencyAccum = shiftAddInput.readLatency + shiftAddCell.readLatency;
					readLatency = readLatencyADC + readLatencyAccum + readLatencyOther;
				}
			} else if (BNNsequentialMode || XNORsequentialMode) {
                /* DICE: REMOVED */
			} else if (BNNparallelMode || XNORparallelMode) {
                /* DICE: REMOVED */
			}
		}
	}
}

/* read operation only*/
void CIMArray::CalculatePower(double columnRes, double weightMatrixRow, double weightMatrixCol, double numBitInput,
                            double inputActiveRatio) {
	if (!initialized) {
		cout << "[Subarray] Error: Require initialization first!" << endl;
	} else {
		readDynamicEnergy = 0;
		writeDynamicEnergy = 0;
		readDynamicEnergyArray = 0;

        double activityRowRead = weightMatrixRow / numRow * inputActiveRatio;
        double activityColRead = weightMatrixCol / numCol;
		
        double numMuxing = ceil(weightMatrixCol / numColPerRead);
        double numCIM = numMuxing * numBitInput;

	    if (cell.memCellType == Type::RRAM || cell.memCellType == Type::FeFET) {
			if (conventionalSequential) {
                /* DICE: REMOVED */
			} else if (conventionalParallel) {

				double capBL = lengthCol * 0.2e-15/1e-6;
                double activityColRead = weightMatrixCol / numCol;
			
				wlSwitchMatrix.CalculatePower(numBitInput, 0, activityRowRead, 0);
				slSwitchMatrix.CalculatePower(1, 0, activityColRead, 0);
				if (numColMuxed > 1) {
					mux.CalculatePower(numCIM);	// Mux still consumes energy during row-by-row read
					muxDecoder.CalculatePower(numCIM, 1);
				}
				if (SARADC) {
					sarADC.CalculatePower(columnRes, weightMatrixCol, numBitInput);
				} else {
					multilevelSenseAmp.CalculatePower(columnRes, weightMatrixCol, numBitInput);
					multilevelSAEncoder.CalculatePower(numCIM);
				}
				if (numReadPulse > 1) {
					shiftAddInput.CalculatePower(numCIM);
				}
				if (numCellPerSynapse > 1) {
                    double numSA = (numCellPerSynapse > numMuxing)? numCellPerSynapse : numMuxing;
					shiftAddCell.CalculatePower(numSA);
				}
				// Read
				readDynamicEnergyArray = 0;
				readDynamicEnergyArray += capBL * cell.readVoltage * cell.readVoltage * weightMatrixCol; // Selected BLs
				readDynamicEnergyArray += capRow2 * tech.vdd * tech.vdd * numRow * activityRowRead; // Selected WL
                readDynamicEnergyArray *= numBitInput;  // repeat numBitInput times

                //// DICE: debugging
                //printf("CIMarray readDynamicEnergyArray Breakdown\n");
                //printf("%-20s %15.4e\n", "capBL", capBL);
                //printf("%-20s %15.4e\n", "capRow2", capRow2);
                //printf("%-20s %15.4f\n", "cell.readVoltage", (double)cell.readVoltage);
                //printf("%-20s %15.4f\n", "tech.vdd", (double)tech.vdd);
                //printf("%-20s %15.4f\n", "weightMatrixCol", (double)weightMatrixCol);
                //printf("%-20s %15.4f\n", "numRow", (double)numRow);
                //printf("%-20s %15.4f\n", "activityRowRead", (double)activityRowRead);
                //printf("%-20s %15.4f\n", "numBitInput", (double)numBitInput);
				
				readDynamicEnergy = 0;
				readDynamicEnergy += wlSwitchMatrix.readDynamicEnergy;
				readDynamicEnergy += slSwitchMatrix.readDynamicEnergy;
				readDynamicEnergy += ( (numColMuxed > 1)==true? (mux.readDynamicEnergy + muxDecoder.readDynamicEnergy) : 0 );
				readDynamicEnergy += multilevelSenseAmp.readDynamicEnergy;
				readDynamicEnergy += multilevelSAEncoder.readDynamicEnergy;
				readDynamicEnergy += shiftAddInput.readDynamicEnergy;
				readDynamicEnergy += shiftAddCell.readDynamicEnergy;
				readDynamicEnergy += readDynamicEnergyArray;
				readDynamicEnergy += sarADC.readDynamicEnergy;
				
				readDynamicEnergyADC = multilevelSenseAmp.readDynamicEnergy + multilevelSAEncoder.readDynamicEnergy 
                                    + sarADC.readDynamicEnergy;
				readDynamicEnergyAccum = shiftAddInput.readDynamicEnergy + shiftAddCell.readDynamicEnergy;
				readDynamicEnergyOther = wlSwitchMatrix.readDynamicEnergy + slSwitchMatrix.readDynamicEnergy
                                    + ( (numColMuxed > 1)==true? (mux.readDynamicEnergy + muxDecoder.readDynamicEnergy) : 0 );

                //// DICE: debugging
                //printf("CIMarray CalculatePower Breakdown\n");
                //printf("%-20s %15.1f\n", "weightMatrixCol", weightMatrixCol);
                //printf("%-20s %15.1f\n", "numMuxing", numMuxing);
                //printf("%-20s %15.1f\n", "numBitInput", numBitInput);
                //printf("%-20s %15.4e pJ\n", "readDynamicEnergy", readDynamicEnergy*1e12);
                //printf("%-20s %15.4e pJ\n", "readDynamicEnergyArray", readDynamicEnergyArray*1e12);
                //printf("%-20s %15.4e pJ\n", "readDynamicEnergyOther", readDynamicEnergyOther*1e12);
                //printf("%-20s %15.4e pJ\n", "readDynamicEnergyADC", readDynamicEnergyADC*1e12);
                //printf("%-20s %15.4e pJ\n", "readDynamicEnergyAccum", readDynamicEnergyAccum*1e12);
                //printf("%-20s %15.4e pJ\n", "wlSwitchMatrix", wlSwitchMatrix.readDynamicEnergy*1e12);
                //printf("%-20s %15.4e pJ\n", "slSwitchMatrix", slSwitchMatrix.readDynamicEnergy*1e12); 
                //printf("%-20s %15.4e pJ\n", "mux", mux.readDynamicEnergy*1e12);
                //printf("%-20s %15.4e pJ\n", "muxDecoder", muxDecoder.readDynamicEnergy*1e12);
                //printf("%-20s %15.4e pJ\n", "multilevelSenseAmp", multilevelSenseAmp.readDynamicEnergy*1e12);
                //printf("%-20s %15.4e pJ\n", "multilevelSAEncoder", multilevelSAEncoder.readDynamicEnergy*1e12);
                //printf("%-20s %15.4e pJ\n", "sarADC", sarADC.readDynamicEnergy*1e12);
                //printf("%-20s %15.4e pJ\n", "shiftAddInput", shiftAddInput.readDynamicEnergy*1e12);
                //printf("%-20s %15.4e pJ\n", "shiftAddCell", shiftAddCell.readDynamicEnergy*1e12);

				
				// Leakage
				leakage = 0;
				leakage += wlSwitchMatrix.leakage;
				leakage += wlNewSwitchMatrix.leakage;
				leakage += slSwitchMatrix.leakage;
				leakage += ((numColMuxed > 1)==true? (mux.leakage):0);
				leakage += ((numColMuxed > 1)==true? (muxDecoder.leakage):0);
				leakage += multilevelSenseAmp.leakage;
				leakage += multilevelSAEncoder.leakage;
				leakage += shiftAddInput.leakage;
				leakage += shiftAddCell.leakage;
				
			} else if (BNNsequentialMode || XNORsequentialMode) {
                /* DICE: REMOVED */
			} else if (BNNparallelMode || XNORparallelMode) {
                /* DICE: REMOVED */
			}
		} 
	}
}

void CIMArray::PrintProperty() {

	if (cell.memCellType == Type::RRAM || cell.memCellType == Type::FeFET) {
		
		cout << endl << endl;
	    cout << "Array:" << endl;
	    cout << "Area = " << heightArray*1e6 << "um x " << widthArray*1e6 << "um = " << areaArray*1e12 << "um^2" << endl;
	    cout << "Read Dynamic Energy = " << readDynamicEnergyArray*1e12 << "pJ" << endl;
	    cout << "Write Dynamic Energy = " << writeDynamicEnergyArray*1e12 << "pJ" << endl;
		cout << "Write Latency = " << writeLatencyArray*1e9 << "ns" << endl;

		if (conventionalSequential) {
            /* DICE: REMOVED */
		} else if (conventionalParallel) {
			if (cell.accessType == CMOS_access) {
				wlNewSwitchMatrix.PrintProperty("wlNewSwitchMatrix");
			} else {
				wlSwitchMatrix.PrintProperty("wlSwitchMatrix");
			}
			slSwitchMatrix.PrintProperty("slSwitchMatrix");
			mux.PrintProperty("mux");
			muxDecoder.PrintProperty("muxDecoder");
			multilevelSenseAmp.PrintProperty("multilevelSenseAmp");
			multilevelSAEncoder.PrintProperty("multilevelSAEncoder");
			if (numReadPulse > 1) {
				shiftAddInput.PrintProperty("shiftAddInput");
			}
			if (numCellPerSynapse > 1) {
				shiftAddCell.PrintProperty("shiftAddCell");
			}
		} else if (BNNsequentialMode || XNORsequentialMode) {
            /* DICE: REMOVED */
		} else if (BNNparallelMode || XNORparallelMode) {
            /* DICE: REMOVED */
		} else {
            /* DICE: REMOVED */
		}
	} 
	FunctionUnit::PrintProperty("CIMArray");
	cout << "Used Area = " << usedArea*1e12 << "um^2" << endl;
	cout << "Empty Area = " << emptyArea*1e12 << "um^2" << endl;
}

