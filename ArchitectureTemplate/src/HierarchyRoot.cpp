#include <iostream>
#include <cmath>
#include "constant.h"
#include "formula.h"
#include "HierarchyRoot.h"
#include "Param.h"

extern Param *param;

HierarchyRoot::HierarchyRoot(InputParameter& _inputParameter, Technology& _tech, MemCell& _cell):
inputParameter(_inputParameter), tech(_tech), cell(_cell) {
	/*** circuit level parameters ***/
	switch(param->memcelltype) {
		case 2:	    cell.memCellType = Type::RRAM; break;
		case 1:	    cell.memCellType = Type::SRAM; break;
		case -1:	break;
		default:	exit(-1);
	}
	switch(param->accesstype) {
		case 1:	    cell.accessType = CMOS_access;  break;
		case -1:	break;
		default:	exit(-1);
	}				
					
	switch(param->transistortype) {
		case 3:	    inputParameter.transistorType = TFET;          break;
		case 2:	    inputParameter.transistorType = FET_2D;        break;
		case 1:	    inputParameter.transistorType = conventional;  break;
		case -1:	break;
		default:	exit(-1);
	}
	
	switch(param->deviceroadmap) {
		case 2:	    inputParameter.deviceRoadmap = LSTP;  break;
		case 1:	    inputParameter.deviceRoadmap = HP;    break;
		case -1:	break;
		default:	exit(-1);
	}
	inputParameter.temperature = param->temp;   // Temperature (K)
	inputParameter.processNode = param->technode;    // Technology node
	tech.Initialize(inputParameter.processNode, inputParameter.deviceRoadmap, inputParameter.transistorType);
	
	cell.resistanceOn = param->resistanceOn;	                                // Ron resistance at Vr in the reported measurement data (need to recalculate below if considering the nonlinearity)
	cell.resistanceOff = param->resistanceOff;	                                // Roff resistance at Vr in the reported measurement dat (need to recalculate below if considering the nonlinearity)
	cell.resistanceAvg = (cell.resistanceOn + cell.resistanceOff)/2;            // Average resistance (for energy estimation)
	cell.readVoltage = param->readVoltage;	                                    // On-chip read voltage for memory cell
	cell.readPulseWidth = param->readPulseWidth;
	cell.accessVoltage = param->accessVoltage;                                       // Gate voltage for the transistor in 1T1R
	cell.resistanceAccess = param->resistanceAccess;
	cell.featureSize = param->featuresize; 
	cell.writeVoltage = param->writeVoltage;

	if (cell.memCellType == Type::SRAM) {   // SRAM
		cell.heightInFeatureSize = param->heightInFeatureSizeSRAM;                   // Cell height in feature size
		cell.widthInFeatureSize = param->widthInFeatureSizeSRAM;                     // Cell width in feature size
		cell.widthSRAMCellNMOS = param->widthSRAMCellNMOS;
		cell.widthSRAMCellPMOS = param->widthSRAMCellPMOS;
		cell.widthAccessCMOS = param->widthAccessCMOS;
		cell.minSenseVoltage = param->minSenseVoltage;
	} else {
		cell.heightInFeatureSize = (cell.accessType==CMOS_access)? param->heightInFeatureSize1T1R : param->heightInFeatureSizeCrossbar;         // Cell height in feature size
		cell.widthInFeatureSize = (cell.accessType==CMOS_access)? param->widthInFeatureSize1T1R : param->widthInFeatureSizeCrossbar;            // Cell width in feature size
	} 

	cimArray = new CIMArray(inputParameter, tech, cell);
		
	/* Create CIMArray object initialization */
	cimArray->XNORparallelMode = param->XNORparallelMode;               
	cimArray->XNORsequentialMode = param->XNORsequentialMode;             
	cimArray->BNNparallelMode = param->BNNparallelMode;                
	cimArray->BNNsequentialMode = param->BNNsequentialMode;              
	cimArray->conventionalParallel = param->conventionalParallel;                  
	cimArray->conventionalSequential = param->conventionalSequential;                 
	cimArray->numRow = param->numRowCIMArray;
	cimArray->numCol = param->numRowCIMArray;
	cimArray->levelOutput = param->levelOutput;
	cimArray->numColMuxed = param->numColMuxed;               // How many columns share 1 read circuit (for neuro mode with analog RRAM) or 1 S/A (for memory mode or neuro mode with digital RRAM)
    cimArray->clkFreq = param->clkFreq;                       // Clock frequency
	cimArray->relaxArrayCellHeight = param->relaxArrayCellHeight;
	cimArray->relaxArrayCellWidth = param->relaxArrayCellWidth;
	cimArray->numReadPulse = param->numBitInput;
	cimArray->avgWeightBit = param->cellBit;
	cimArray->numCellPerSynapse = param->numColPerSynapse;
	cimArray->SARADC = param->SARADC;
	cimArray->currentMode = param->currentMode;
	cimArray->validated = param->validated;
	cimArray->spikingMode = NONSPIKING;
	
	numRow = param->numRowCIMArray;
	numCol = param->numColCIMArray;
	
	if (cimArray->numColMuxed > numCol) {                      // Set the upperbound of numColMuxed
		cimArray->numColMuxed = numCol;
	}

	cimArray->numReadCellPerOperationFPGA = numCol;	           // Not relevant for IMEC
	cimArray->numWriteCellPerOperationFPGA = numCol;	       // Not relevant for IMEC
	cimArray->numReadCellPerOperationMemory = numCol;          // Define # of SRAM read cells in memory mode because SRAM does not have S/A sharing (not relevant for IMEC)
	cimArray->numWriteCellPerOperationMemory = numCol/8;       // # of write cells per operation in SRAM memory or the memory mode of multifunctional memory (not relevant for IMEC)
	cimArray->numWriteCellPerOperationNeuro = numCol;	       // For SRAM or analog RRAM in neuro mode
    cimArray->maxNumWritePulse = MAX(cell.maxNumLevelLTP, cell.maxNumLevelLTD);

    maxConductance = param->maxConductance;                    // max conductance level of the memory cell
    minConductance = param->minConductance;                    // min conductance level of the memroy cell

    //** NOTE: we should support bitplane-wise ratio in the future **//
    inputActiveRatio = param->inputActiveRatio;                 // ratio of active value in each bitplane of input
    weightLevelRatioVector = param->weightLevelRatioVector;     // ratio of each weight-cell level

	/*** initialize modules ***/
	cimArray->Initialize(numRow, numCol, param->unitLengthWireResistance);        // initialize cimArray

    numOutBit = cimArray->numOutBit;
    // calculate area as the area information might be used on the childern object
	CalculateArea();
    // calculate Leakage
    CalculateLeakage();

    // set initialized flag
    initialized = false;
}

/* Initialize hierarchy root - CIM array */
void HierarchyRoot::Initialize() {

    GetColumnResistance();
    // set initialized flag
    initialized = true;

}

void HierarchyRoot::GetColumnResistance() {

    int cellRange = pow(2, param->cellBit);
    double cellG, cellR;
    double columnG = 0;

    // get conductance values of memory cells for each level
    if (cell.memCellType == Type::RRAM) {
        if (cell.accessType == CMOS_access) {
            for (int i=0; i < cellRange; i++ ) {
                cellG = i / (cellRange - 1) * (maxConductance - minConductance) + minConductance;
                cellR = (double) 1.0 / cellG + cell.resistanceAccess;
                columnG += (double) 1.0 / cellR * (numRow * inputActiveRatio * weightLevelRatioVector[i]);
            }
        } else{
            cerr << "[Error] We support CMOS_access for RRAM only." << endl;
        }
    } else if (cell.memCellType == Type::SRAM) {
                cellR = (double) (cimArray->resCellAccess + param->wireResistanceCol);
                columnG += (double) 1.0 / cellR * numRow * inputActiveRatio;
    } else {
        cerr << "[Error] We support RRAM & SRAM memory cell only, not memCellType: " << cell.memCellType << endl;
    } 
   
    columnRes = 1.0 / columnG;

}

/* Calculate area of hierarchy root - CIM array */
void HierarchyRoot::CalculateArea() {

	cimArray->CalculateArea();

    usedArea = cimArray->usedArea;
	area = cimArray->area;
    height = cimArray->height;
    width = cimArray->width;
	
    // clear vector before update
    areaVector.clear();

	areaVector.push_back(area);
    // breakdown
	areaVector.push_back(cimArray->areaArray + cimArray->areaOther); // array
	areaVector.push_back(cimArray->areaADC);  // ADC
	areaVector.push_back(cimArray->areaAccum); // accum
    areaVector.push_back(0); // buffer
    areaVector.push_back(0); // ic
    areaVector.push_back(0); // other digital

}

/* Calculate Latency of hierarchy root - CIM array */
void HierarchyRoot::CalculateLatency(vector<double> infoReadCIM) {
    double weightMatrixRow = infoReadCIM[0];
    double weightMatrixCol = infoReadCIM[1];
    double numBitInput = infoReadCIM[2];
    double numCellPerSynapse = infoReadCIM[3];
    double numRead = ceil(infoReadCIM[4]);

    double dummy;

    CalculateLatency(false, weightMatrixRow, weightMatrixCol,
                    numBitInput, numCellPerSynapse,
                    numRead, &dummy);

}

/* Calculate Latency of hierarchy root - CIM array */
void HierarchyRoot::CalculateLatency(bool CalculateclkFreq, double weightMatrixRow, double weightMatrixCol, 
                                    double numBitInput, double numCellPerSynapse,
                                    double numRead, double *clkPeriod) {

    double colR = columnRes * ( weightMatrixRow / numRow ); 

    cimArray->CalculateLatency(colR, weightMatrixRow, weightMatrixCol, numBitInput, CalculateclkFreq);

    if(CalculateclkFreq && (*clkPeriod < cimArray->readLatency)){
    	*clkPeriod = cimArray->readLatency;					//clk freq is decided by the longest sensing latency
    }							
    
    if(!CalculateclkFreq){
        // clear vector before update
        latencyVector.clear();

        numRead = ceil(numRead); // for latency --> ceil!!

    	latencyVector.push_back(cimArray->readLatency * numRead);		
        // breakdown
    	latencyVector.push_back( (cimArray->readLatencyOther + cimArray->readLatencyADC) * numRead ); // array (include ADC) 
    	latencyVector.push_back(cimArray->readLatencyAccum * numRead);		// accum
        latencyVector.push_back(0); // buffer
        latencyVector.push_back(0); // ic
        latencyVector.push_back(0); // other digital
    }

}

/* Calculate Power of hierarchy root - CIM array */
void HierarchyRoot::CalculatePower(vector<double> infoReadCIM) {
    double weightMatrixRow = infoReadCIM[0];
    double weightMatrixCol = infoReadCIM[1];
    double numBitInput = infoReadCIM[2];
    double numCellPerSynapse = infoReadCIM[3];
    double numRead = ceil(infoReadCIM[4]);

    CalculatePower(weightMatrixRow, weightMatrixCol,
                    numBitInput, numCellPerSynapse,
                    numRead);

}

/* Calculate readDynamicEnergy & leakage of hierarcy root - CIM array */
void HierarchyRoot::CalculatePower(double weightMatrixRow, double weightMatrixCol, 
                                double numBitInput, double numCellPerSynapse, double numRead) {

    double colR = columnRes * ( weightMatrixRow / numRow ); 
    cimArray->CalculatePower(colR, weightMatrixRow, weightMatrixCol, numBitInput, inputActiveRatio);

    // leakage
    leakage = cimArray->leakage;

    // clear vector before update
    readDynamicEnergyVector.clear();

    // dynamic energy
    readDynamicEnergyVector.push_back(cimArray->readDynamicEnergy * numRead);
    // breakdown of dynamic energy
    readDynamicEnergyVector.push_back( (cimArray->readDynamicEnergyArray + cimArray->readDynamicEnergyOther) * numRead); // array
    readDynamicEnergyVector.push_back(cimArray->readDynamicEnergyADC * numRead); // ADC
    readDynamicEnergyVector.push_back(cimArray->readDynamicEnergyAccum * numRead); // accum
    readDynamicEnergyVector.push_back(0); // buffer
    readDynamicEnergyVector.push_back(0); // ic
    readDynamicEnergyVector.push_back(0); // other

}

/* Calculate leakage of hierarcy root - CIM array */
void HierarchyRoot::CalculateLeakage() {

    cimArray->CalculatePower(columnRes, (double)cimArray->numRow, (double)cimArray->numCol, (double)cimArray->numReadPulse, 1);

    // leakage
    leakage = cimArray->leakage;

}

void HierarchyRoot::PrintProperty() {
    cimArray->PrintProperty();
}




