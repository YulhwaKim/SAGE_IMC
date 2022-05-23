#include <cmath>
#include <iostream>
#include "constant.h"
#include "formula.h"
#include "HierarchyObject.h"
#include "Param.h"

extern Param *param;

HierarchyObject::HierarchyObject(
            InputParameter& _inputParameter, Technology& _tech, MemCell& _cell,
            const int _hlevel, const HierarchyRoot* _rootObject, const HierarchyObject* _subObject,
            const vector<double> _designObject):
inputParameter(_inputParameter), tech(_tech), cell(_cell), 
hlevel(_hlevel), rootObject(_rootObject), subObject(_subObject), 
designObject(_designObject){

    // get subObject placement design
    numSubObjectRow = (int)designObject[0];
    numSubObjectCol = (int)designObject[1];
    numSubObject = numSubObjectRow * numSubObjectCol;

    //// DICE: debugging
    //printf("numSubObjectRow: %d, numSubObjectCol: %d\n", numSubObjectRow, numSubObjectCol);

    // get DE, IC, BU design
    designDE.assign(designObject.begin() + 2 , designObject.begin() + 10);
    designIC.assign(designObject.begin() + 10, designObject.begin() + 15);
    designBU.assign(designObject.begin() + 15, designObject.begin() + 22);

    digitalElements = new DigitalElements(inputParameter, tech, designDE);
    interConnect = new InterConnect(inputParameter, tech, designIC);
    bufferUnit = new BufferUnit(inputParameter, tech, cell/*dummy*/, designBU);

    // get row/col count info
    if ( hlevel == 1 ) {
        numRowSubObject = rootObject->numRow;
        numColSubObject = rootObject->numCol;
        numOutBitSubObject = rootObject->numOutBit;
        // channel dim
        numInCSubObject = rootObject->numRow;
        numOutCSubObject = rootObject->numCol;
    } else {
        numRowSubObject = subObject->numRow;
        numColSubObject = subObject->numCol;
        numOutBitSubObject = subObject->numOutBit;
        // channel dim
        numInCSubObject = subObject->numInC;
        numOutCSubObject = subObject->numOutC;
    }
    numRow = numRowSubObject * numSubObjectRow;
    numCol = numColSubObject * numSubObjectCol;
    // channel dim
    if ( interConnect->outType == 1 ) { /* Linear Array (Operate as Systolic Array) */
        numOutC = numOutCSubObject;
    } else {
        numOutC = numOutCSubObject * numSubObjectCol;
    }
    if ( interConnect->inType == 1 ) { /* Linear Array (Operate as Systolic Array) */
        numInC = numInCSubObject;
    } else {
        numInC = numInCSubObject * numSubObjectRow;
    }

    numOutBit = (int)designDE[1] + 1; // numBit of adderTree + 1 (assume that merge does not increase the outbit more than 1)

    //// DICE: debugging
    //printf("desingDE\n");
    //for (int i = 0; i < designDE.size(); i++ ) {
    //    cout << designDE[i] << " ";
    //}
    //cout << endl;

    //printf("desingIC\n");
    //for (int i = 0; i < designIC.size(); i++ ) {
    //    cout << designIC[i] << " ";
    //}
    //cout << endl;

    //printf("desingBU\n");
    //for (int i = 0; i < designBU.size(); i++ ) {
    //    cout << designBU[i] << " ";
    //}
    //cout << endl;


    inputBuffer = (bufferUnit->inBUSize > 0)? true : false;

    // not initialized
    initialized = false ;

}

/* Initialize hierarchy object */
void HierarchyObject::Initialize(double _clkFreq) {

    clkFreq = _clkFreq;

    if ( hlevel < 1 ) {
        cerr << "hlevel of HiearchyObject should not be smaller than 1, but we got " << hlevel << endl;
    } else if ( hlevel == 1) { /* sub-object is rootObject */
        subWidth = rootObject->width;
        subHeight = rootObject->height;
    } else {
        subWidth = subObject->width;
        subHeight = subObject->height;
    }

    digitalElements->Initialize((designBU[4]/*inBUSize*/==0)/*fixedDataFlow*/, clkFreq);
    bufferUnit->Initialize(param->unitLengthWireResistance, clkFreq);
    interConnect->Initialize(numSubObjectRow, numSubObjectCol, subHeight, subWidth, bufferUnit->inBUSize, clkFreq);


    // calculate area as the area information might be used on the childern object
    CalculateArea();
    // calculate leakage info
    CalculateLeakage();

    // set initialization flag
    initialized = true;

}


/* Calculate area of hierarchy object */
/* Asssume that the area of subObject is already calculated */
void HierarchyObject::CalculateArea() {

    // obtain subobject area/width/height info
    vector<double> subAreaVector;
    double subArea;
    if ( hlevel < 1 ) {
        cerr << "hlevel of HiearchyObject should not be smaller than 1, but we got " << hlevel << endl;
    } else if ( hlevel == 1) { /* sub-object is rootObject */
        subAreaVector = rootObject->areaVector;
        subArea = rootObject->usedArea;
    } else {
        subAreaVector = subObject->areaVector;
        subArea = subObject->area;
    }
    
    // obtain area information of the elements
    digitalElements->CalculateArea(numSubObjectCol * subWidth);
    interConnect->CalculateArea();
    bufferUnit->CalculateArea(numSubObjectRow * subHeight, numSubObjectCol * subWidth);

    area = subArea * numSubObject +
        digitalElements->areaVector[0] + 
        bufferUnit->area +
        interConnect->area;

    height = sqrt(area);
    width = area/height;

    // clear vector before update
    areaVector.clear();
    areaVector2.clear();

    //breakdown type1
    areaVector.push_back(area); // total
    areaVector.push_back(subAreaVector[1] * numSubObject); // array
    areaVector.push_back(subAreaVector[2] * numSubObject); // ADC
    areaVector.push_back(subAreaVector[3] * numSubObject + digitalElements->areaVector[1]); // accum
    areaVector.push_back(subAreaVector[4] * numSubObject + bufferUnit->area); // buffer
    areaVector.push_back(subAreaVector[5] * numSubObject + interConnect->area); // ic
    areaVector.push_back(subAreaVector[6] * numSubObject + 
                        digitalElements->areaVector[2] + digitalElements->areaVector[3]); // other digital

    // breakdown type2
    areaVector2.push_back(area); // total
    areaVector2.push_back(subArea * numSubObject); // subObject
    areaVector2.push_back(digitalElements->areaVector[1]); // accum
    areaVector2.push_back(bufferUnit->area); // buffer
    areaVector2.push_back(interConnect->area); // ic
    areaVector2.push_back(digitalElements->areaVector[2] + digitalElements->areaVector[3]); // other digital (reLu, Max)

}


/* Calculate Latency of hierarchy object */
/* Asssume that the latency of subObject is already calculated */
void HierarchyObject::CalculateLatency(const vector<double> infoRead, const vector<double> subLatencyVector) {

    // get infoRead
    vector<double> infoReadDE {&infoRead[0],  &infoRead[4]};
    vector<double> infoReadIC {&infoRead[4],  &infoRead[10]};
    vector<double> infoReadBU {&infoRead[10], &infoRead[16]};

    // sub latency information (sync, cycle counting)
    digitalElements->CalculateLatency(infoReadDE);
    interConnect->CalculateLatency(infoReadIC);
    bufferUnit->CalculateLatency(infoReadBU);

    double latency = subLatencyVector[0] +
                    digitalElements->latencyVector[0] +
                    bufferUnit->latency +
                    interConnect->readLatency;

    // clear vector before update
    latencyVector.clear();
    latencyVector2.clear();

    // breakdown type1
    latencyVector.push_back(latency); // total
    latencyVector.push_back(subLatencyVector[1]); // array
    latencyVector.push_back(subLatencyVector[2] + digitalElements->latencyVector[1]); // accum
    latencyVector.push_back(subLatencyVector[3] + bufferUnit->latency); // buffer
    latencyVector.push_back(subLatencyVector[4] + interConnect->readLatency); // ic
    latencyVector.push_back(subLatencyVector[5] + 
                            digitalElements->latencyVector[2] + digitalElements->latencyVector[3]); // other digital

    // breakdown type2
    latencyVector2.push_back(latency); // total
    latencyVector2.push_back(subLatencyVector[0]); // subObject
    latencyVector2.push_back(digitalElements->latencyVector[1]); // accum
    latencyVector2.push_back(bufferUnit->latency); // buffer
    latencyVector2.push_back(interConnect->readLatency); //ic
    latencyVector2.push_back(digitalElements->latencyVector[2] + digitalElements->latencyVector[3]); // other digital

    //// DICE: debugging
    //printf("%d-th Object Latency breakdown (cycle)\n", hlevel);
    //printf("%-20s %15.1f\n", "Object", latencyVector[0]);
    //printf("%-20s %15.1f\n", "subArray", latencyVector[1]);
    //printf("%-20s %15.1f\n", "accum", latencyVector[2]);
    //printf("%-20s %15.1f\n", "buffer", latencyVector[3]);
    //printf("%-20s %15.1f\n", "ic", latencyVector[4]);
    //printf("%-20s %15.1f\n", "digit", latencyVector[5]);

}

/* Calculate readDynamicEnergy of hierarcy object */
/* Asssume that the energy of subObject is already calculated */
void HierarchyObject::CalculatePower(const vector<double> infoRead,  
                                    const vector<double> subReadDynamicEnergyVector) {

    // get infoRead
    vector<double> infoReadDE {&infoRead[0],  &infoRead[4]};
    vector<double> infoReadIC {&infoRead[4],  &infoRead[10]};
    vector<double> infoReadBU {&infoRead[10], &infoRead[16]};
   
    // sub power information
    digitalElements->CalculatePower(infoReadDE);
    interConnect->CalculatePower(infoReadIC);
    bufferUnit->CalculatePower(infoReadBU);

    double readDynamicEnergy = subReadDynamicEnergyVector[0] +
                            digitalElements->readDynamicEnergyVector[0] +
                            bufferUnit->dynamicEnergy +
                            interConnect->readDynamicEnergy;

    // clear vector before update
    readDynamicEnergyVector.clear();
    readDynamicEnergyVector2.clear();

    //breakdown type1
    readDynamicEnergyVector.push_back(readDynamicEnergy); // total
    readDynamicEnergyVector.push_back(subReadDynamicEnergyVector[1]); // array
    readDynamicEnergyVector.push_back(subReadDynamicEnergyVector[2]); // ADC
    readDynamicEnergyVector.push_back(subReadDynamicEnergyVector[3] + 
                                                digitalElements->readDynamicEnergyVector[1]); // accum
    readDynamicEnergyVector.push_back(subReadDynamicEnergyVector[4] + bufferUnit->dynamicEnergy); // buffer
    readDynamicEnergyVector.push_back(subReadDynamicEnergyVector[5] + interConnect->readDynamicEnergy); // ic
    readDynamicEnergyVector.push_back(subReadDynamicEnergyVector[6] + 
                                                digitalElements->readDynamicEnergyVector[2] +
                                                digitalElements->readDynamicEnergyVector[3]); // other digital

    //breakdown type2
    readDynamicEnergyVector2.push_back(readDynamicEnergy); // total
    readDynamicEnergyVector2.push_back(subReadDynamicEnergyVector[0]); // subObject
    readDynamicEnergyVector2.push_back(digitalElements->readDynamicEnergyVector[1]); // accum
    readDynamicEnergyVector2.push_back(bufferUnit->dynamicEnergy); // buffer
    readDynamicEnergyVector2.push_back(interConnect->readDynamicEnergy); // ic
    readDynamicEnergyVector2.push_back(digitalElements->readDynamicEnergyVector[2] +
                                    digitalElements->readDynamicEnergyVector[3]); // other digital

    ////// DICE: debugging
    //if ( hlevel == 2 ) {
    //    printf("%-20s %2d\n", "accum numRead", (int)infoReadDE[0]);
    //    printf("%-20s %2d\n", "accum numUnitAdd", (int)infoReadDE[1]);
    //    printf("%-20s %15.4e pJ\n", "accum", readDynamicEnergyVector2[2]*1e12);
    //    //printf("%-20s %15.4e pJ\n", "sub accum", subReadDynamicEnergyVector[3]*1e12);
    //}
    //if ( hlevel == 3 ) {
    //    printf("%-20s %15.4e pJ\n", "3-sub accum", subReadDynamicEnergyVector[3]*1e12);
    //}
    //printf("%d-th Object readDynamicEnergy breakdown\n", hlevel);
    //printf("%-20s %15.4e pJ\n", "Object", readDynamicEnergyVector[0]*1e12);
    //printf("%-20s %15.4e pJ\n", "subArray", readDynamicEnergyVector[1]*1e12);
    //printf("%-20s %15.4e pJ\n", "ADC", readDynamicEnergyVector[2]*1e12);
    //printf("%-20s %15.4e pJ\n", "accum", readDynamicEnergyVector[3]*1e12);
    //printf("%-20s %15.4e pJ\n", "buffer", readDynamicEnergyVector[4]*1e12);
    //printf("%-20s %15.4e pJ\n", "ic", readDynamicEnergyVector[5]*1e12);
    //printf("%-20s %15.4e pJ\n", "other", readDynamicEnergyVector[6]*1e12);
    
}

/* Calculate Latency of hierarchy object IC */
void HierarchyObject::CalculateICLatency(const vector<double> infoReadIC) {
    interConnect->CalculateLatency(infoReadIC);
}

/* Calculate readDynamicEnergy of hierarcy object IC*/
void HierarchyObject::CalculateICPower(const vector<double> infoReadIC) {
    interConnect->CalculatePower(infoReadIC);
}

/* Calculate leakage of hierarcy object */
/* Asssume that leakage of subObject is already calculated */
void HierarchyObject::CalculateLeakage() {

    leakage = 0;

    double subLeakage;
    // get subObject leakage
    if ( hlevel < 1 ) {
        cerr << "hlevel of HiearchyObject should not be smaller than 1, but we got " << hlevel << endl;
    } else if ( hlevel == 1) { /* sub-object is rootObject */
        subLeakage = rootObject->leakage * numSubObject;
    } else {
        subLeakage = subObject->leakage * numSubObject;
    }

    // sub leakage information
    digitalElements->CalculateLeakage();
    interConnect->CalculateLeakage();
    bufferUnit->CalculateLeakage();

    leakage = subLeakage + 
            digitalElements->leakage +
            bufferUnit->leakage +
            interConnect->leakage;


}

//void HierarchyObject::PrintProperty(const char* str) {
//}




