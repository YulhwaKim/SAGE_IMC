#include <iostream>
#include "constant.h"
#include "formula.h"
#include "DigitalElements.h"


DigitalElements::DigitalElements(const InputParameter& _inputParameter, const Technology& _tech, const vector<double> _designDE):
inputParameter(_inputParameter), tech(_tech), designDE(_designDE) {
    // place adderTree
    if (designDE[0] < 1) {
        placeAdderTree = false;
    } else {
        placeAdderTree = true;
        adderTree = new AdderTree(inputParameter, tech);
    }
    // place reLu
    if (designDE[3] < 1) {
        placeReLu = false;
    } else {
        placeReLu = true;
        reLu = new BitShifter(inputParameter, tech);
    }
    // place maxPooling 
    if (designDE[5] < 1) {
        placeMaxPooling = false;
    } else {
        placeMaxPooling = true;
        maxPooling = new MaxPooling(inputParameter, tech);
    }
    // not initialized
    initialized = false;
}

/* Initialize digital module */
void DigitalElements::Initialize(bool _fixedDataFlow, double _clkFreq) {
    // set dataflow option
    fixedDataFlow = _fixedDataFlow;
    // set clock frequency
    clkFreq = _clkFreq;
    // set digital circuit components
    if ( placeAdderTree ) {
        adderTree->Initialize((int)designDE[2]/*numAdd*/, (int)designDE[1]/*numBit*/, (int)designDE[0]/*numUnit*/, clkFreq);
    }
    if ( placeReLu ) {
        reLu->Initialize((int)designDE[3]/*numUnit*/, (int)designDE[4]/*numBit*/, clkFreq);
    }
    if ( placeMaxPooling ) {
        maxPooling->Initialize((int)designDE[6]/*numBit*/, (int)designDE[7]/*window*/, (int)designDE[5]/*numUnit*/, clkFreq);
    }
    // set initialized flag
    initialized = true;
}

/* Calculate area of digital module */
void DigitalElements::CalculateArea(double newWidth) {
    double area = 0;
    double areaAdderTree = 0;
    double areaReLu = 0;
    double areaMaxPooling = 0;

    areaVector.clear();

    if ( placeAdderTree ) {
        adderTree->CalculateArea(NULL, newWidth, NONE);
        areaAdderTree = adderTree->area;
    }
    if ( placeReLu ) {
        reLu->CalculateArea(NULL, newWidth, NONE);
        areaReLu = reLu->area;
    }
    if ( placeMaxPooling ) {
        maxPooling->CalculateUnitArea(NONE);
        maxPooling->CalculateArea(newWidth);
        areaMaxPooling = maxPooling->area;
    }

    area = areaAdderTree + areaReLu + areaMaxPooling;
    areaVector.push_back(area);
    areaVector.push_back(areaAdderTree);
    areaVector.push_back(areaReLu);
    areaVector.push_back(areaMaxPooling);

}

/* Calculate Latency of digital module */
void DigitalElements::CalculateLatency(const vector<double> infoReadDE) {
    double latency = 0;
    double latencyAdderTree = 0;
    double latencyReLu = 0;
    double latencyMaxPooling = 0;

    latencyVector.clear();

    if ( placeAdderTree ) {
        adderTree->CalculateLatency(infoReadDE[0]/*numRead*/, (int)infoReadDE[1]/*numUnitAdd*/, 0);
        latencyAdderTree = adderTree->readLatency;
    }
    if ( placeReLu ) {
        reLu->CalculateLatency(infoReadDE[2]/*numRead*/);
        latencyReLu = reLu->readLatency;
    }
    if ( placeMaxPooling ) {
        maxPooling->CalculateLatency(1e20, 0, infoReadDE[3]/*numRead*/);
        latencyMaxPooling = maxPooling->readLatency;
    }

    latency = latencyAdderTree + latencyReLu + latencyMaxPooling;
    latencyVector.push_back(latency);
    latencyVector.push_back(latencyAdderTree);
    latencyVector.push_back(latencyReLu);
    latencyVector.push_back(latencyMaxPooling);
}

/* Calculate readDynamicEnergy & leakage of digital module */
void DigitalElements::CalculatePower(const vector<double> infoReadDE) {
    double readDynamicEnergy = 0;
    double readDynamicEnergyAdderTree = 0;
    double readDynamicEnergyReLu = 0;
    double readDynamicEnergyMaxPooling = 0;

    readDynamicEnergyVector.clear();
    leakage = 0;

    if ( placeAdderTree ) {
        adderTree->CalculatePower(infoReadDE[0]/*numRead*/, infoReadDE[1]/*numUnitAdd*/);
        leakage += adderTree->leakage;
        readDynamicEnergyAdderTree = adderTree->readDynamicEnergy;
    }
    if ( placeReLu ) {
        reLu->CalculatePower(infoReadDE[2]/*numRead*/);
        leakage += reLu->leakage;
        readDynamicEnergyReLu = reLu->readDynamicEnergy;
    }
    if ( placeMaxPooling ) {
        maxPooling->CalculatePower(infoReadDE[3]/*numRead*/);
        leakage += maxPooling->leakage;
        readDynamicEnergyMaxPooling = maxPooling->readDynamicEnergy;
    }

    readDynamicEnergy = readDynamicEnergyAdderTree + readDynamicEnergyReLu + readDynamicEnergyMaxPooling;
    readDynamicEnergyVector.push_back(readDynamicEnergy);
    readDynamicEnergyVector.push_back(readDynamicEnergyAdderTree);
    readDynamicEnergyVector.push_back(readDynamicEnergyReLu);
    readDynamicEnergyVector.push_back(readDynamicEnergyMaxPooling); 
}

/* Calculate leakage of digital module */
void DigitalElements::CalculateLeakage() {
    leakage = 0;

    if ( placeAdderTree ) {
        adderTree->CalculatePower(0, adderTree->numSubcoreRow);
        leakage += adderTree->leakage;
    }
    if ( placeReLu ) {
        reLu->CalculatePower(0);
        leakage += reLu->leakage;
    }
    if ( placeMaxPooling ) {
        maxPooling->CalculatePower(0);
        leakage += maxPooling->leakage;
    }

}

void DigitalElements::PrintProperty(const char* str) {
    if ( placeAdderTree ) {
        adderTree->PrintProperty(str);
    }
    if ( placeReLu ) {
        reLu->PrintProperty(str);
    }
    if ( placeMaxPooling ) {
        maxPooling->PrintProperty(str);
    }
}





