#include <iostream>
#include "constant.h"
#include "formula.h"
#include "InterConnect.h"

InterConnect::InterConnect(const InputParameter& _inputParameter, const Technology& _tech, const vector<double> _designIC):
inputParameter(_inputParameter), tech(_tech), designIC(_designIC) {

    delaytolerance = designIC[0];
    outType = (int)designIC[1];
    outBusWidth = designIC[2];
    inType = (int)designIC[3];
    inBusWidth = designIC[4];


    // outBus
    if ( outType == 0 /* Bus */) { 
        outBus = new Bus(inputParameter, tech);    
    } else if ( outType == 1 /* LinearArray */) {
        outLinear = new LinearArray(inputParameter, tech);
    } else if ( outType == 2 /* 2D Mesh */) {
        mesh = new Mesh(inputParameter, tech);
        inType = 2; // when outType is 2D Mesh, set inType as 2D Mesh, too.
        flitSize =  (int)outBusWidth; 
        numPort = (int)designIC[3];
        inBusWidth = outBusWidth; // Rx
    } else if ( outType == 3 /* hierarchical Bus */ ) {
        outHBus = new HBus(inputParameter, tech);
    } else {
        cerr << "[InterConnect] Error: outType should be one of [0 (Bus), 1 (Linear), 2 (mesh), 3 (HBus)]!" << endl;
        exit(-1);
    }

    // inBus 
    if ( inType == 0 /* Bus */) { 
        inBus = new Bus(inputParameter, tech);    
    } else if ( inType == 1 /* LinearArray */) {
        inLinear = new LinearArray(inputParameter, tech);
    } else if ( inType == 2 /* 2D Mesh */) {
        /* NOTE: in the 2D Mesh case, input & output share a single 'mesh' for communication*/
    } else if ( inType == 3 /* hierarchical Bus */) {
        inHBus  = new HBus(inputParameter, tech);
    } else {
        cerr << "[InterConnect] Error: inType should be one of [0 (Bus), 1 (Linear), 2 (mesh), 3 (HBus)]!" << endl;
        exit(-1);
    }
    inBusMode = HORIZONTAL;
    // not initialized
    initialized = false;

}

/* Initialize interconnect */
void InterConnect::Initialize(int numRow, int numCol, double _unitHeight, double _unitWidth, int inBUSize, double _clkFreq) {
    // set properties
    clkFreq = _clkFreq;
    unitHeight = _unitHeight;
    unitWidth = _unitWidth;    

    // outBus
    if ( outType == 0 /* Bus */) { 
        outBus->Initialize(VERTICAL, false, numRow, numCol, delaytolerance, outBusWidth, unitHeight, unitWidth, 1, clkFreq);
    } else if ( outType == 1 /* LinearArray */) {
        outLinear->Initialize(VERTICAL, false, numRow, numCol, delaytolerance, outBusWidth, unitHeight, unitWidth, 1, clkFreq);
    } else if ( outType == 2 /* 2D Mesh */) {
        mesh->Initialize(numPort, flitSize, numRow, numCol, delaytolerance, unitHeight, unitWidth, 1, clkFreq);
    } else if ( outType == 3 /* hiearchical Bus */) {
        outHBus->Initialize(VERTICAL, false, numRow, numCol, delaytolerance, outBusWidth, unitHeight, unitWidth, 1, clkFreq);
    }

    // inBus
    if ( inType == 0 /* Bus */) { 
        inBusMode = (inBUSize > 0)? HORIZONTAL : VERTICAL;
        inBus->Initialize(inBusMode, false, numRow, numCol, delaytolerance, inBusWidth, unitHeight, unitWidth, 1, clkFreq);
    } else if ( inType == 1 /* LinearArray */) {
        inBusMode = HORIZONTAL;
        bool inputWire = (inBUSize > 0)? false : true;
        inLinear->Initialize(inBusMode, inputWire, numRow, numCol, delaytolerance, inBusWidth, 
                            unitHeight, unitWidth, 1, clkFreq);
    } else if ( inType == 2 /* 2D Mesh */) {
        /* NOTE: in the 2D Mesh case, input & output share a single 'mesh' for communication*/
        inBusMode = HORIZONTAL;
    } else if ( inType == 3 /* hierarchical Bus */) {
        inBusMode = HORIZONTAL;
        bool inputWire = (inBUSize > 0)? false : true;
        //bool inputWire = false;
        inHBus->Initialize(inBusMode, inputWire, numRow, numCol, delaytolerance, inBusWidth, unitHeight, unitWidth, 1, clkFreq);
        //inBusMode = (inBUSize > 0)? HORIZONTAL : VERTICAL;
        //inHBus->Initialize(inBusMode, false, numRow, numCol, delaytolerance, inBusWidth, unitHeight, unitWidth, 1, clkFreq);
    }

    initialized = true;
}

/* Calculate area of interconnect */
void InterConnect::CalculateArea() {

    area = 0;

    // outBus
    if ( outType == 0 /* Bus */) { 
        outBus->CalculateArea();
        area += outBus->area;
    } else if ( outType == 1 /* LinearArray */) {
        outLinear->CalculateArea();
        area += outLinear->area;
    } else if ( outType == 2 /* 2D Mesh */) {
        mesh->CalculateArea();
        area += mesh->area;
    } else if ( outType == 3 /* hierarchical Bus */) {
        outHBus->CalculateArea();
        area += outHBus->area;
    }

    // inBus
    if ( inType == 0 /* Bus */) { 
        inBus->CalculateArea();
        area += inBus->area;
    } else if ( inType == 1 /* LinearArray */) {
        inLinear->CalculateArea();
        area += inLinear->area;
    } else if ( inType == 2 /* 2D Mesh */) {
        /* NOTE: in the 2D Mesh case, input & output share a single 'mesh' for communication*/
    } else if ( inType == 3 /* hierarchical Bus */) {
        inHBus->CalculateArea();
        area += inHBus->area;
    }
}

/* Calculate Latency of interconnect */
void InterConnect::CalculateLatency(const vector<double> infoReadIC) {

    int state = (int)infoReadIC[0];
    int dataType = (int)infoReadIC[1];

    if ( state == 0 ) { return; }  // state 0 - pass
    else if ( state == 1 ) { readLatency = 0; } // state 1 - reset, else - cumulate

    if ( dataType == 0 /*output*/) {
        if ( outType == 0 /* Bus */) { 
            outBus->CalculateLatency(infoReadIC[2]/*numRead*/);
            readLatency += outBus->readLatency;
        } else if ( outType == 1 /* LinearArray */) {
            outLinear->CalculateLatency(infoReadIC[3]/*numActiveRow*/, 
                                        infoReadIC[4]/*numActiveCol*/, infoReadIC[2]/*numRead*/);
            readLatency += outLinear->readLatency;
        } else if ( outType == 2 /* 2D Mesh */) {
            mesh->CalculateLatency(infoReadIC[3]/*numHopsRow*/,
                                   infoReadIC[4]/*numHopsCol*/, infoReadIC[2]/*numRead*/);
            readLatency += mesh->readLatency;
        } else if ( outType == 3 /* hierarchical Bus */) {
            outHBus->CalculateLatency(infoReadIC[3]/*numHops*/, infoReadIC[4]/*numBusAccess*/,
                                    infoReadIC[2]/*numRead*/);
            readLatency += outHBus->readLatency;
        }
    } else { /*input*/
        if ( inType == 0 /* Bus */) { 
            inBus->CalculateLatency(infoReadIC[2]/*numRead*/);
            readLatency += inBus->readLatency;
        } else if ( inType == 1 /* LinearArray */) {
            inLinear->CalculateLatency(infoReadIC[3]/*numActiveRow*/, 
                                       infoReadIC[4]/*numActiveCol*/, infoReadIC[2]/*numRead*/);
            readLatency += inLinear->readLatency;
        } else if ( inType == 2 /* 2D Mesh */) {
            mesh->CalculateLatency(infoReadIC[3]/*numHopsRow*/,
                                   infoReadIC[4]/*numHopsCol*/, infoReadIC[2]/*numRead*/);
            readLatency += mesh->readLatency;
        } else if ( inType == 3 /* hierarchical Bus */) {
            inHBus->CalculateLatency(infoReadIC[3]/*numHops*/, infoReadIC[4]/*numBusAccess*/,
                                    infoReadIC[2]/*numRead*/);
            readLatency += inHBus->readLatency;
        }
    }

}

/* Calculate readDynamicEnergy of interconnect */
void InterConnect::CalculatePower(const vector<double> infoReadIC) {

    int state = (int)infoReadIC[0];
    int dataType = (int)infoReadIC[1];
    
    if ( state == 0 ) { return; }  // state 0 - pass
    else if ( state == 1 ) { readDynamicEnergy = 0; } // state 1 - reset, else - cumulate

    //printf("%10.4f, %10.4f, %10.4f\n", infoReadIC[3], infoReadIC[4], infoReadIC[2]);

    if ( dataType == 0 /*output*/) {
        if ( outType == 0 /* Bus */) { 
            outBus->CalculatePower(outBus->busWidth/*numBitAccess*/, infoReadIC[2]/*numRead*/);
            readDynamicEnergy += outBus->readDynamicEnergy;
        } else if ( outType == 1 /* LinearArray */) {
            outLinear->CalculatePower(infoReadIC[3]/*numActiveRow*/, infoReadIC[4]/*numActiveCol*/, 
                                    outLinear->busWidth/*numBitAccess*/, infoReadIC[2]/*numRead*/);
            readDynamicEnergy += outLinear->readDynamicEnergy;
        } else if ( outType == 2 /* 2D Mesh */) {
            mesh->CalculatePower(infoReadIC[3]/*numHopsRow*/, infoReadIC[4]/*numHopsCol*/, 
                                mesh->flitSize/*numBitAccess*/, infoReadIC[2]/*numRead*/);
            readDynamicEnergy += mesh->readDynamicEnergy;
        } else if ( outType == 3 /* hierarchical Bus */) {
            outHBus->CalculatePower(infoReadIC[3]/*numHops*/, infoReadIC[4]/*numBusAccess*/,
                                    outHBus->busWidth/*numBitAccess*/, infoReadIC[2]/*numRead*/);
            readDynamicEnergy += outHBus->readDynamicEnergy;
        }
    } else { /*input*/
        if ( inType == 0 /* Bus */) { 
            inBus->CalculatePower(inBus->busWidth /*numBitAccess*/, infoReadIC[2]/*numRead*/);
            readDynamicEnergy += inBus->readDynamicEnergy;
        } else if ( inType == 1 /* LinearArray */) {
            inLinear->CalculatePower(infoReadIC[3]/*numActiveRow*/, infoReadIC[4]/*numActiveCol*/, 
                                    inLinear->busWidth/*numBitAccess*/, infoReadIC[2]/*numRead*/);
            readDynamicEnergy += inLinear->readDynamicEnergy;
        } else if ( inType == 2 /* 2D Mesh */) {
            mesh->CalculatePower(infoReadIC[3]/*numHopsRow*/, infoReadIC[4]/*numHopsCol*/, 
                                mesh->flitSize/*numBitAccess*/, infoReadIC[2]/*numRead*/);
            readDynamicEnergy += mesh->readDynamicEnergy;
        } else if ( inType == 3 /* hierarchical Bus */) {
            inHBus->CalculatePower(infoReadIC[3]/*numHops*/, infoReadIC[4]/*numBusAccess*/,
                                    inHBus->busWidth/*numBitAccess*/, infoReadIC[2]/*numRead*/);
            readDynamicEnergy += inHBus->readDynamicEnergy;
        }
    }

}

/* Calculate leakage of interconnect */
void InterConnect::CalculateLeakage() {

    leakage = 0;

    // outBus
    if ( outType == 0 /* Bus */) { 
        outBus->CalculatePower(0, 0);
        leakage += outBus->leakage;
    } else if ( outType == 1 /* LinearArray */) {
        outLinear->CalculatePower(0, 0, 0, 0);
        leakage += outLinear->leakage;
    } else if ( outType == 2 /* 2D Mesh */) {
        mesh->CalculatePower(0, 0, 0, 0);
        leakage += mesh->leakage;
    } else if ( outType == 3 /* hierarchical Bus */) {
        outHBus->CalculatePower(0, 0, 0, 0);
        leakage += outHBus->leakage;
    }

    // inBus
    if ( inType == 0 /* Bus */) { 
        inBus->CalculatePower(0, 0);
        leakage += inBus->leakage;
    } else if ( inType == 1 /* LinearArray */) {
        inLinear->CalculatePower(0, 0, 0, 0);
        leakage += inLinear->leakage;
    } else if ( inType == 2 /* 2D Mesh */) {
        /* NOTE: in the 2D Mesh case, input & output share a single 'mesh' for communication*/
    } else if ( inType == 3 /* hierarchical Bus */) {
        inHBus->CalculatePower(0, 0, 0, 0);
        leakage += inHBus->leakage;
    }

}

void InterConnect::PrintProperty(const char* str) {

    if ( outType == 0 /* Bus */) { 
        outBus->PrintProperty(str);
        if ( inBusWidth > 0 ) {
            inBus->PrintProperty(str);
        }
    } else if ( outType == 1 /* LinearArray */) {
        printf("hi\n");
    }

}





