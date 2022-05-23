#include <cmath>
#include <iostream>
#include <stdlib.h>
#include "LayerScheduler.h"
#include "Param.h"

extern Param *param;

/* HierarchyObject Scheduling 
* Scheduler type01 - compact mapping */
vector<vector<double>> LayerScheduler::HObjectScheduling_01(const HierarchyObject *hObject,
                                double idxRow, double idxCol, 
                                double numRowObjectAvailable, double numColObjectAvailable,
                                double wRow, double wCol, double wInC, double wOutC) {

    /* Check Mapping option of the kernel - kernel split in spatial dim or channel dim */
    // Row-dim
    bool spatialRowMapping, spatialColMapping;
    if ( hlevelRowSystolic < 0 ) { // no spatial row-dim mapping for the entire mapping
        spatialRowMapping = false; 
        wInC = wRow * wInC;
        wRow = 1;
    } else if ( hObject->hlevel == hlevelRowSystolic ) {
        spatialRowMapping = true;
    } else if ( (hObject->hlevel == hlevelRowSystolic + 1) && (numHObjectRowSAExt > 1) ) {
        /* NOTE: Systolic spill is allowed for one more higher hObect */
        spatialRowMapping = true;
    } else {
        spatialRowMapping = false;
    }
    // Col-dim
    if ( hlevelColSystolic < 0 ) { // no spatial col-dim mapping for the entire mapping
        spatialColMapping = false; 
        wInC = wCol * wInC;
        wCol = 1;
    } else if ( hObject->hlevel == hlevelColSystolic ) {
        spatialColMapping = true;
    } else if ( (hObject->hlevel == hlevelColSystolic + 1) && (numHObjectColSAExt > 1) ) {
        /* NOTE: Systolic spill is allowed for one more higher hObect */
        spatialColMapping = true;
    } else {
        spatialColMapping = false;
    }

    /* Check and update residual information for the subObjects */
    // idxOffset info : idxSubTop, outC_split info (2), inC_split info (2) (total length: 5)
    // split info     : numRowAvailable, numColAvailable
    vector<int> idxOffset = {0, 0, 0, 0, 0};

    double idxOffsetRow = 0;
    double idxOffsetCol = 0;
    double numRowSubObjectAvailable = 0;
    double numColSubObjectAvailable = 0;
    int split_type = -1;

    if ( hObject->hlevel == hlevelTop ) {
        if ( idxOffsetVector.size() > 0 ) {
            // get offset info    
            idxOffset = idxOffsetVector.back();        
            idxOffsetVector.pop_back();
            idxOffsetRow = idxOffset[0];

            // update split type info ( split info with larger inResidual is selected )
            if ( idxOffset[1] > idxOffset[3] ) {
                split_type = 0; // outC split
                numRowSubObjectAvailable = idxOffset[1];
                numColSubObjectAvailable = idxOffset[2];
            } else {
                split_type = 1; // inC split
                numRowSubObjectAvailable = idxOffset[3];
                numColSubObjectAvailable = idxOffset[4];
            }
        } else { 
            numRowSubObjectAvailable = hObject->numInCSubObject;
            numColSubObjectAvailable = hObject->numOutCSubObject;
        }
    } else {
        double numAvailableSubObjectRow = ceil( numRowObjectAvailable / hObject->numRowSubObject );
        double numAvailableSubObjectCol = ceil( numColObjectAvailable / hObject->numColSubObject );
        idxOffsetRow = hObject->numSubObjectRow - numAvailableSubObjectRow;
        idxOffsetCol = hObject->numSubObjectCol - numAvailableSubObjectCol;
        numRowSubObjectAvailable = numRowObjectAvailable - (numAvailableSubObjectRow - 1) * hObject->numRowSubObject; 
        numColSubObjectAvailable = numColObjectAvailable - (numAvailableSubObjectCol - 1) * hObject->numColSubObject; 
    }

    //printf("h: %d, numObjectAvaialable: (%d x %d), numSubObjectAvailable: (%d x %d), idxOffset: (%d, %d)\n",
    //       hObject->hlevel, (int)numRowObjectAvailable, (int)numColObjectAvailable,
    //       (int)numRowSubObjectAvailable, (int)numColSubObjectAvailable, (int)idxOffsetRow, (int)idxOffsetCol);
    //printf("h: %d, numObjectAvaialable: (%d x %d), wInC: %d, wOutC: %d, idxOffset: (%d, %d)\n",
    //       hObject->hlevel, (int)numRowObjectAvailable, (int)numColObjectAvailable,
    //       (int)wInC, (int)wOutC, (int)idxOffsetRow, (int)idxOffsetCol);

    // check the type of subObject address, 0: row/col for input/output extension, 1: linear address
    int idxType = ( hObject->inputBuffer || spatialColMapping || spatialRowMapping )? 0: 1;

    // check if IC info
    const InterConnect *IC = hObject->interConnect;
    int icOutType = IC->outType;
    int icInType = IC->inType;
    BusMode icInBusMode = IC->inBusMode;

    // generate vector
    vector<vector<double>> infoReadVector;

    //*** [STEP1] Scheduling Sub-HierarchyObject ***//

    // get & update sub-object scheduling results
    int subObjectCounter = 0;
    int idxSubObjectRow = 0;
    int idxSubObjectCol = 0;
    double wRowSubObject, wColSubObject, wInCSubObject, wOutCSubObject;
    double wInC_loop, wOutC_loop;
    double wInCSubObjectAvailable, wOutCSubObjectAvailable;
    double wInCSubObjectResidual, wOutCSubObjectResidual;
    // for mesh scheduling
    vector<double> infoReadICIn, infoReadICOut;
    vector<vector<double>> infoReadICInVector, infoReadICOutVector;
    double numOutBitSubObject, numInBitSubObject;
    int prevIdxSubObjectRow;
    double icState = 1; /*reset before update*/
    int numSubObjectCol = 0;
    int numSubObjectRow = 0;
    int subObjectRow = 0;
    int subObjectCol = 0;
    bool next_offset = false;

    wInC_loop = wInC;

    while ( wInC_loop > 0 ) { // NOTE: not working for the systolic array
        // update size of weight matrix assigned for each subObject
        // Row-dim
        //// NOTE: modify for spatialMapping with Systolic Array Extension
        //// NOTE: modify for inter-layer SubObject splitting
        if ( spatialRowMapping ) {
            wRowSubObject = 1; // assume only 1 pixel is assigned for spatial dim mapping
            wInCSubObject = wInC;
            wInCSubObjectResidual = 0;
        } else {
            wRowSubObject = wRow;
            if ( (subObjectRow == 0) | (hObject->hlevel == hlevelTop) ) {
                wInCSubObjectAvailable = numRowSubObjectAvailable;
            } else {
                wInCSubObjectAvailable = hObject->numInCSubObject;
            }
            wInCSubObject = MIN(wInCSubObjectAvailable, wInC_loop); 
            wInCSubObjectResidual = wInCSubObjectAvailable - wInCSubObject;
        }
        //printf("wInCSubObject: %d, wInCSubObjectAvaialable: %d, wInCSubObjectResidual: %d\n", 
        //        (int)wInCSubObject, (int)wInCSubObjectAvailable, (int)wInCSubObjectResidual);

        // update remained wInC_loop
        wInC_loop = wInC_loop - wInCSubObject;
        // update wOutC_loop
        wOutC_loop = wOutC;
        subObjectCol = 0;
        numSubObjectCol = 0;

        // for output extension
        // for (int subObjectCol=0; subObjectCol < (int)numSubObjectCol; subObjectCol++) {
        while ( wOutC_loop > 0 ) {
            // update size of weight matrix assigned for each subObject
            // Col-dim
            if ( spatialColMapping ) {
                wColSubObject = 1; // assume only 1 pixel is assigned for spatial dim mapping
                wOutCSubObject = wOutC;
            } else {
                wColSubObject = wCol;
                if ( (subObjectCol == 0) | (hObject->hlevel == hlevelTop) ) {
                    wOutCSubObjectAvailable = numColSubObjectAvailable;
                } else {
                    wOutCSubObjectAvailable = hObject->numOutCSubObject;
                }
                wOutCSubObject = MIN(wOutCSubObjectAvailable, wOutC_loop); 
                wOutCSubObjectResidual = wOutCSubObjectAvailable - wOutCSubObject;


            }
            wOutC_loop = wOutC_loop - wOutCSubObject;

            // get index 
            if ( idxType == 0 ) {
                idxSubObjectRow = idxOffsetRow + subObjectRow;
                idxSubObjectCol = idxOffsetCol + subObjectCol;
            } else { // idxType == 1
                div_t divResult = div((subObjectCounter + idxOffsetRow), hObject->numSubObjectCol);
                idxSubObjectRow = divResult.quot;
                idxSubObjectCol = divResult.rem;
            }

            //// debugging
            //printf("hl: %d, idxOffRow: %d, idxSubObjectRow: %d, idxSubObjectCol: %d, wInCSubObject: %d, wOutCSubObject: %d\n", 
            //        hObject->hlevel, (int)idxOffsetRow, (int)idxSubObjectRow, (int)idxSubObjectCol, 
            //        (int)wInCSubObject, (int)wOutCSubObject);
            //printf("hl: %d, numSubObjectRow: %d, numSubObjectCol: %d, wInCSubObject: %d, wOutCSubObject: %d\n", 
            //        hObject->hlevel, numSubObjectRow, numSubObjectCol, 
            //        (int)wInCSubObject, (int)wOutCSubObject);

            // check if the idx is valid
            if ( (idxSubObjectRow < 0) || (idxSubObjectRow >= hObject->numSubObjectRow) ||
                 (idxSubObjectCol < 0) || (idxSubObjectCol >= hObject->numSubObjectCol) ) {
                //// debugging
                //printf("hlevel: %d, idxSubObjectRow: %d, idxSubObjectCol: %d, subObjectCounter: %d\n", 
                //        hObject->hlevel, idxSubObjectRow, idxSubObjectCol, subObjectCounter);
                cerr << "[Error] Invalid idxSubObjectRow/Col!" << endl;
                exit(-1);
            }

            // scheduling sub-object 
            if ( hObject->hlevel == 1 ) {
                if ( (wRowSubObject > 1) || (wColSubObject > 1) ) {
                    cerr << "[Error] wRow/Col for hRootObject should be 1!" << endl;
                    exit(-1);
                }
                vector<double> infoReadSub = HRootScheduling( hObject->rootObject, idxSubObjectRow, idxSubObjectCol, 
                                                            wInCSubObject, wOutCSubObject);
                infoReadVector.push_back(infoReadSub); // update scheduling result
                infoReadSub.clear(); // clear updated result
            } else {
                vector<vector<double>> infoReadVectorSub 
                                            = HObjectScheduling_01( hObject->subObject, idxSubObjectRow, idxSubObjectCol, 
                                                            wInCSubObjectAvailable, wOutCSubObjectAvailable, 
                                                            wRowSubObject, wColSubObject, wInCSubObject, wOutCSubObject);
                infoReadVector.insert(infoReadVector.end(), infoReadVectorSub.begin(), infoReadVectorSub.end()); // update
                infoReadVectorSub.clear(); // clear updated result
            }

            // update infoReadIC
            // (layerIdx, hlevel, idxRow, idxCol, state, dataType, numRead, numHopsRow, numHopsCol, dummy)
            if ((icOutType == 2 /*2D Mesh*/) || (icOutType == 3 /*HBus*/)) {
                numOutBitSubObject = numConv * wOutCSubObject / numCellPerSynapse * hObject->numOutBitSubObject;
                infoReadICOut = {(double)layerIdx, (double)hObject->hlevel, idxRow, idxCol, icState, 0/*dataType*/,
                                numOutBitSubObject / IC->outBusWidth, (double)idxSubObjectRow, 0, 0}; // no col dim move of out
                infoReadICOutVector.push_back(infoReadICOut);
                icState = 2; /*cumulate*/
            } // infoReadICOut

            if ((icInType == 2 /*2D Mesh*/) || (icInType == 3 /*HBus*/)) {
                if ( subObjectCol == 0 ) { // input switch -> get feeding info of new input
                    if ( hlevelColSystolic > 0 ) {
                        if ( hObject->hlevel >= hlevelColSystolic ) {
                            numInBitSubObject = outH * ( wRowSubObject * inW * wInCSubObject) * numBitInput;
                        } else {
                            numInBitSubObject = outH * ( wRowSubObject * outW * wInCSubObject) * numBitInput;
                        }
                    } else {
                        numInBitSubObject = numConv * ( wRowSubObject * wColSubObject * wInCSubObject) * numBitInput;
                    }
                    infoReadICIn = {(double)layerIdx, (double)hObject->hlevel, idxRow, idxCol, 2/*icState*/, 1/*dataType*/,
                                numInBitSubObject / IC->inBusWidth, (double)idxSubObjectRow, (double)idxSubObjectCol};
                } else { // input sharing
                    if ( prevIdxSubObjectRow == idxSubObjectRow ) {
                        infoReadICIn.back() = MAX(infoReadICIn.back(), idxSubObjectCol); // horizontal & same row -> update the col idx
                    } else { // rowIdx of input changed -> update prev input info & get new row input feeding info
                        if ( (icInType == 3 /*HBus*/) && (prevIdxSubObjectRow % 2 == 0)  ) { // row bus sharing for HBus
                            infoReadICIn.rbegin()[1] = (double)idxSubObjectRow;
                            infoReadICIn.back() = MAX(infoReadICIn.back(), idxSubObjectCol);
                        }
                        else {
                            infoReadICIn.push_back(0/*dummy*/);
                            infoReadICInVector.push_back(infoReadICIn); // NOTE: BUG: Latency should be overlapped
                            double rowMove = (double)(idxSubObjectRow - prevIdxSubObjectRow);
                            if ( icInType == 2 /*2D Mesh*/ ) {
                                infoReadICIn = {(double)layerIdx, (double)hObject->hlevel, idxRow, idxCol, 2, 1/*dataType*/,
                                        numInBitSubObject / IC->inBusWidth, rowMove, (double)idxSubObjectCol};
                            } else {
                                infoReadICIn = {(double)layerIdx, (double)hObject->hlevel, idxRow, idxCol, 2, 1/*dataType*/,
                                        numInBitSubObject / IC->inBusWidth, rowMove + 1, (double)idxSubObjectCol};
                            }
                        }
                    }
                }
            } // infoReadICIn

            // when top, check if the residual subObject is used or not, and removed from list if it is used
            if ( hObject->hlevel == hlevelTop ) {
                //printf("Hi!\n");
                /* Update nextIdxOffset list if the subObject is residual */
                if ( (wOutCSubObjectResidual > 0) | (wInCSubObjectResidual > 0) ) {
                    // idxOffset info : idxSubTop, outC_split info (2), inC_split info (2) (total length: 5)
                    // split info     : numRowAvailable, numColAvailable

                    // idx of the subObject
                    int idxSubObject_linear = subObjectCounter + idxOffsetRow; // linear index  
                    vector<int> nextIdxOffset = {idxSubObject_linear};

                    // update outC_split offset
                    //printf("%d %d %d %d %d\n", 
                    //        idxOffset[0], idxOffset[1], idxOffset[2], idxOffset[3], idxOffset[4]);
                    if ( wOutCSubObjectResidual > 0 ) {
                        if ( split_type == 0 /*outC_split*/ ) {
                            nextIdxOffset.push_back(idxOffset[1]);
                            nextIdxOffset.push_back(wOutCSubObjectResidual);
                        } else if ( split_type == 1 /*inC_split*/ ){
                            if ( (idxOffset[1] > 0) & (idxOffset[2] > 0) ) {// exist outC_split info
                                nextIdxOffset.push_back(idxOffset[1]);
                                if ( wOutCSubObjectResidual < idxOffset[2] ) {// outC residual is smaller than previous outC info
                                    nextIdxOffset.push_back(wOutCSubObjectResidual);
                                } else {
                                    nextIdxOffset.push_back(idxOffset[2]);
                                }
                            } else {
                                nextIdxOffset.push_back(idxOffset[3]);
                                nextIdxOffset.push_back(wOutCSubObjectResidual);
                            }
                        } else {
                            nextIdxOffset.push_back(hObject->numInCSubObject);
                            nextIdxOffset.push_back(wOutCSubObjectResidual);
                        }
                    } else {
                        nextIdxOffset.push_back(0);
                        nextIdxOffset.push_back(0);
                    }

                    // update inC_split offset
                    if ( wInCSubObjectResidual > 0 ) {
                        if ( split_type == 0 /*outC_split*/ ) {
                            if ( (idxOffset[3] > 0) & (idxOffset[4] > 0) ) {// exist inC_split info
                                if ( wInCSubObjectResidual < idxOffset[3] ) {// inC residual is smaller than previous inC info
                                    nextIdxOffset.push_back(wInCSubObjectResidual);
                                } else {
                                    nextIdxOffset.push_back(idxOffset[3]);
                                }
                                nextIdxOffset.push_back(idxOffset[4]);
                            } else {
                                nextIdxOffset.push_back(wInCSubObjectResidual);
                                nextIdxOffset.push_back(idxOffset[2]);
                            }
                        } else if ( split_type == 1 /*inC_split*/ ){
                            nextIdxOffset.push_back(wInCSubObjectResidual);
                            nextIdxOffset.push_back(idxOffset[4]);
                        } else {
                            nextIdxOffset.push_back(wInCSubObjectResidual);
                            nextIdxOffset.push_back(hObject->numOutCSubObject);
                        }
                    } else {
                        nextIdxOffset.push_back(0);
                        nextIdxOffset.push_back(0);
                    }

                    nextIdxOffsetVector.push_back(nextIdxOffset);
                    nextIdxOffset.clear();
                }

                // check if the next subObject is residual subObject
                // when top, check if the residual subObject is used or not, and removed from list if it is used
                idxOffset = {0, 0, 0, 0, 0};
                bool offsetSubObject = false;

                if ( ( idxOffsetVector.size() > 0 ) /* there is offset that we have to check */
                     & ((wInC_loop > 0) | (wOutC_loop > 0)) /* another subObject should be assigned */) {

                    //printf("In!\n");

                    while ( (!offsetSubObject) & (idxOffsetVector.size() > 0) ) {
                        // get offset info    
                        idxOffset = idxOffsetVector.back();        
                        idxOffsetVector.pop_back();
                        idxOffsetRow = idxOffset[0];
                        subObjectCounter = -1; // 1 is added at the end

                        // update split type info ( split info with larger inResidual is selected )
                        if ( idxOffset[1] > idxOffset[3] ) {
                            split_type = 0; // outC split
                            numRowSubObjectAvailable = idxOffset[1];
                            numColSubObjectAvailable = idxOffset[2];
                        } else {
                            split_type = 1; // inC split
                            numRowSubObjectAvailable = idxOffset[3];
                            numColSubObjectAvailable = idxOffset[4];
                        }

                        // if inC of next subObject is not than the requirement of wOutC_loop
                        // or move to next row, then accept the subObject
                        if ( (wOutC_loop <= 0) | (numRowSubObjectAvailable >= wInCSubObject) ) {
                            offsetSubObject = true;
                        } else {
                            idxOffsetRow += 1;
                        }
                    }

                }

                if ( !offsetSubObject ) {
                    split_type = -1;
                    numRowSubObjectAvailable = hObject->numRowSubObject;
                    numColSubObjectAvailable = hObject->numColSubObject;
                //} else {
                //    printf("Bye! (idxOffsetRow: %d)\n", (int)idxOffsetRow);
                //    printf("h: %d, numSubObjectAvailable: (%d x %d), idxOffset: (%d, %d)\n",
                //           hObject->hlevel,
                //           (int)numRowSubObjectAvailable, (int)numColSubObjectAvailable, (int)idxOffsetRow, (int)idxOffsetCol);
                //    
                }
                wInCSubObjectAvailable = numRowSubObjectAvailable;
                wOutCSubObjectAvailable = numColSubObjectAvailable;
                wInCSubObjectResidual = wInCSubObjectAvailable - wInCSubObject;

            } // hlevelTop

            // move to next col
            subObjectCol = subObjectCol + 1;
            numSubObjectCol = numSubObjectCol + 1;

            // move to next subObject
            prevIdxSubObjectRow = idxSubObjectRow;
            subObjectCounter += 1;

        } // numSubObjectCol

        // update infoReadICIn
        if ((icInType == 2 /*2D Mesh*/) || (icInType == 3 /*HBus*/)) {
            infoReadICIn.push_back(0/*dummy*/);
            infoReadICInVector.push_back(infoReadICIn);
        }

        // move to next row
        subObjectRow = subObjectRow + 1;
        numSubObjectRow = numSubObjectRow + 1;

    }


    // update numUsed SubObject_top
    if ( hObject->hlevel == hlevelTop ) {
        numUsedSubObject_top = subObjectCounter + idxOffsetRow;
    } else {
        // check if there was any overflow in the subObject index
        if ( idxType == 0 ) {
            if ( ( (idxSubObjectRow + 1) > hObject->numSubObjectRow ) || ( (idxSubObjectCol + 1) > hObject->numSubObjectCol) ) {
                cerr << "[ERROR] index of subObject row overflow" << endl;
                cerr << "last idxSubObjectRow: " << idxSubObjectRow <<  ", subObjectRow: " << hObject->numSubObjectRow << endl;
                cerr << "last idxSubObjectCol: " << idxSubObjectCol <<  ", subObjectCol: " << hObject->numSubObjectCol << endl;
                exit(-1);
            }
        } else { // idxType == 1
            int numUsedSub = subObjectCounter + idxOffsetRow;
            if ( numUsedSub > hObject->numSubObject ) {
                cerr << "[ERROR] number of subObject overflow" << endl;
                cerr << "used numSubObject: " << numUsedSub << ", available numSubObject: " << hObject->numSubObject << endl;
                exit(-1);
            }
        }
    }


    // update nextIdxOffsetVector
    if ( hObject->hlevel == hlevelTop ) {
        if ( idxOffsetVector.size() > 0 ) { // there is partially used subObject that was not checked in this layer
            nextIdxOffsetVector.insert(nextIdxOffsetVector.end(), 
                                       idxOffsetVector.begin(), idxOffsetVector.end());
        } else if ( nextIdxOffsetVector.size() == 0 ) { 
            // update the subObject that we have to start with in the next layer
            int idxSubObject_linear = subObjectCounter + idxOffsetRow; // linear index
            vector<int> nextIdxOffset = {idxSubObject_linear,
                                         hObject->numInCSubObject, hObject->numOutCSubObject,
                                         hObject->numInCSubObject, hObject->numOutCSubObject};
            nextIdxOffsetVector.push_back(nextIdxOffset);
            nextIdxOffset.clear();
        } // if there is nextIdxOffset, there is the correct starting point already.
    }

    if ( doingAct ) {
        doneAct = true;
    }
    if ( doingMaxPool ) {
        doneMaxPool = true;
    }

    //*** [STEP1] Scheduling Current HierarchyObject ***//
    vector<double> infoRead; // storage for update scheduling results
    
    // update object info
    infoRead.push_back((double)layerIdx);
    infoRead.push_back((double)hObject->hlevel); // hlevel
    infoRead.push_back(idxRow); // object row index
    infoRead.push_back(idxCol); // object col index
    infoRead.push_back(wRow * wCol * wInC); 
    infoRead.push_back(wOutC); 

    // update infoReadCIM
    infoRead.push_back(0);  // numBitInput
    infoRead.push_back(0);  // numCellPerSynapse
    infoRead.push_back(0);  // numRead Array

    //*** [STEP2] Scheduling DE ***//
    const DigitalElements *DE = hObject->digitalElements;
    double numRead, numUnitAdd;

    // number of total output 
    double numOut = numConv * wOutC / numCellPerSynapse; // NOTE: What for systolic array?  -- seems like the same for adderTree, but different in subObject addition & the size of adderTree in the object

    // adderTree
    if ( DE->placeAdderTree && !doneAct ) {
        numRead = numOut / DE->adderTree->numAdderTree;
        numUnitAdd = numSubObjectRow;
    } else {
        numRead = 0;
        numUnitAdd = 0;
    }
    infoRead.push_back(numRead);
    infoRead.push_back(numUnitAdd);

    // reLu
    if ( DE->placeReLu && (wRow * wCol * wInC == fanIn) && !doneAct) { // has reLu unit & finished MAC operation
        numRead = numOut / DE->reLu->numUnit;
        doingAct = true;
    } else {
        numRead = 0;
    }
    infoRead.push_back(numRead);

    // maxPool
    if ( DE->placeMaxPooling && (doingAct || doneAct) && !doneMaxPool) { // has maxPooling unit & finished MAC operation
        numRead = numOut / DE->maxPooling->window / DE->maxPooling->numMaxPooling;
        doingMaxPool = true;
    } else {
        numRead = 0;
    }
    infoRead.push_back(numRead);
    

    //*** [STEP3] Scheduling IC ***//
    double totalOutBit, totalInBit;
    if ( (hlevelRowSystolic > 0 ) & (hObject->hlevel == hlevelRowSystolic) ) {
        totalOutBit = numOut * hObject->numOutBitSubObject;
    } else {
        totalOutBit = numSubObjectRow * numOut * hObject->numOutBitSubObject;
    }

    if ( hlevelColSystolic > 0 ){ // NOTE:BUG: with stride, add additional FF access count for just passing!!
        if ( hObject->hlevel >= hlevelColSystolic ) {
            totalInBit = outH * (wRow * inW * wInC) * numBitInput;
        } else {
            totalInBit = outH * (wRow * outW * wInC) * numBitInput;
        }
    } else {
        totalInBit = numConv * (wRow * wCol * wInC) * numBitInput;
    }

    
    // update infoReadICVector for output
    if ((icOutType == 2 /*2D Mesh*/) || (icOutType == 3 /*HBus*/)) {
        infoReadVector.insert(infoReadVector.end(), infoReadICOutVector.begin(), infoReadICOutVector.end());
        infoReadICOutVector.clear();
    } else {
        // numOutRead
        double numOutRead = totalOutBit / IC->outBusWidth;
        
        // update IC infoRead // NOTE:BUG: assume single layer on hObject w/ Linear Array
        vector<double> tmpInfoReadIC; 
        if ( IC->outType == 1 /*LinearArray*/) {
            tmpInfoReadIC = {(double)layerIdx, (double)hObject->hlevel, idxRow, idxCol, 1/*icState*/, 0/*dataType*/,
                            numOutRead, (double)numSubObjectRow, (double)numSubObjectCol, 0/*dummy*/};
        } else{
            tmpInfoReadIC = {(double)layerIdx, (double)hObject->hlevel, idxRow, idxCol, 1/*icState*/, 0/*dataType*/,
                            numOutRead, 0, 0, 0/*dummy*/};
        }
        infoReadVector.push_back(tmpInfoReadIC);
        tmpInfoReadIC.clear();
    }

    // update infoReadICVector for input
    if ((icInType == 2 /*2D Mesh*/) || (icInType == 3 /*HBus*/)) {
        infoReadVector.insert(infoReadVector.end(), infoReadICInVector.begin(), infoReadICInVector.end());
        infoReadICInVector.clear();
    } else {
        double numInRead;

        // numInRead
        if ( IC->inType == 0 /*Bus*/) {
            numInRead = totalInBit / IC->inBusWidth;
        } else { /*Linear Array*/
            numInRead = ( totalInBit / wRow ) / IC->inBusWidth;
        }

        // update IC infoRead // NOTE:BUG: assume single layer on hObject w/ Linear Array
        vector<double> tmpInfoReadIC; 
        if ( IC->inType == 1 /*LinearArray*/) {
            tmpInfoReadIC = {(double)layerIdx, (double)hObject->hlevel, idxRow, idxCol, 2/*icState*/, 1/*dataType*/,
                            numInRead, (double)numSubObjectRow, (double)numSubObjectCol, 0/*dummy*/};
        } else{
            tmpInfoReadIC = {(double)layerIdx, (double)hObject->hlevel, idxRow, idxCol, 2/*icState*/, 1/*dataType*/,
                            numInRead, 0, 0, 0/*dummy*/};
        }
        infoReadVector.push_back(tmpInfoReadIC);
        tmpInfoReadIC.clear();
    }

    // empty IC infoRead
    vector<double> tmpInfoReadIC = {0, 0, 0, 0, 0, 0};
    infoRead.insert(infoRead.end(), tmpInfoReadIC.begin(), tmpInfoReadIC.end());
    tmpInfoReadIC.clear();


    //*** [STEP4] Scheduling BU ***//
    const BufferUnit *BU = hObject->bufferUnit;
    double numOutRead, numOutWrite, numInRead, numInWrite, outParallelism, inParallelism;

    // get output data size stored in buffer unit (BU)
    double numOutBitBU, numOutBU;
    if ( doneMaxPool || doingMaxPool ) {
        numOutBU = numOut / DE->maxPooling->window;
        numOutBitBU = numBitInput;
    } else {
        numOutBU = numOut;
        if ( doneAct || doingAct ) {
            numOutBitBU = numBitInput;
        } else {
            if ( DE->placeAdderTree ) {
                numOutBitBU = DE->adderTree->numOutBit;
            } else {
                numOutBitBU = hObject->numOutBitSubObject;
            }
        }
    }
    double totalOutBitBU = numOutBU * numOutBitBU;

    // get number of read 
    if ( BU->inBUSize > 0 ) {
        numOutRead = totalOutBitBU / BU->outBUCoreBW;
        numInRead = totalInBit / BU->inBUCoreBW;
    } else {
        numOutRead = (totalOutBitBU + totalInBit) / BU->outBUCoreBW;
        numInRead = 0;
    }

    // get number of write & parallelism
    if ( BU->buType == 0 /*D Flip-Flop*/ ) {
        numOutWrite = 0;
        numInWrite = 0;
        outParallelism = 1;
        inParallelism = 1;
    } else {
        numOutWrite = numOutRead;
        numInWrite = numInRead;
        // NOTE: separate IC if separate buffer unit for in/out
        if ( BU->inBUSize > 0 ) { // have input buffer
            outParallelism = MIN(BU->numOutBUCore, IC->outBusWidth/BU->outBUCoreBW);
            inParallelism = MIN(BU->numInBUCore, IC->inBusWidth/BU->inBUCoreBW);
        } else { // only output buffer
            outParallelism = MIN(BU->numOutBUCore, (IC->outBusWidth + IC->inBusWidth)/BU->outBUCoreBW);
            inParallelism = 0;
        }
    }
    infoRead.push_back(numOutRead);
    infoRead.push_back(numOutWrite);
    infoRead.push_back(outParallelism);
    infoRead.push_back(numInRead);
    infoRead.push_back(numInWrite);
    infoRead.push_back(inParallelism);
    
    //*** [STEP5] Update scheduling results ***//
    infoReadVector.push_back(infoRead);

    return infoReadVector;
}

