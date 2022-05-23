#ifndef LINEARARRAY_H_
#define LINEARARRAY_H_

#include "typedef.h"
#include "InputParameter.h"
#include "Technology.h"
#include "MemCell.h"
#include "FunctionUnit.h"
#include "DFF.h"

class LinearArray: public FunctionUnit {
public:
	LinearArray(const InputParameter& _inputParameter, const Technology& _tech);
	virtual ~LinearArray() {}
	const InputParameter& inputParameter;
	const Technology& tech;

	/* Functions */
	void PrintProperty(const char* str);
	void Initialize(BusMode _mode, bool _inputWire, int _numRow, int _numCol, double _delaytolerance, double _busWidth, double unitHeight, double unitWidth, double foldedratio, double _clkFreq);
	void CalculateArea();
	void CalculateLatency(int numActiveRow, int numActiveCol, double numRead);
	void CalculatePower(int numActivRow, int numActiveCol, double numBitAccess, double numRead);
	double GetUnitLengthRes(double wireLength);

	/* Properties */
	bool initialized;	/* Initialization flag */
	double widthInvN, widthInvP, wInv, hInv, capInvInput, capInvOutput;
	double widthMinInvN, widthMinInvP, wMinInv, hMinInv, capMinInvInput, capMinInvOutput, wRep, hRep, capRepInput, capRepOutput;
	double AR, Rho, unitLengthWireResistance, minDist, minDelay, resOnRep;
	int numRow, numCol, numRepeater, repeaterSize, numNode, numBranch;
	double unitHeight, unitWidth, wireWidth;
	double busWidth, delaytolerance, unitLengthWireCap, wireLength;
	double unitLatencyRep, unitLatencyWire, repeaterLeakage, unitLengthEnergyRep, unitLengthEnergyWire;
    BusMode mode; /* vertical or horizontal bus */
	double clkFreq;
    double foldedratio;
    double switchingRatio;
   
    bool inputWire;
    int numRepeaterInputWire;
    double inputWireLength;

    DFF dff;
};

#endif /* LINEARARRAY_H_ */
