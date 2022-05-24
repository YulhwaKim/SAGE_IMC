# SAGE-IMC
An open-source framework for the IMC architecture design.

The SAGE-IMC consists of 
1. Architecture Template (./ArchitectureTemplate)
2. DNN Mapper (./Scheduler)
3. Performance Analysis Engine (./CircuitModule)

Note that the performance analysis engine (./CircuitModule) is based on [NeuroSim](https://github.com/neurosim/MLP_NeuroSim_V3.0)[1].

## Requirements
Conda environment required for the experiment is updated in 'environment.yml'.

## Makefile & Program Usage
First, compile the C++ source code.

    cmake CMakeLikst.txt
    make
    
Second, run python script to evaluate previous IMC accelerator design.

    python script_evaluate_previous_arch.py
    
Third, run python script to explore design space of IMC accelerators.

    python script_dse.py

## Reference
[1] P.-Y. Chen, X. Peng, S. Yu, “NeuroSim: A circuit-level macro model for benchmarking neuro-inspired architectures in online learning,” IEEE Trans. CAD,vol. 37, no. 12, pp. 3067-3080, 2018. Source code is Available at: https://github.com/neurosim/MLP_NeuroSim_V3.0.
