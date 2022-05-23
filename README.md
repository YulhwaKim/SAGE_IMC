# SAGE_IMC
An open-source framework for the IMC architecture design.

## Requirements
Conda environment required for the experiment is updated in 'environment.yml'.

## Makefile & Program Usage
First, compile the C++ source code.

    cmake CMakeLikst.txt
    make
    
Second, run python script to evaluate previous IMC accelerator desing

    pythons script_evaluate_previous_arch.py
    
Third, run python script to explore design space of IMC accelerators

    pythons script_dse.py
