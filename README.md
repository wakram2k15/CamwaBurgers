# Feedback Stabilization and Finite Element Error Analysis of Viscous Burgers Equation around Non-Constant Steady State

**Author:** Wasim Akram (2025)

## Description
This repository contains MATLAB code developed to reproduce the numerical results presented in the manuscript:

> *Feedback Stabilization and Finite Element Error Analysis of Viscous Burgers Equation around Non-Constant Steady State*, submitted to *Computers and Mathematics with Applications (CAMWA)*.

The repository includes implementations of:
1. Finite element discretization of the viscous Burgers equation around a prescribed non-constant steady state.  
2. Computation of feedback gains via the discrete algebraic Riccati equation.  
3. Stabilization of both linear and nonlinear systems using feedback control.  
4. Error and convergence analysis corresponding to the theoretical results in the manuscript (Theorems 3.14–4.14).

---

## Repository Structure

Note: Please reduce or increase the value of nl (number of refinement) to get fewer results in short 
      time or to get results in better refinement, respectively, in each main code. 
The repository contains **three main folders**, corresponding to the examples presented in the paper:

```
├── Example_1_Linear_System                    	# Example 1: Linear feedback stabilization test
│   ├── Main_Figure_2_LinVis_NoCtrl.m   	# To check unstable solution (without control)
│   ├── Main_Figure_2c_LinVis_NoCtrl.m  	# To check rates for unstable solution (without control)
│   ├── Main_Figure_3_LinVis_Stab.m     	# To check stabilizability using feedback gain
│   ├── Main_Table_1_Figure_3c_LinVisStab.m     # To see errors and rates for stabilized solution/control
│   ├── *Supporting scripts*  
│
├── Example_2_Nonlinear_NoCtrl       		# Example 2
│   ├── Main_Table_3_Figure_4.m       		# To observe errors and rates for test example (Table 3 and Figure 4)
│   ├── *Supporting scripts*  			
│
├── Example_3_Stab_NonLinear_System   		# Example 3: Stabilizability, Error estimation and convergence
│   ├── Main_Figure_5_NoCtrl.m       		# Behavious of unstable solution (Figure 5)
│   ├── Main_Table_4_Figure_6_NoCtrl.m		# Errors and Rates for unstable solution (Table 4, Figure 6)
│   ├── Main_Figures_7_8_Stab.m			# Behaviour of stabilized solution (Figures 7 & 8)
│   ├── Main_Table_5_8a_Figure_8c_Stab.m	# Errors in stabilized solution/control/tame cost (Tables 5 & 8a, Figures 8c)
│   ├── Main_Figures_9_10_CtrlSubDom.m		# Behaviour of stabilized solution using prop sub-control-domain (Figures 9 & 10)
│   ├── Main_Table_6_CtrlSubDom.m		# Errors and rates using proper sub-control-domain (Table 6)
│   ├── Main_Figure_11_ProjBased_Stab.m		# Behaviour of stabilized solution using projection based algorithm (Figure 11)
│   ├── Main_Table_7_8b_ProjBased_stab.m	# error, rates, and time-cost using projection based algorithm (Tables 7 & 8b)
│   ├── *Supporting scripts*  			# (error computation, rate estimation, data export)
│
└── README.md
```

Each folder contains:
- The **main script** (`main_exampleX.m`) that reproduces the results from the corresponding section in the paper.  
- Several **supporting functions/scripts** used internally for assembling matrices, solving Riccati equations, computing errors, and generating plots.  

---

## Software Requirements
- MATLAB R2015b (or later)  
- OS: Windows / Linux / macOS

---

## How to Run
1. Open any of the example folders (`Example1_Linear`, `Example2_Nonlinear`, or `Example3_ErrorAnalysis`).  
2. Run the main driver script (`main_exampleX.m`) to reproduce the results.  
3. The results (stabilized solutions, error plots, and convergence tables) will be shown in command window or a pop-up window.

---

## Expected Output
- Time evolution plots showing exponential stabilization/unstable of the system.  
- Tables of L² and H¹ errors showing quadratic and linear, respectively, convergence (as reported in the paper).  
- Numerical values may differ slightly (≈5-7%) from the manuscript due to mesh resolution, solver tolerances, and *random initial guess*.

---

## Reproducibility Note
Minor numerical deviations (≤5%) may occur depending on:
- Mesh density, time step size  
- random initial guess  
- Round-off errors and system precision  

These do **not** affect the qualitative behavior or convergence rates reported in the paper.

---

---

## Citation
If you use this code, please cite:
> Wasim Akram (2025). *Feedback Stabilization and Finite Element Error Analysis of Viscous Burgers Equation around Non-Constant Steady State*, 
	Computers and Mathematics with Applications, arXiv:2406.01553, https://arxiv.org/abs/2406.01553

---

## Code DOI
`[To be added after upload — e.g., https://doi.org/10.5281/zenodo.xxxxxxx]`
