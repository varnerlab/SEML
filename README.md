# SEML
SEML: A Simplified English Modeling Language for Constructing Biological Models in Julia

## Introduction to SEML
SEML compiler is a code generator that transforms simplified English description of biological network into modeling written in a runnable programming language, such as the [Julia](http://julialang.org), [Python](https://www.python.org) and [Matlab](https://www.mathworks.com/products/matlab.html).
See our [Github page](https://github.com/varnerlab/SEML.git)



## Requirement
In order to use SEML, the user needs to [install Julia](https://julialang.org/downloads/platform.html) first. This version is compatible with [Julia v1.1.0](https://julialang.org/downloads/index.html).



## Installation
Within [Julia](http://http://julialang.org), use the `add` command of the package manager to download and install the SEML repository:

```
using Pkg
]
(v1.1) pkg> add https://github.com/varnerlab/SEML.git
```

To delete the SEML package use the command:

```
]
(v1.1) pkg> rm SEML
```



## Usage
To use SEML in your project simply issue the command:

```
using SEML
make_model(in_path[, out_path=arg2, host=arg3, model=arg4, lang=arg5])
```

The ``make_model()`` command takes five arguments:

Argument | Required | Description | Options | Default
:--- | :--- | :--- | :--- | :---
in_path | Yes | path to input file | | none
out_path | No	| path to where files are written | | ./autogeneratedmodel
host | No	| host type  | bacterial \| mammalian | bacterial
model | No | model type | FBA \| FVA \| Kinetics | Kinetics
lang | No | programming language | julia \| python \| python2 \| python3 \| matlab | julia



## Output Files
### in Julia
FBA/FVA model files:

File | Description
:--- | :---
DataDictionary.jl | contains all data of the model, including flux bounds, species concentration bounds, and objective coefficients, etc.
FluxDriver.jl  |  interface with GLPK solver to run FBA
FVA.jl  |  implementation of fastFVA from [gudmundsson2010computationally](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-489)
InputFile | SEML description of the model
include.jl | contains all the include statements for the project.
Solve.jl | interface to run the simulation
stoichiometry.dat  | Stoichiometric matrix of the model
Utility.jl  |  provide some auxiliary functions

Kinetic model files:

File | Description
:--- | :---
Balances.jl | encodes mass balance of the model
DataDictionary.jl | contains all data of the model, including initial condition, kinetic constants and Monod affinity constants, etc.
InputFile | SEML description of the model
include.jl | contains all the include statements for the project.
Kinetics.jl | calculates kinetic rates
SolveBalances.jl | interface to run the simulation
stoichiometry.dat  | Stoichiometric matrix of the model

### in Python
FBA/FVA model files:

File | Description
:--- | :---
DataDictionary.py | contains all data of the model, including flux bounds, species concentration bounds, and objective coefficients, etc.
FluxDriver.py  |  interface with scipy.optimize.linprog to run FBA
FVA.py  |  implementation of fastFVA from [gudmundsson2010computationally](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-489)
InputFile | SEML description of the model
Solve.py | interface to run the simulation
stoichiometry.dat  | Stoichiometric matrix of the model

Kinetic model files:

File | Description
:--- | :---
Balances.py | encodes mass balance of the model
DataDictionary.py | contains all data of the model, including initial condition, kinetic constants and Monod affinity constants, etc.
InputFile | SEML description of the model
Kinetics.py | calculates kinetic rates
SolveBalances.py | interface to run the simulation
stoichiometry.dat  | Stoichiometric matrix of the model

### in Matlab
FBA/FVA model files:

File | Description
:--- | :---
DataDictionary.m | contains all data of the model, including flux bounds, species concentration bounds, and objective coefficients, etc.
FluxDriver.m  |  interface with GLPK solver to run FBA
FVA.m  |  implementation of fastFVA from [gudmundsson2010computationally](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-489)
InputFile | SEML description of the model
maximizeProductDictionary.m | modify the objective coefficient array.
Solve.m | interface to run the simulation
stoichiometry.dat  | Stoichiometric matrix of the model

Kinetic model files:

File | Description
:--- | :---
Balances.m | encodes mass balance of the model
DataDictionary.m | contains all data of the model, including initial condition, kinetic constants and Monod affinity constants, etc.
InputFile | SEML description of the model
Kinetics.m | calculates kinetic rates
SolveBalances.m | interface to run the simulation
stoichiometry.dat  | Stoichiometric matrix of the model



## Reproduce Examples
The following Julia packages are required to run the demontration examples in Julia:
* [GLPK](https://github.com/JuliaOpt/GLPK.jl)
* [ODE](https://github.com/JuliaDiffEq/ODE.jl)
* [PyPlot](https://github.com/JuliaPy/PyPlot.jl)

They can be installed by running in Julia:

```
using Pkg
]
(v1.1) pkg> add GLPK ODE PyPlot
```
### reproduce directly
1. Download the `example` folder
2. Open a terminal window, `cd` to the `example` folder
3. To run the FBA model (Fig. 2):
  - `cd` to `fauto` folder
  - run `julia Solve.jl`
3. To run the kinetic model (Fig. 5):
  - `cd` to `kiauto` folder
  - run `julia Simulation.jl`
  (Note that while running the kinetic model, some figure(s) will be generated and pause the program, close the unwanted figure windows to get the final color figures as Fig. 5.)

### reproduce from scratch
1. Download the `example/testcase` folder to get two models in SEML
2. Open a terminal window, go into `julia` and run `using SEML`
3. To build the FBA model (Fig. 2):
  - run `make_model("PathTo/fbacase.txt", model="FBA", out_path="outpath1")`
  - change the following parameters to corresponding values in the generated `DataDictionary.jl`
  ```
  default_bounds_array = [
	...
	0 1.0; # 2 1.0*m_A_c<catalyze:>1.0*m_B_c
	0 1.0; # 3 1.0*m_B_c<catalyze:>1.0*m_A_c
	0 2.0; # 4 1.0*m_A_c<catalyze:>1.0*m_C_c
        ...
  ]
  ...
  species_bounds_array = [
	-10.0 10.0; # 1 m_A_e
	-10.0 10.0; # 2 m_B_e
	-10.0 10.0; # 3 m_C_e
	 ...
  ]
  ```
  - `cd` to `outpath1` folder
  - run `julia Solve.jl`

4. To build the kinetic model (Fig. 4 and 5):
  - run `make_model("PathTo/kinetcase.dat", out_path="outpath2")`
  - change the following parameters to corresponding values in the generated `DataDictionary.jl`
  ```
  kcat_signaling = ones(7)  # kcat[#reaction]: reaction name
  kcat_signaling[1] = 1.1e-3  # kcat: 1.0*m_A_e<uptake:1.0*p_TA_c>1.0*m_A_c
  kcat_signaling[2] = 8e-4  # kcat: 1.0*m_B_c<secrete:1.0*p_TB_c>1.0*m_B_e
  kcat_signaling[3] = 9e-4  # kcat: 1.0*m_C_c<secrete:1.0*p_TC_c>1.0*m_C_e
  kcat_signaling[4] = 1.8e-4  # kcat: 1.0*m_A_c<catalyze:1.0*p_E1_c>1.0*m_B_c
  ...
  W_value_dict = Dict{String, Float64}()
  W_value_dict["W~m_B_c~mRNA_E4_c"] = 0.1
  W_value_dict["W~m_A_c~mRNA_E2_c"] = 0.7
  W_value_dict["W~m_C_c~mRNA_E3_c"] = 0.2
  W_value_dict["W~m_A_e~mRNA_TA_c"] = 0.1
  ...
  W_value_dict["W~m_A_c~mRNA_E1_c"] = 0.6
  W_value_dict["W~m_C_c~mRNA_TC_c"] = 1.1
  ```
  - `cd` to `outpath2` folder
  - run `julia Simulation.jl`
  (Note that while running the kinetic model, some figure(s) will be generated and pause the program, close the unwanted figure windows to get the final color figures as Fig. 5.)



## Support or Contact

Having trouble at installation or function? Feel free to contact the [authors](https://github.com/varnerlab).
