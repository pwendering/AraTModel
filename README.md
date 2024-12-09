## Metabolic modeling of temperature effects on growth and metabolism of _Arabidopsis thaliana_

### Requirements

The code was tested on both Windows 10 and Ubuntu 20.04.3 LTS.

To execute the code, the following software must be installed

* Matlab (tested with 2020a/b)
* [Gurobi solver](https://www.gurobi.com/) (tested with version 9.1.1)
* [COBRA toolbox](https://opencobra.github.io/cobratoolbox/stable/) for Matlab
* (only for statistics of growth experiment) R programming language (tested with version 4.1.2)

### Setup

Clone this repository by running the following code in the command line
```
git clone https://github.com/pwendering/AraTModel
```

Change into the `AraTModel`directory and run `config` to add the functions to the Matlab search path.

### Constraint-based analysis of temperature effects

The function `simulateTempEffects.m` adds temperature-dependent constraints on 
* k<sub>cat</sub> values,
* total protein content,
* photosynthesis (FvCB model) and CO<sub>2</sub> uptake,

and solves the resulting quadratically-constrained optimization problem(s).

The function takes two required input arguments:
* `TGEM`, which a COBRA-formated metabolic model with additional fields required for the addition of temperature-dependent constaints
	- the configured model based on the enzyme-constrained AraCore model can be loaded by running
	```
	load(config('tgemFile'))
	```
* `I`, the light intensity that should be used in the simulation

By default, a simulations will be carried out at temperatures between 10 °C and 40 °C.

Optional input arguments can be given as key-value pairs. For instance, a custom temperature range (20 °C to 30 °C) can be used:
```
simulationResult = simulateTempEffects(TGEM, I, 'tempRange', celsius2kelvin(20:30))
```

The function returns up to four results:
* `simulationResult`: growth rates (`mu`) and net CO<sub>2</sub> assimilation rates (`A`) and the underlying flux distributions from with and without minimization of total flux
	- For `mu` and `A`, the minimum, average, and maximum values from 20 pool solutions are returned. For the publication, the maximum values were used.
* `photParams`: temperature-adjusted values that were used to parametrize the constraints based on the FvCB model
* `fvcbParams`: temperature-adjusted values of obtained by the original FvCB model (Farquhar et al. 1980, [10.1007/BF00386231](https://doi.org/10.1007/BF00386231))
* `gurobiProblem`: pFBA problem with temperature-dependent constaints formatted for the Gurobi solver
	- this is the problem solved for the _last_ temperature if a range is given as input
	
### Configuration of a metabolic model for temperature-dependent analysis

1. Create an enzyme-constrained model using [GECKO 2.0.2](https://github.com/SysBioChalmers/GECKO/tree/v2.0.2) (not compatible with later versions). See script `create_a_th_ec_model.m` for specific changes introduced to the AraCore model. Specific adjustments made to functions of the GECKO toolbox can be obtained upon request.

2. Adjust `config` file if necessary
* input model file
* UniProt protein information file
* photosynthesis parameters
* protein content data
* IDs of some relevant exchange reactions
* (predicted) optimal temperatures
	- can be generated using the command line application at [https://github.com/pwendering/topt-predict](https://github.com/pwendering/topt-predict)
* file containing TPP data (see `protein-stability/ath-protein-stability.csv` for the format)
	- if downloaded from the [Meltome Atlas](https://meltomeatlas.proteomics.wzw.tum.de/master_meltomeatlasapp/), adapt and use `code/bash/get-protein-stability-data-ath` to create the file

3. Add additional fields using the `createTGEM` function (see also `code/matlab/thermo_model/create_etc_model.m` for an example)

### Reproduction of published results

All scripts that generated the results presented in the publication are located at `code/matlab/analysis_scripts`, `code/matlab/plotting`, and `code/R`.

### Reference
