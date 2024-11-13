# mito-agents
ATP reaction-diffusion with mitochondrial agents

Experiments: (0) uniform mitos, uniform consumption; (1) clustered mitos, uniform consumption;
                    (2) uniform mitos, clustered consumption; (3) clustered mitos, clustered consumption

`mito-agents.c` uses a crude PDE solver to compute ATP concentration landscapes in the cell under these conditions, with different, biologically-plausible parameterisations of ATP consumption, production, and more. It produces a set of summary statistics over time for each parameterisation, and some explicit snapshots of the cellular state for each parameterisation. `plots.R` analyses the output, asking under which circumstances pronounced ATP gradients can occur, and looking at the dependencies in the system.

The code parallelises over experiments, taking a command-line parameter to specify the experiment. The clustered consumption experiments equilibrate rapidly; the others take longer (perhaps hours). Run `run.sh` to compile and run the code (or, if not using Bash, look at the file to see how to do so).

Nuances: physical values for cellular ATP concentration (given 2D picture, need to assume a thickness, maybe 10um), reasonable values for total ATP consumption rate, diffusion; how much ATP gradient exists at equilibrium? Need concentration around 1mM; consumption rate around 10^9 molecules / s. Control parameters: consumption rate constant, production rate constant per mito. Hard to find a combination that supports a large gradient while hitting those values.
