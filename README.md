# mito-agents
ATP reaction-diffusion modelling with mitochondrial agents, exploring the conditions under which substantial ATP gradients emerge in the cell

![image](https://github.com/user-attachments/assets/f6fc5c09-b3fa-4866-8cbd-3487a80f3f9c)

Default pipeline
-----

Run `run.sh` with one of the arguments below to compile and run the code for a particular experiment (or, if not using Bash, look at the file to see how to do so):

* `default`     -- default model structure
* `smaller`     -- smaller cell size
* `thinner`     -- thinner vertical cell section
* `fatter`      -- fatter vertical cell section
* `more`        -- more mitochondria
* `fewer`       -- fewer mitochondria
* `coarser`     -- coarser grid elements in numerical scheme
* `postrev`     -- following peer review, ATP-independent diffusion and "fibrous" ATP consumption

After the code completes (an hour or two on a modern machine, if 8 cores are available to run the 8 experimental designs in parallel), run `plots.R` to analyse and visualise the output. `plots-postrev.R` produces plots showing the `postrev` content above. In future these pipelines will be pulled together for parsimony.

Details
----

`mito-agents.c` uses a simple PDE solver to compute ATP concentration landscapes in the cell under different conditions, with different, biologically-plausible parameterisations of ATP consumption, production, and more. It produces a set of summary statistics over time for each parameterisation, and some explicit snapshots of the cellular state for each parameterisation. `plots.R` analyses the output, asking under which circumstances pronounced ATP gradients can occur, and looking at the dependencies in the system. `plot-specific.R` is (old) quicker code for visualising a single snapshot.

Compile the code with, for example (using `gcc`):

`gcc mito-agents.c -lm -o mito-agents.ce`

The code parallelises over experiments, taking command-line parameters to specify several parameters. In order, these are cell width, cell depth, number of mitochondria (max 1000), and how many simulation elements correspond to 1μm, the experiment to run (see below), and (for some experiments) a motion scale parameter. 

The label determining the experiment to run is a little complicated at the moment thanks to the chronology of the project. We'll simplify this in future.

(0) uniform mitos, uniform consumption; (1) uniform mitos, clustered consumption; (-1) uniform mitos, fibrous consumption; (2) clustered mitos, uniform consumption; (3) clustered mitos, clustered consumption; (-3) clustered mitos, fibrous consumption; (4) mobile-1 mitos, uniform consumption; (5) mobile-1 mitos, clustered consumption; (-5) mobile-1 mitos, fibrous consumption; (6) mobile-2 mitos, uniform consumption; (7) mobile-2 mitos, clustered consumption; (-7) mobile-2 mitos, fibrous consumption; (8) mobile-3 mitos, uniform consumption; (9) mobile-3 mitos, clustered consumption; (-9) mobile-3 mitos, fibrous consumption.

mobile-1 means ATP-dependent diffusion; mobile-2 means directed down the ATP gradient; mobile-3 means ATP-independent diffusion. These three cases require a motion parameter scaling either the diffusion kernel or the gradient influence on motion. 

For example:

`./mito-agents.ce 50 10 100 2 -9 0.0125`

runs the default case of a 50x50x10μm cell with 100 mitochondria and simulation elements of width 0.5μm. The experiment run is ATP-independent mitochondrial diffusion with fibrous ATP consumption [-9], with a motion parameter of 0.0125 giving the scale of the diffusion kernel. See `run.sh` for the default parameterisations used in the project.


