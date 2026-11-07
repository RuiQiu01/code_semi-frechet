This folder contains the reproducibility materials for the manuscript "Efficient Semi-supervised Frechet Regression with Low-dimensional Manifold Structures". It includes simulation scripts, real data scripts, saved result files, plotting scripts, and two R Markdown workflow files.

The workflow files summarize how the saved results and scripts correspond to the manuscript figures.

# Folder Structure

```
.
|-- Readme.md
|-- simulation_workflow.Rmd
|-- realdata_workflow.Rmd
|-- simulation/
|   |-- code/       # original R scripts for simulation experiments
|   |-- plot/       # plotting scripts for simulation figures
|   `-- results/    # saved simulation CSV files
`-- real data/
    |-- code/       # original R scripts and face_data.mat
    |-- plot/       # plotting scripts for real data figures
    `-- results/    # saved real data CSV files
```

# Software

The experiments were run in R version 4.5.0.

The R Markdown workflow files use base R for reading saved result files and summarizing reported numerical values. Running the experiment and plotting scripts requires the packages loaded by those scripts:

```
parallel
pbapply
igraph
Matrix
matrixStats
shapes
RiemBase
tmvtnorm
Directional
Riemann
spherepc
FNN
pracma
R.matlab
ggplot2
ggpubr
plot3D
scatterplot3d
rmarkdown
knitr
```

Install missing packages before running full computations or rendering the notebooks.

# Quick Start

1. Open R from the root of this folder.
2. Render the simulation workflow notebook:

```r
rmarkdown::render("simulation_workflow.Rmd")
```

3. Render the real data workflow notebook:

```r
rmarkdown::render("realdata_workflow.Rmd")
```

4. Use the rendered files to identify which saved result files, scripts, and settings correspond to each manuscript figure.

# Reproduction Workflow

The first step is to render the two R Markdown files. They list the saved result files, provide figure-to-script mappings, summarize the reported real data numerical results, and give commands for regenerating figures using the manuscript plotting scripts.

For figure regeneration from saved results, use the plotting scripts in:

```
simulation/plot/
real data/plot/
```

For reruns from scratch, use the scripts in:

```
simulation/code/
real data/code/
```

The R Markdown tables list the figure-to-script mapping and the relevant parameter settings. Rerun chunks are set to `eval=FALSE` because these computations can be slow. The `sessionInfo()` output records the R environment used to render each workflow file; the saved experimental results were generated using R version 4.5.0.

# Runtime Notes

Most full simulation and real data scripts use parallel computation with up to 25 cores:

```r
cores <- min(25, parallel::detectCores() - 1)
```

Adjust this value manually if your machine has fewer cores or limited memory.

The saved CSV files are included so that the manuscript figures can be checked without rerunning all Monte Carlo experiments.

# Version Control

A `.gitignore` file is included for common local R, R Markdown, operating-system, and temporary files. It does not ignore saved CSV results or manuscript figure files.

If creating a new repository, initialize version control from this folder and commit the reproducibility materials, including the saved results and workflow notebooks.
