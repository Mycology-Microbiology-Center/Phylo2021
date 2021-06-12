## Phylogenetic reconstruction workshop
**University of Tartu, 2021** 

[information is being updated]

Welcome! With this workshop we want to share knowledge on how to run phylogenetic software on an HPC computing cluster, and then how to plot resulting trees in R environment. We will not cover theory of phylogenetic tree estimations, nor go into how to select parameters and priors for phylogenetic software. Instead, the focus is on the most technical aspects of the work. The intended audience is taxonomists and people who need to produce and visualize phylogenetic trees more efficiently. We appreciate if you bring **your own data** formatted similarly to the examples.

**Where:** Tartu, Ravila 14A (Chemicum), room 1100.

**When:** 15.06.2021, 11:00 till 16:00 the latest (we expect the main part to be done in ~3 hours).

### Prerequisites
- Bring your own laptop.
- We'll run analysis on the [Rocket Cluster](https://hpc.ut.ee/en/resources/rocket-cluster-en/) of the University of Tartu. Beforehand, you need to request an account from support@hpc.ut.ee (provide your University of Tartu username). If you are a student, cite and cc: your supervisor.
- HPC part: For Windows users, we suggest to use [WinSCP](https://winscp.net/eng/downloads.php) and [PuTTY](https://winscp.net/eng/downloads.php#putty).
Alternatively, on Windows 10 one may try to use WSL - Windows Subsystem for Linux (see installation guide [here](https://docs.microsoft.com/en-us/windows/wsl/install-win10)).
For Linux and Mac users, it would be possible to use a terminal emulator.
- `R` part: It's most convenient to work from [RStudio](https://www.rstudio.com/products/rstudio/download/) or a similar application of your choice. Please, make sure you have a reasonably fresh version of R itself. You will need the following R libraries: `ggtree`, `ape`,`phylotools`, `ggplot2`, `Hmisc`, `gdata`, `pals`, `geiger` (can be installed during the workshop).
- TreeGraph 2 viewer might be needed for operations with Bayesian trees. Please install it from http://treegraph.bioinfweb.info/

### Course materials
1. [Environment setup on HPC](00.Environment_setup.md)
2. [Scheduling jobs on HPC](01.SLURM.md)
3. [Phylogenetic tree building on HPC](02.Phylo_on_HPC.md)
4. TBA

### [Example files](/data)

### Contact
For any organizational questions contact [Iryna Yatsiuk](mailto:iryna.yatsiuk@ut.ee).
