## Phylogenetic reconstruction workshop
**University of Tartu, 2021** 


Welcome! With this workshop we want to share knowledge on how to run phylogenetic software on an HPC computing cluster, and then how to plot resulting trees in R environment. We will not cover theory of phylogenetic tree estimations, nor go into how to select parameters and priors for phylogenetic software. Instead, the focus is on the most technical aspects of the work. The intended audience is taxonomists and people who need to produce and visualize phylogenetic trees more efficiently. We appreciate if you bring **your own data** formatted similarly to the [examples](examples.md).

**Where:** Tartu, Ravila 14A (Chemicum), room 1121.

**When:** 15.06.2021, 11:00 till 16:00 the latest (we expect the main part to be done in ~3 hours).

### Prerequisites
- Bring your own laptop.
- We'll run analysis on the [Rocket Cluster](https://hpc.ut.ee/en/resources/rocket-cluster-en/) of the University of Tartu. Beforehand, you need to request an account from support@hpc.ut.ee (provide your University of Tartu username). If you are a student, cite and cc: your supervisor.
- HPC part: For Windows users, we suggest to use [WinSCP](https://winscp.net/eng/downloads.php) and [PuTTY](https://winscp.net/eng/downloads.php#putty).
Alternatively, on Windows 10 one may try to use WSL - Windows Subsystem for Linux (see installation guide [here](https://docs.microsoft.com/en-us/windows/wsl/install-win10)).
For Linux and Mac users, it would be possible to use a terminal emulator.
- R part: It's most convenient to work from [RStudio](https://www.rstudio.com/products/rstudio/download/) or a similar application of your choice. Please, make sure you have a reasonably fresh version of R itself. You will need the following R libraries: `ggtree`, `treeio`, `ape`, `ggplot2`, `geiger` (can be installed during the workshop).

### Course materials

1. [Environment setup on HPC cluster](00.Environment_setup.md)
2. [Scheduling jobs on HPC cluster](01.SLURM.md)
3. [Phylogenetic tree building on HPC cluster](02.Phylo_on_HPC.md)
4. [Visualization of phylogenetic trees: preparation](03.Tree_viz.md)
5. [Visualization of phylogenetic trees: annotation](04.Tree_plotting.md)

Visualization parts are also available as an R Script file: [`scripts_tree_operations.R`](https://raw.githubusercontent.com/Mycology-Microbiology-Center/Phylo2021/main/data/scripts_tree_operations.R)

### Example files
[Description of the input data](examples.md).<br/>
All example files could be found [here](https://github.com/Mycology-Microbiology-Center/Phylo2021/tree/main/data).

_Example of an ouput:_ <br/>
<img src="img/tree_example_high.png" width="700" title="Tree example"/>


### Contact
For any organizational questions contact [Iryna Yatsiuk](mailto:iryna.yatsiuk@ut.ee). <br/>
Workshop prepared by Vladimir Mikryukov, Anton Savchenko, and Iryna Yatsiuk.
