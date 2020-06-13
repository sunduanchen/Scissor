# Scissor: Single-Cell Identification of Subpopulations with bulk RNA-Seq phenotype cORrelation #

### Introduction ###
`Scissor` is a novel single cell data analysis approach. By leveraging bulk data and phenotype information, Scissor automatically selects a subset of cells from the single-cell data that is most responsible for the differences of phenotypes through a graph-regularized sparse regression model.

<p align="center">
<img src=Figure_Method.jpg height="602" width="500">
</p>

### News ###
Jun, 2020: Scissor version 1.0.0 is launched.

### Installation ###
* Prerequisites:
Scissor is developed under R (*version >= 3.6.1*). The [Seurat](https://satijalab.org/seurat/) packages (*version >= 3.0.0*) is used for loading data and preprocessing.


* Latest version: The latest developmental version of Scissor can be downloaded from GitHub and installed from source by
`# install.packages("devtools")`  
`devtools::install_github('sunduanchen/Scissor')`

### Manual ###
Please see https://sunduanchen.github.io/Scissor/vignettes/Scissor_Tutorial.html for details. In the R terminal, please use the command `?Scissor` to read the help documents.

### Example ###
Our tutorial https://sunduanchen.github.io/Scissor/vignettes/Scissor_Tutorial.html provided an example of Scissor, which identified a lung cancer cell subpopulation from lung cancer single cell data, guided by the TCGA-LUAD 471 bulk RNA-seq samples and their corresponding survival information.

### How to cite `Scissor` ###
Please cite the following publication:

> *Phenotype-guided subpopulation identification from single-cell sequencing data  
Duanchen Sun, Xiangnan Guan, Amy E Moran, David Z Qian, Pepper J Schedin, Andrew Adey, Paul T Spellman, and Zheng Xia*<br />


### License ###
Scissor is licensed under the GNU General Public License v3.0.
