# Scissor: Single-Cell Identification of Subpopulations with bulk Sample phenOtype coRrelation #

### Introduction ###
`Scissor` is a novel approach that utilizes the phenotypes, such as disease stage, tumor metastasis, treatment response, and survival outcomes, collected from bulk assays to identify the most highly phenotype-associated cell subpopulations from single-cell data. The workflow of Scissor is shown in the following Figure:

<p align="center">
<img src=Figure_Method.jpg height="702" width="600">
</p>

### News ###
Jun, 2020: Scissor version 1.0.0 is launched.
Apr, 2021: Scissor version 2.0.0 is launched.

### Installation ###
* Prerequisites:
Scissor is developed under R (*version >= 3.6.1*). The [Seurat](https://satijalab.org/seurat/) package (*version >= 3.2.0*) is used for loading data and preprocessing.

* Latest version: The latest developmental version of Scissor can be downloaded from GitHub and installed from source by
`devtools::install_github('sunduanchen/Scissor')`

### Manual ###
Please see https://sunduanchen.github.io/Scissor/vignettes/Scissor_Tutorial.html for details. In the R terminal, please use the command `?Scissor` to access the help documents.

### Example ###
Our tutorial https://sunduanchen.github.io/Scissor/vignettes/Scissor_Tutorial.html provides  examples of Scissor, which identified a lung cancer cell subpopulation from lung cancer single cell data, guided by the TCGA-LUAD 471 bulk RNA-seq samples and their corresponding survival information.

### How to cite `Scissor` ###
Please cite the following manuscript:

> *Phenotype-guided subpopulation identification from single-cell sequencing data  
Duanchen Sun and Zheng Xia*<br />

### License ###
Scissor is licensed under the GNU General Public License v3.0.
