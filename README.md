# Scissor: Single-Cell Identification of Subpopulations with bulk Sample phenOtype coRrelation #

### Introduction ###
`Scissor` is a novel approach that utilizes the phenotypes, such as disease stage, tumor metastasis, treatment response, and survival outcomes, collected from bulk assays to identify the most highly phenotype-associated cell subpopulations from single-cell data. The workflow of Scissor is shown in the following Figure:

<p align="center">
<img src=Figure_Method.jpg height="702" width="600">
</p>

### News ###
* May, 2021: Scissor version 2.1.0 is updated.  
    + Add utilities for cell level evaludations including correlation check and bootstrap (function: evaluate.cell)
* Feb, 2021: Scissor version 2.0.0 is launched.  
    + Optimize the inputs and outputs in Scissor main function
    + Add utilities for the reliability significance test (function: reliability.test)
* Jun, 2020: Scissor version 1.0.0 is launched.

### Installation ###
* Prerequisites:
Scissor is developed under R (*version >= 3.6.1*). The [Seurat](https://satijalab.org/seurat/) package (*version >= 3.2.0*) is used for loading data and preprocessing.

* Latest version: The latest developmental version of Scissor can be downloaded from GitHub and installed from source by
`devtools::install_github('sunduanchen/Scissor')`

### Manual ###
Please see https://sunduanchen.github.io/Scissor/vignettes/Scissor_Tutorial.html for details. In the R terminal, please use the command `?Scissor` to access the help documents.

### Examples ###
In our [Scissor Tutorial](https://sunduanchen.github.io/Scissor/vignettes/Scissor_Tutorial.html), we use several applications on the Lung Adenocarcinoma (LUAD) scRNA-seq cancer cells as examples to show how to execute Scissor in real applications.

### How to cite `Scissor` ###
Please cite the following manuscript:

> *Phenotype-guided subpopulation identification from single-cell sequencing data  
Duanchen Sun and Zheng Xia*<br />

### License ###
Scissor is licensed under the GNU General Public License v3.0.

Improvements and new features of Scissor will be updated on a regular basis. Please post on the [GitHub discussion page](https://github.com/sunduanchen/Scissor/discussions) with any questions.
