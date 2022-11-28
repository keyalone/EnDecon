# EnDecon
EnDecon integrates multiple base deconvolution results using a weighted optimization model to generate a more accurate result. EnDecon mainly includes two steps: (1) running each base deconvolution method individually to obtain the base cell type deconvolution results, and (2) integrating these base deconvolution results into a better deconvolution result using a new proposed ensemble strategy. EnDecon obtains the ensemble result by alternatively updating the ensemble result as a weighted median of the base deconvolution results and the weights of base results based on their distance from the ensemble result. R package applies ensemble learning for the deconvolution of spatial transcriptomic data. 
![alt
text](https://github.com/keyalone/EnDecon/blob/master/docs/Figure1.png?raw=true)

The EnDecon package has the main following R-package dependencies: SCDC, spacexr, MuSiC, DeconRNASeq, DWLS, Seurat, SPOTlight, Giotto, STdeconvolve, spatstat.geom, CARD, parallel, doParallel, foreach, reticulate and several python packages: scvi-tools, cell2location, scanpy, anndata. For the R-package dependencies, you can load on most of R dependencies packages on your R when install the EnDecon R package by run the code:
 ``` buildoutcfg
 devtools::install_github("Zhangxf-ccnu/EnDecon")
```
However, if the dependencies are not installed correctly, please install them by yourself by the following instruction. We check all the codes and examples in the package on the computer with the system of ubuntu 18.04, 64 GB RAM, i7-10700 CPU, RTX 3080 GPU. We also test the R package on the linux sever with A100 GPU. 

## Install individual dependencies
* **Install python dependencies**
 ``` buildoutcfg
 ### construct EnDecon python environment with pytorch GPU version 
 conda env create -f requirments_GPU.yml
 ### construct EnDecon python environment with pytorch CPU version
 conda env create -f requirment_EnDecon_CPU.yml
```
If you want to run the DWLS, SpatialDWLS, Stereoscope and cell2location for the ensemble learning, we advise that the user should install [anaconda]( https://www.anaconda.com/) and run the upper command on the terminal (ubuntu)/CMD (windows) to install the python dependencies for running the methods. In our application, due to the computer with RTX3080 GPU, we install the [pytorch with cudatookit]( https://pytorch.org/). If you don’t want to use the *.yml provided. You can install the python dependencies by the following code.
```buildoutcfg
 pip install scvi-tools
 pip install cell2location
 pip install scanpy
 pip install anadata
 pip install igraph
 pip install networkx
 pip install leidenalg
 pip install community
 pip install  smfishHmrf
 pip install scikit-learn
# install pytorch with CPU or GPU version
```
After install the python dependencies, the user need to get the path of environment of conda and set the path to the python_env variable in the function of EnDecon_individual_methods in our package. The path is similar to "\~/.conda/envs/EnDecon\_env/bin/python" on the ubuntu and "\~/anaconda3/envs/EnDecon\_env/python.ext" on Windows.
* **Install R dependencies**

SCDC
```buildoutcfg
install.packages("remotes")
  remotes::install_github("renozao/xbioc")
install.package("devtools")
  devtools::install_github("meichendong/SCDC")
```
RCTD
```buildoutcfg
devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
```
MuSiC
```buildoutcfg
devtools::install_github('xuranw/MuSiC')
```
DeconRNASeq
```buildoutcfg
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DeconRNASeq")
```
DWLS
```buildoutcfg
remotes::install_github("sistia01/DWLS")
```
Seurat
```buildoutcfg
install.packages("Seurat")
```
SPOTlight (Version 0.1.7)
```buildoutcfg
devtools::install_github("https://github.com/MarcElosua/SPOTlight/tree/spotlight-0.1.7")
```
Giotto
```buildoutcfg
devtools::install_github('RubD/Giotto')
```
spatstat.geom
```buildoutcfg
install.packages("spatstat.geom")
```
CARD
```buildoutcfg
devtools::install_github('YingMa0107/CARD')
```
STdeconvolve
```buildoutcfg
require(remotes)
remotes::install_github('JEFworks-Lab/STdeconvolve')
```
parallel and doParallel
```buildoutcfg
install.packages("parallel")
install.packages("doParallel")
```
reticulate
```buildoutcfg
install.packages('reticulate')
```
## Run the example
```buildoutcfg
data("breast.sc.ref")
data("breast.sc.cell.label")
data("breast.st")
data("breast.st.loc")
##### path on ubuntu platform on our computer
python_env <- "~/.conda/envs/EnDecon_GPU/bin/python"
Results.dec.mouse <- EnDecon_individual_methods(sc_exp = breast.sc.ref,
sc_label = breast.sc.cell.label, spot_exp = breast.st,
spot_loc = breast.st.loc, python_env = python_env,
use_gpu = TRUE,gene_det_in_min_cells_per = 0.01,
RCTD.CELL_MIN_INSTANCE = 5, saving_results = FALSE)
ensemble.results <- solve_ensemble(Results.dec.mouse[[1]])
```
## Recommendation for the selection of base deconvolution methods

![alt
text](https://github.com/keyalone/EnDecon/blob/master/docs/Figure2.png?raw=true)
![alt
text](https://github.com/keyalone/EnDecon/blob/master/docs/Figure4.png?raw=true)
![alt
text](https://github.com/keyalone/EnDecon/blob/master/docs/Figure3.png?raw=true)
![alt
text](https://github.com/keyalone/EnDecon/blob/master/docs/Figure5.png?raw=true)
For a computational method, the accuracy is important, but the running time also needs to be considered. Therefore, we also report the computational time requirement for the deconvolution methods. To obtain the running time, we run the deconvolution methods on a workstation with Intel core i7-10700 CPU (2.90GHz*16), 64 RAM and RTX 3080 GPU. We use the simulation data in Scenario 1 and Senirao2 and SRT with different downsampling rates of scRNA-seq data to show the accuracy and running time of deconvolution methods. The boxplots show that all individual deconvolution methods can be finished in less than 50 minutes. Cell2location, DestVI, DWLS and Stereoscope require more time than other methods. For the tables, we find that the running time of DestVI and Stereoscope is sensitive to sample size of reference scRNA-seq data. As the number of cells in scRNA-seq increases, the longer time requirement of the methods than other deconvolution methods. Note that after running the individual deconvolution methods, EnDecon can integrate the results from individual methods in a short time. Besides, we also provide the overview of the deconvolution methods in terms of PCC, 1-RMSE, 1-JSD and running time on the dataset for users to select appropriate deconvolution methods for ensemble learning.

## Tutorials
- [Analysis of Scenario 1 data with `EnDecon`](https://github.com/keyalone/EnDecon/blob/master/docs/simulation.md)

- [Analysis of breast cancer data with `EnDecon`](https://github.com/keyalone/EnDecon/blob/master/docs/EnDecon.md)

Please do not hesitate to contact Prof. Zhang at zhangxf@ccnu.edu.cn to seek any clarifications regarding any content or operation of the archive.

 
