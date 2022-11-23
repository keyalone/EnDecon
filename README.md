# EnDecon
EnDecon integrates multiple base deconvolution results using a weighted optimization model to generate a more accurate result. EnDecon mainly includes two steps: (1) running each base deconvolution method individually to obtain the base cell type deconvolution results, and (2) integrating these base deconvolution results into a better deconvolution result using a new proposed ensemble strategy. EnDecon obtains the ensemble result by alternatively updating the ensemble result as a weighted median of the base deconvolution results and the weights of base results based on their distance from the ensemble result. R package applies ensemble learning for the deconvolution of spatial transcriptomic data. 
![alt
text](https://github.com/keyalone/EnDecon/blob/master/docs/Figure1.png?raw=true)
![alt
text](https://github.com/keyalone/EnDecon/blob/master/docs/Figure2.png?raw=true)
![alt
text](https://github.com/keyalone/EnDecon/blob/master/docs/Figure3.png?raw=true)
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
devtools::install_github("MarcElosua/SPOTlight/tree/spotlight-0.1.7")
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
'''buildoutcfg
devtools::install_github('YingMa0107/CARD')
'''
parallel and doParallel
'''buildoutcfg
install.packages("parallel")
install.packages("doParallel")
'''
reticulate
```buildoutcfg
Install.packages('reticulate')
```
## Run the example
```buildoutcfg
data("MVC.reference")
data("MVC.reference.cell.label")
data("MVC.ST")
data("MVC.ST.coor")
##### path on ubuntu platform on our computer
python_env <- "~/.conda/envs/EnDecon_GPU/bin/python"
##### we use all the genes for the deconvolutioin
Results.Deconv <- EnDecon_individual_methods(MVC.reference, MVC.reference.cell.label,
                  MVC.ST, MVC.ST.coor,  python_env = python_env, use_gpu = TRUE,
                  RCTD.CELL_MIN_INSTANCE = 10, gene_det_in_min_cells_per = 0,
                  expression_threshold = 0, nUMI = 1, DWLS.is_select_DEGs = FALSE,
                  SpatialDWLS.is_select_DEGs = FALSE)
ensemble.results <- solve_ensemble(Results.Deconv)
```
Please do not hesitate to contact Prof. Zhang at zhangxf@mails.ccnu.edu.cn to seek any clarifications regarding any content or operation of the archive.

 
