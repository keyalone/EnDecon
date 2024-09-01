#' This function focuses on cleaning scRNA-seq and stRNA-seq datasets.
#'
#'
#' @param sc_exp scRNA-seq matrix, genes * cells. The format should be raw-counts. The matrix need include gene names and cell names.
#' @param sc_label cell type information. The cells are need be divided into multiple category.
#' @param spot_exp stRNA-seq matrix, genes * spots. The format should be raw counts. The matrix need include gene names and spot names.
#' @param spot_loc coordinate matrix, spots * coordinates. The matrix need include spot names and coordinate name (x, y).
#' @param gene_det_in_min_cells_per a floor variable. minimum percent of genes that need to be detected in a cell.
#' @param expression_threshold a floor variable. Threshold to consider a gene expressed.
#' @param nUMI a floor variable. 	minimum of read count that need to be detected in a cell or spot.
#' @param verbose a logical variable that defines whether to print the processing flow of data process.
#' @param plot a logical variable that defines whether to plot the selected genes and selected cell expression.
#'
#' @return a list includes processed scRNA-seq matrix, cell type, stRNA-seq matrix.
#'
#' @examples
#' data("breast.sc.ref")
#' data("breast.sc.cell.label")
#' data("breast.st")
#' data("breast.st.loc")
#' database <- data_process(breast.sc.ref, breast.sc.cell.label, breast.st, breast.st.loc)
#'
#' @export

data_process <- function(sc_exp, sc_label, spot_exp, spot_loc,
                         gene_det_in_min_cells_per = 0.01, expression_threshold = 1,
                         nUMI = 100, verbose = FALSE, plot= FALSE){
  if(ncol(sc_exp) != length(sc_label))
    stop("Require cell labels!")

  if(ncol(spot_exp) != nrow(spot_loc))
    stop("Require x , y coordinations")

  #### scRNA-seq data process
  sc_matrix = t(cleanCounts(t(sc_exp), gene_det_in_min_cells_per = gene_det_in_min_cells_per,
                            expression_threshold = expression_threshold, nUMI = nUMI,
                            verbose = verbose, plot = plot))
  sc_matrix= as.matrix(sc_matrix)
  ind = match(colnames(sc_matrix), colnames(sc_exp))
  sc_label = sc_label[ind]
  # cell_type = sort(unique(sc_label))

  #### ST data process
  st_matrix = t(cleanCounts(t(spot_exp), gene_det_in_min_cells_per = gene_det_in_min_cells_per,
                            expression_threshold = expression_threshold, nUMI = nUMI,
                            verbose = verbose, plot = plot))
  st_matrix= as.matrix(st_matrix)
  ind_sp = match(colnames(st_matrix), colnames(spot_exp))
  spot_loc = spot_loc[ind_sp, ]

  #### find common genes
  com_gene = intersect(rownames(sc_matrix),rownames(st_matrix))
  sc_exp = sc_matrix[com_gene,]
  st_exp = st_matrix[com_gene,]

  ### rechecking nUMI
  index_sc <- colSums(sc_exp) >= nUMI
  sc_exp_filter <- sc_exp[,index_sc]
  sc_label_filter <- sc_label[index_sc]

  index_st <- colSums(st_exp) >= nUMI
  st_exp_filter = st_exp[,index_st]
  spot_loc_filter <- spot_loc[index_st,]

  database <- list(sc_exp = sc_exp_filter, sc_label = sc_label_filter,
                   spot_exp = st_exp_filter, spot_loc = spot_loc_filter)
  return(database)
}

cleanCounts <- function (counts, gene_det_in_min_cells_per = 0.01,
                         expression_threshold = 1 ,
                         nUMI = 100,
                         verbose = FALSE, plot= FALSE) {
  n = nrow(counts)
  ##### select of the genes
  filter_index_genes = Matrix::colSums(counts >= expression_threshold) >=
    gene_det_in_min_cells_per*n

  #### filter the cell
  filter_index_cells = Matrix::rowSums(counts[,filter_index_genes] >=
                                         expression_threshold) >= nUMI

  counts = counts[filter_index_cells, filter_index_genes]

  if (verbose) {
    message("Resulting matrix has ", nrow(counts), " cells and ", ncol(counts), " genes")
  }
  if (plot) {
    par(mfrow=c(1,2), mar=rep(5,4))
    hist(log10(Matrix::rowSums(counts)+1), breaks=20, main='Genes Per Dataset')
    hist(log10(Matrix::colSums(counts)+1), breaks=20, main='Datasets Per Gene')
  }
  return(counts)
}

####### Utils of data process of DWLS ##########
create_group_exp <- function(sc_exp,sc_label) {

  #sc_exp single cell gene expression datasets
  #sc_label  cell annotation of the single cells of the reference

  ##group cells
  # reference matrix (C) + refProfiles.var from TRAINING dataset
  cell_type = sort(unique(sc_label))
  group = list()
  for(i in 1:length(cell_type)){
    temp_use <- which(sc_label == cell_type[i])
    names(temp_use) <- NULL
    group[[i]] <- temp_use
  }
  sc_group_exp = sapply(group,function(x) Matrix::rowMeans(sc_exp[,x]))
  #sapply
  sc_group_exp = as.matrix(sc_group_exp)
  colnames(sc_group_exp) = cell_type
  return(sc_group_exp)
}


############## apply SCDC #####################
SCDC_run = function(database, iter.max = 1000){

  ### extract SC, CL, ST from database
  sc_exp <- database$sc_exp
  sc_label <- database$sc_label
  spot_exp <- database$spot_exp
  cell_type <- sort(unique(sc_label))

  sc_metadata <- data.frame(row.names = c('cluster', 'sample'))
  sc_pData <- data.frame(as.character(sc_label), as.character(colnames(sc_exp)))
  colnames(sc_pData) <- c('cluster', 'sample')
  rownames(sc_pData) <- colnames(sc_exp)
  cat(paste0('SCDC scref num: ',length(rownames(sc_pData))))
  sc_phenoData <- methods::new("AnnotatedDataFrame", data=sc_pData, varMetadata=sc_metadata)
  sc_exp_es <- Biobase::ExpressionSet(as.matrix(sc_exp), phenoData=sc_phenoData, byrow=FALSE)

  st_metadata <- data.frame(row.names=c('sample'))
  st_pData <- data.frame(as.character(colnames(spot_exp)))
  colnames(st_pData) <- c('sample')
  rownames(st_pData) <- colnames(spot_exp)
  st_phenoData <- methods::new("AnnotatedDataFrame", data=st_pData,varMetadata=st_metadata)
  spot_matrix_es <- Biobase::ExpressionSet(as.matrix(spot_exp), phenoData=st_phenoData, byrow=FALSE)

  bulkseger.scseger <- suppressMessages(SCDC::SCDC_prop(bulk.eset = spot_matrix_es,
                                                        sc.eset = sc_exp_es,
                                                        ct.varname = "cluster",
                                                        ct.sub = cell_type,
                                                        iter.max = iter.max))

  SCDC_results <- bulkseger.scseger$prop.est.mvw
  SCDC_results <- SCDC_results[,cell_type]
  return(SCDC_results)
}


############## apply RCTD #####################
RCTD_run = function(database, CELL_MIN_INSTANCE = 20){
  ### extract SC, CL, ST from database
  sc_exp <- database$sc_exp
  sc_label <- database$sc_label
  spot_exp <- database$spot_exp
  cell_type <- sort(unique(sc_label))
  spot_loc <- database$spot_loc

  sparse_sc_exp <- as(sc_exp, "sparseMatrix")
  sparse_spot_exp <- as(spot_exp, "sparseMatrix")

  ## the reference scRNA-seq data
  cellnames <- colnames(sc_exp)
  cell_types <- as.factor(sc_label)
  names(cell_types) <- cellnames
  sc_nUMI <- as.numeric(colSums(sc_exp))
  names(sc_nUMI) <- cellnames
  reference <- spacexr::Reference(sparse_sc_exp, cell_types, nUMI = sc_nUMI)

  ### Create SpatialRNA object
  coords <- as.data.frame(spot_loc)
  # coords <- as.data.frame(matrix(1,dim(spot_exp)[2],2))
  rownames(coords) <- as.character(colnames(spot_exp))
  nUMI <- colSums(spot_exp)
  puck <- spacexr::SpatialRNA(coords, counts=sparse_spot_exp, nUMI=nUMI)

  myRCTD <- spacexr::create.RCTD(puck, reference, max_cores = 1,CELL_MIN_INSTANCE = CELL_MIN_INSTANCE, UMI_min=5, counts_MIN = 1)
  myRCTD <- spacexr::run.RCTD(myRCTD, doublet_mode = 'full')
  results <- myRCTD@results

  temp <- as.matrix(results$weights)
  norm_weights_temp <- sweep(temp, 1, rowSums(temp), '/')
  RCTD_results <- norm_weights_temp[,cell_type]

  return(RCTD_results)
}


############## apply MuSiC #####################
MuSiC_run = function(database, iter.max = 1000, nu = 1e-04, eps = 0.01){
  ### extract SC, CL, ST from database
  sc_exp <- database$sc_exp
  sc_label <- database$sc_label
  spot_exp <- database$spot_exp
  cell_type <- sort(unique(sc_label))


  sc_label_num <- sc_label
  for (k in seq(length(unique(sc_label)))) {
    sc_label_num[sc_label == sort(unique(sc_label))[k]] <- k
  }

  # sc_metadata <- data.frame(labelDescription=c('Sample ID','Subject Name','Cell Type ID ', 'Cell Type Name'),
  #                           row.names=c('sampleID','SubjectName','cellTypeID ','cellType'))


  sc_pData <- data.frame(as.numeric(c(1:dim(sc_exp)[2])), as.factor(colnames(sc_exp)),
                         as.numeric(sc_label_num), as.factor(sc_label))
  colnames(sc_pData) <- c('sampleID', 'SubjectName', 'cellTypeID ', 'cellType')
  rownames(sc_pData) <- colnames(sc_exp)
  # sc_phenoData <- methods::new("AnnotatedDataFrame", data = sc_pData,
  #                              varMetadata=sc_metadata)
  # sc_exp_es <- Biobase::ExpressionSet(as.matrix(sc_exp),
  #                                     phenoData=sc_phenoData, byrow=FALSE)

  # normcounts <- normalized(sc_exp, method = norm.form)
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = sc_exp))
  gene_df <- data.frame(gene.name = rownames(sce))
  # cell_df <- data.frame(label = sc_label, cell = colnames(sce))
  rownames(gene_df) <- gene_df$Gene
  # rownames(cell_df) <- cell_df$cell
  cat(paste0('MuSiC scref num: ', dim(sc_exp)))
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = sc_exp),
                                                      colData = sc_pData,
                                                      rowData = gene_df)


  st_metadata <- data.frame(labelDescription=c('Sample ID','Subject Name'), row.names=c('sampleID','SubjectName'))
  st_pData <- data.frame(as.integer(c(1:dim(spot_exp)[2])), as.factor(colnames(spot_exp)))
  colnames(st_pData) <- c('sampleID', 'SubjectName')
  rownames(st_pData) <- colnames(spot_exp)
  st_phenoData <- methods::new("AnnotatedDataFrame", data=st_pData,varMetadata=st_metadata)
  spot_exp_es <- Biobase::ExpressionSet(as.matrix(spot_exp),phenoData= st_phenoData, byrow=FALSE)
  bulk.mtx = exprs(spot_exp_es)

  #### MuSiC version = 0.2.0
  Est.prop.GSE = MuSiC::music_prop(bulk.mtx = bulk.mtx, sc.sce = sce,
                                   clusters = 'cellType',
                                   samples = 'sampleID', select.ct = cell_type, verbose = F,
                                   iter.max = iter.max, nu = nu, eps = eps)

  Music_weight = as.matrix(Est.prop.GSE$Est.prop.weighted)
  Music_allgene = as.matrix(Est.prop.GSE$Est.prop.allgene)

  Music_weight <- Music_weight[,cell_type]
  Music_allgene <- Music_allgene[,cell_type]

  Music_results = list(Music_weight = Music_weight, Music_allgene = Music_allgene)

  return(Music_results)
}


############## apply DeconRNASeq #####################
# perc variable use to filter genes
DeconRNASeq_run = function(database, perc = 0.05){
  ### extract SC, CL, ST from database
  sc_exp <- database$sc_exp
  sc_label <- database$sc_label
  spot_exp <- database$spot_exp
  cell_type <- sort(unique(sc_label))


  # perc <- 0.05
  n <- dim(sc_exp)[2]
  id <- rowSums(sc_exp) >= perc * n
  sc_exp <- sc_exp[id,]

  sc_exp <- t(median(colSums(sc_exp))*(t(sc_exp)/colSums(sc_exp)))
  sc_exp <- log(sc_exp + 1)

  spot_exp <- t(median(colSums(spot_exp))*(t(spot_exp)/colSums(spot_exp)))
  spot_exp <- log(spot_exp+1)

  type_exp <- NULL
  for (k in seq(length(cell_type))) {
    exp <- rowSums(sc_exp[,sc_label == cell_type[k]])/sum(sc_label == cell_type[k])
    type_exp <- cbind(type_exp,exp)
  }
  type_exp <- as.matrix(type_exp)
  colnames(type_exp) <- cell_type

  decon_out <- DeconRNASeq::DeconRNASeq(as.data.frame(spot_exp), as.data.frame(type_exp), checksig=FALSE, use.scale = TRUE)

  DeconRNASeq_results <- decon_out$out.all
  rownames(DeconRNASeq_results) <- colnames(spot_exp)
  DeconRNASeq_results <- DeconRNASeq_results[,cell_type]

  return(DeconRNASeq_results)
}


########## run DWLS and SVR ##############
DWLS_run = function(database, parallel = TRUE,  is_select_DEGs = TRUE, python_env){

  ### extract SC, CL, ST from database
  sc_exp <- database$sc_exp
  sc_label <- database$sc_label
  spot_exp <- database$spot_exp
  cell_type <- sort(unique(sc_label))


  if(is_select_DEGs){
    instrs = Giotto::createGiottoInstructions(python_path = python_env)

    sc_cortex <- Giotto::createGiottoObject(raw_exprs = sc_exp,
                                            instructions = instrs)

    sc_cortex <- Giotto::normalizeGiotto(gobject = sc_cortex)

    sc_cortex@cell_metadata$leiden_clus <- as.character(sc_label)

    gini_markers_subclusters = Giotto::findMarkers_one_vs_all(gobject = sc_cortex,
                                                              method = 'gini',
                                                              expression_values = 'normalized',
                                                              cluster_column = 'leiden_clus',
                                                              verbose = FALSE)
    topgenes_mast = gini_markers_subclusters[, head(.SD, 100), by = 'cluster']

    sc_exp <- sc_exp[topgenes_mast$genes, ]
    spot_exp <- spot_exp[topgenes_mast$genes, ]
  }
  cellAnnots <- data.frame(CellID = colnames(sc_exp),
                           cellType = sc_label)

  cell_type_exp <- create_group_exp(sc_exp,sc_label)
  cat(paste0('DWLS cell_type_exp dim : ', dim(cell_type_exp)))

  if(parallel){
    #### because the cores can't large than 2 in using R CMD check processing
    chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    if (nzchar(chk) && chk == "TRUE") {
      # use 2 cores in CRAN/Travis/AppVeyor
      num.cores <- 2L
    } else {
      # use all cores in devtools::test()
      num.cores <- parallel::detectCores() - 1
    }
    myCluster <- parallel::makeCluster(num.cores)
    doParallel::registerDoParallel(myCluster)

    pred.DWLS <- Sys.time()
    DampenedWLS <- foreach::foreach(i=1:ncol(spot_exp), .combine = 'rbind', .inorder = TRUE) %dopar% {
      DWLS::solveDampenedWLS(cell_type_exp, spot_exp[,i])
    }
    end.DWLS <- Sys.time()
    time.DWLS <- difftime(end.DWLS, pred.DWLS, units = "mins")

    pred.SVR <- Sys.time()
    SVR <- foreach::foreach(i=1:ncol(spot_exp), .combine = 'rbind', .inorder = TRUE) %dopar% {
      cat(i)
      DWLS::solveSVR(cell_type_exp, spot_exp[,i])
    }
    end.SVR <- Sys.time()
    time.SVR <- difftime(end.SVR, pred.SVR, units = "mins")

    parallel::stopCluster(myCluster)

  }else{
    SVR = DampenedWLS = matrix(NA, ncol(spot_exp), length(cell_type))

    pred.DWLS <- Sys.time()
    for (i in seq(ncol(spot_exp))) {
      DampenedWLS[i,] <- DWLS::solveDampenedWLS(cell_type_exp, spot_exp[,i])
    }
    end.DWLS <- Sys.time()
    time.DWLS <- difftime(end.DWLS, pred.DWLS, units = "mins")

    pred.SVR <- Sys.time()
    for (i in seq(ncol(spot_exp))) {
      SVR[i,] <- DWLS::solveSVR(cell_type_exp, spot_exp[,i])
    }
    end.SVR <- Sys.time()
    time.SVR <- difftime(end.SVR, pred.SVR, units = "mins")
  }

  rownames(DampenedWLS) <- colnames(spot_exp)
  rownames(SVR) <- colnames(spot_exp)

  colnames(DampenedWLS) <- cell_type
  colnames(SVR) <- cell_type

  DWLS <- list(DampenedWLS = DampenedWLS, SVR = SVR, time.DWLS = time.DWLS, time.SVR = time.SVR)
  return(DWLS)
}



############## apply SPOTlight #####################
SPOTlight_run = function(database, cl_n = 100, hvg = 1000, min_cont = 0.001){
  ### extract SC, CL, ST from database
  sc_exp <- database$sc_exp
  sc_label <- database$sc_label
  spot_exp <- database$spot_exp
  cell_type <- sort(unique(sc_label))


  sc_ref <- Seurat::CreateSeuratObject(counts = sc_exp)
  sc_ref@meta.data$subclass <- as.factor(sc_label)

  sc_ref <- Seurat::SCTransform(sc_ref, verbose = FALSE)

  Seurat::Idents(object = sc_ref) <- sc_ref@meta.data$subclass
  cluster_markers_all <- Seurat::FindAllMarkers(object = sc_ref,
                                                assay = "SCT",
                                                slot = "data",
                                                verbose = FALSE,
                                                only.pos = TRUE)

  # Downsample scRNAseq to select gene set and number of cells to train the model
  # with default parameters
  se_sc_down <- SPOTlight::downsample_se_obj(se_obj = sc_ref,
                                             clust_vr = "subclass",
                                             cluster_markers = cluster_markers_all,
                                             cl_n = cl_n,
                                             hvg = hvg)

  nmf_mod_ls <- SPOTlight::train_nmf(cluster_markers = cluster_markers_all,
                                     se_sc = se_sc_down,
                                     mtrx_spatial = spot_exp,
                                     clust_vr = "subclass",
                                     ntop = NULL,
                                     hvg = hvg,
                                     transf = "uv",
                                     method = "nsNMF")

  nmf_mod <- nmf_mod_ls[[1]]

  # get basis matrix W
  w <- NMF::basis(nmf_mod)

  # get coefficient matrix H
  h <- NMF::coef(nmf_mod)

  # Subset to genes used to train the model
  temp_index <- ! is.na(pmatch(rownames(spot_exp), rownames(w)))
  spot_exp_train <- spot_exp[temp_index, ]

  ct_topic_profiles <- SPOTlight::topic_profile_per_cluster_nmf(h = h,
                                                                train_cell_clust = nmf_mod_ls[[2]])

  decon_mtrx <- SPOTlight::mixture_deconvolution_nmf(nmf_mod = nmf_mod,
                                                     mixture_transcriptome = spot_exp_train,
                                                     transf = "uv",
                                                     reference_profiles = ct_topic_profiles,
                                                     min_cont = min_cont)

  SPOTlight_results <- as.matrix(decon_mtrx[,- ncol(decon_mtrx)])
  rownames(SPOTlight_results) <- colnames(spot_exp)
  colnames(SPOTlight_results) <- cell_type

  return(SPOTlight_results)
}



############## apply SpatialDWLS #####################
spatialDWLS_run = function(database, my_python_path, is_select_DEGs = FALSE,
                           findmarker_method = "gini",
                           ncp_spa = 100, dimensions_to_use = 10, k = 10,
                           resolution = 0.4, n_iterations = 1000, n_cell = 50){
  ### extract SC, CL, ST from database
  sc_exp <- database$sc_exp
  sc_label <- database$sc_label
  spot_exp <- database$spot_exp
  cell_type <- sort(unique(sc_label))

  instrs = Giotto::createGiottoInstructions(python_path = my_python_path)

  sc_cortex <- Giotto::createGiottoObject(raw_exprs = sc_exp,
                                          instructions = instrs)

  sc_cortex <- Giotto::normalizeGiotto(gobject = sc_cortex)

  sc_cortex@cell_metadata$leiden_clus <- as.character(sc_label)

  if(is_select_DEGs){
    gini_markers_subclusters = Giotto::findMarkers_one_vs_all(gobject = sc_cortex,
                                                              method = findmarker_method,
                                                              expression_values = 'normalized',
                                                              cluster_column = 'leiden_clus',
                                                              verbose = FALSE)
    topgenes_gini = gini_markers_subclusters[, head(.SD, 100), by = 'cluster']
    sc_norm_exp <- 2^(sc_cortex@norm_expr)-1
    ExprSubset <- sc_norm_exp[as.character(topgenes_gini$genes),]
  }else{
    sc_norm_exp <- 2^(sc_cortex@norm_expr)-1
    ExprSubset <- sc_norm_exp
  }


  Sig<-NULL
  for (i in as.character(unique(sc_label))){
    Sig<-cbind(Sig,(apply(ExprSubset,1,function(y) mean(y[which(sc_label==i)]))))
  }
  colnames(Sig)<-as.character(unique(sc_label))



  grid_seqFish <- Giotto::createGiottoObject(raw_exprs = spot_exp,instructions = instrs)
  grid_seqFish <- Giotto::normalizeGiotto(gobject = grid_seqFish)
  grid_seqFish <- Giotto::calculateHVG(gobject = grid_seqFish,show_plot = FALSE,
                                       return_plot =  FALSE)
  gene_metadata = Giotto::fDataDT(grid_seqFish)
  featgenes = gene_metadata[hvg == 'yes']$gene_ID
  grid_seqFish <- Giotto::runPCA(gobject = grid_seqFish, genes_to_use = featgenes, scale_unit = F, ncp = ncp_spa)
  grid_seqFish <- Giotto::createNearestNetwork(gobject = grid_seqFish, dimensions_to_use = 1:dimensions_to_use, k = k)
  grid_seqFish <- Giotto::doLeidenCluster(gobject = grid_seqFish, resolution = resolution, n_iterations = n_iterations)

  grid_seqFish <- Giotto::runDWLSDeconv(gobject = grid_seqFish, sign_matrix = Sig, n_cell = n_cell)

  spatialDWLS_results <- as.matrix(grid_seqFish@spatial_enrichment$DWLS[,-1])
  spatialDWLS_results <- spatialDWLS_results[,cell_type]
  rownames(spatialDWLS_results) <- colnames(spot_exp)

  return(spatialDWLS_results)
}

############## apply Stereoscope #####################


Stereoscope_run <- function(database, python_env = NULL, use_gpu = FALSE, select_HVG = TRUE, HVG_num = 2000,
                            sc_training_plot = FALSE, sc_training_save_trained_model = FALSE, sc_max_epochs = 10000, sc_lr = 0.01,
                            st_training_plot = FALSE, st_training_save_trained_model = FALSE, st_max_epochs = 10000, st_lr = 0.01){

  if(is.null(python_env)){
    cat("if you want to use the environment of your own ptython, please set the path of your python environment
        in yout computer. Or the code will install Miniconda and scvi-tools to run Stereoscope. We recommond
        the user to construct anaconda env with the requirments.yaml by run the code <conda env create -f requirments.yaml> we provided.")
  }else{
    reticulate::use_python(python_env, require = T)
    reticulate::py_config()
  }

  ### extract SC, CL, ST from database
  sc_exp <- database$sc_exp
  sc_label <- database$sc_label
  spot_exp <- database$spot_exp
  cell_type <- sort(unique(sc_label))
  # st_expr <- spot_exp
  cell_label_sc <- cell_type


  if(nrow(sc_exp) != nrow(spot_exp))
    stop("The features of scRNA-seq and stRNA-seq need be equal!")
  # pip install scvi-tools
  scvi <- reticulate::import("scvi")
  np <- reticulate::import("numpy")
  sc <- reticulate::import("scanpy")
  # pip install anndata
  anndata <- reticulate::import("anndata")
  # plt <- reticulate::import("matplotlib.pyplot")
  #### Stereoscppe function in scvi
  #### scvi.external.stereoscope.RNAStereoscope
  #### scvi.external.stereoscope.SpatialStereoscope
  #### scvi.external.stereoscope.RNAStereoscope.setup_anndata
  #### scvi.external.stereoscope.SpatialStereoscope.setup_anndata
  ### training for scRNA-seq
  # dataprocess
  # adata_sc <- anndata$AnnData(X = t(sc_exp))
  if((nrow(sc_exp) >= HVG_num)){
    if(select_HVG){
      # adata_sc <- anndata$AnnData(X = t(sc_exp))
      # sc$pp$normalize_total(adata_sc, target_sum = as.integer(1e5))
      # sc$pp$log1p(adata_sc)
      # sc$pp$highly_variable_genes(adata_sc, n_top_genes = as.integer(HVG_num), flavor = "seurat_v3")
      Seurat.sc <- Seurat::CreateSeuratObject(counts = sc_exp)
      Seurat.sc <- Seurat::NormalizeData(object = Seurat.sc, verbose = FALSE)
      Seurat.sc <- Seurat::FindVariableFeatures(object = Seurat.sc, nfeatures = HVG_num, verbose = FALSE)
      var.features <- Seurat.sc@assays$RNA@var.features

      sc_exp_hvg <- sc_exp[var.features,]
      spot_exp_hvg <- spot_exp[var.features,]
      adata_sc <- anndata$AnnData(X = t(sc_exp_hvg))
      st_expr <- spot_exp_hvg
    }else{
      message("We advise to select HVG by setting the select_HVG = TRUE")
      adata_sc <- anndata$AnnData(X = t(sc_exp))
      st_expr <- spot_exp
    }
  }else{
    adata_sc <- anndata$AnnData(X = t(sc_exp))
    st_expr <- spot_exp
  }


  adata_sc$obs["cell_label"] <- sc_label
  scvi$external$stereoscope$RNAStereoscope$setup_anndata(adata_sc, labels_key =  "cell_label")
  # training setting
  sc_model <- scvi$external$stereoscope$RNAStereoscope(adata_sc)
  ## advised by cell2location: using the default of Stereoscope.
  scvi$external$stereoscope$RNAStereoscope$train(sc_model, use_gpu = use_gpu, check_val_every_n_epoch = as.integer(1),
                                                 early_stopping = FALSE, early_stopping_monitor = "elbo_validation",
                                                 max_epochs = as.integer(sc_max_epochs), lr = sc_lr)

  if(sc_training_plot){
    loss_r <- np$asarray(sc_model$history["elbo_train"]$elbo_train)
    loss_plot <- data.frame(epoches = 40:nrow(loss_r), elbo_train = unlist(loss_r[40:nrow(loss_r)]))
    plot(loss_plot, type = "l", main = "scRNA-seq training")
  }
  ###### You can reload the trained model by
  ##### scvi$external$stereoscope$RNAStereoscope$load("scmodel", adata_sc)
  if(sc_training_save_trained_model){
    scvi$external$stereoscope$RNAStereoscope$save(sc_model, "scmodel")
  }

  ######### training for sc_data
  adata_st <- anndata$AnnData(X = t(st_expr))
  scvi$external$stereoscope$SpatialStereoscope$setup_anndata(adata_st)
  st_model <- scvi$external$stereoscope$SpatialStereoscope$from_rna_model(adata_st, sc_model, prior_weight = "minibatch")

  scvi$external$stereoscope$SpatialStereoscope$train(st_model, max_epochs = as.integer(st_max_epochs), check_val_every_n_epoch = as.integer(1),
                                                     early_stopping_monitor = "elbo_train", use_gpu = use_gpu,
                                                     max_epochs = as.integer(st_max_epochs), lr = st_lr)

  if(st_training_plot){
    loss_r <- np$asarray(st_model$history["elbo_train"]$elbo_train)
    loss_plot <- data.frame(epoches = 40:nrow(loss_r), elbo_train = unlist(loss_r[40:nrow(loss_r)]))
    plot(loss_plot, type = "l", main = "stRNA-seq training")
  }

  ### get the prob
  weight_label <- scvi$external$stereoscope$SpatialStereoscope$get_proportions(st_model)
  weight_label <- as.matrix(weight_label)
  rownames(weight_label) <- colnames(st_expr)
  weight_label <- weight_label[,cell_type]
  return(weight_label)
}


############## apply cell2location #####################

cell2location_run <- function(database, python_env = NULL,
                              sc_max_epoches = 3000, sc_lr = 0.002, sc_use_gpu = FALSE,
                              st_N_cells_per_location = 30, st_detection_alpha = 200.00,
                              st_max_epoches = 3000){

  if(is.null(python_env)){
    cat("if you want to use the environment of your own ptython, please set the path of your python environment
        in yout computer. Or the code will install Miniconda and scvi-tools, cell2location to run cell2location.  We recommond
        the user to construct anaconda env with the requirments.yaml by run the code <conda env create -f requirments.yaml> we provided.")
  }else{
    reticulate::use_python(python_env, require = T)
    reticulate::py_config()
  }
  ### extract SC, CL, ST from database
  sc_exp <- database$sc_exp
  sc_label <- database$sc_label
  spot_exp <- database$spot_exp
  cell_type <- sort(unique(sc_label))


  st_expr <- spot_exp
  cell_label_sc <- cell_type


  if(nrow(sc_exp) != nrow(st_expr))
    stop("The features of scRNA-seq and stRNA-seq need be equal!")

  scvi <- reticulate::import("scvi")
  np <- reticulate::import("numpy")
  anndata <- reticulate::import("anndata")
  # pip install cell2location
  cell2loc <- reticulate::import("cell2location")
  pd <- reticulate::import("pandas")

  ############## training for scRNA-seq


  adata_sc <- input_for_py(expr = sc_exp, cell_type = sc_label)
  cell2loc$models$RegressionModel$setup_anndata(adata_sc, labels_key = "cell_type")
  mod <- cell2loc$models$RegressionModel(adata_sc)
  
  cat(paste0('ncol(sc_exp)=', ncol(sc_exp)))
  if(ncol(sc_exp) > 10000){
    batch_size <- as.integer(2500)
  }else{
    batch_size <- as.integer(floor(ncol(sc_exp)/2))
  }

  cell2loc$models$RegressionModel$train(mod, max_epochs = as.integer(sc_max_epoches), batch_size = as.integer(batch_size), lr = sc_lr, use_gpu = sc_use_gpu)

  #### extract estimated cell abundance
  parameters_setting <- list(batch_size = as.integer(batch_size), use_gpu = sc_use_gpu)
  adata_sc <- cell2loc$models$RegressionModel$export_posterior(mod, adata_sc,
                                                               sample_kwargs = parameters_setting)
  inf_aver <- data_frame_extracted(adata_sc)


  ############# training for stRNA-seq
  adata_st <- input_for_py(expr = st_expr)

  cell2loc$models$Cell2location$setup_anndata(adata_st)
  mod_st <- cell2loc$models$Cell2location(adata_st, cell_state_df = inf_aver, N_cells_per_location = as.integer(st_N_cells_per_location),
                                          detection_alpha = st_detection_alpha)

  cell2loc$models$Cell2location$train(mod_st, max_epochs = as.integer(st_max_epoches), train_size = 1.00, use_gpu = sc_use_gpu)

  parameters_setting_st <- list(batch_size = as.integer(ncol(st_expr)), use_gpu = sc_use_gpu)
  adata_st <- cell2loc$models$Cell2location$export_posterior(mod_st, adata_st, sample_kwargs = parameters_setting_st)

  ### extracted high confident cell abundance
  adata_temp <- reticulate::py_to_r(adata_st)
  data_frame_total <- adata_temp$obsm["q05_cell_abundance_w_sf"]
  # slice_idx_tmp <- adata_temp$uns['mod']['factor_names']$factor_names
  slice_idx_tmp <- adata_temp$uns$mod$factor_names
  slice_idx <- paste0("q05cell_abundance_w_sf_", slice_idx_tmp)
  data_frame_extracted <- data_frame_total[slice_idx]
  names(data_frame_extracted) <- slice_idx_tmp

  weight_matrix <- matrix(unlist(data_frame_extracted), ncol = length(data_frame_extracted))
  colnames(weight_matrix) <- slice_idx_tmp
  rownames(weight_matrix) <- rownames(data_frame_extracted)
  weight_matrix <- weight_matrix[,cell_type]

  weight_matrix_norm <- sweep(weight_matrix, 1, rowSums(weight_matrix), "/")
  return(weight_matrix_norm)
}

####### apply STdeconvolve #########
### opt_method : kneed, min
STdeconvolve_run <- function(database, min.lib.size = 100, min.reads = 1, nTopOD = 1000, betaScale = 1000){
  ### extract SC, CL, ST from database
  sc_exp <- database$sc_exp
  sc_label <- database$sc_label
  spot_exp <- database$spot_exp
  cell_type <- sort(unique(sc_label))

  counts <- STdeconvolve::cleanCounts(counts = spot_exp, min.lib.size = min.lib.size,
                                      min.reads = min.reads)
  corpus <- STdeconvolve::restrictCorpus(counts = counts, nTopOD = nTopOD)

  ldas <- STdeconvolve::fitLDA(counts = t(as.matrix(corpus)), Ks = c(length(cell_type)),
                               plot = FALSE,verbose = FALSE)

  optLDA <- STdeconvolve::optimalModel(models = ldas, opt = c(length(cell_type)))

  results <- STdeconvolve::getBetaTheta(lda = optLDA, corpus = t(as.matrix(corpus)),
                                        betaScale = betaScale,verbose = FALSE)
  deconProp <- results$theta

  #### annotation
  deconGexp <- results$beta * 1000
  mobProxyTheta2 <- model.matrix(~0 + sc_label)
  rownames(mobProxyTheta2) <- colnames(sc_exp)

  mobProxyTheta2 <- sc_exp %*% mobProxyTheta2

  corMtx_beta <- STdeconvolve::getCorrMtx(m1 = as.matrix(deconGexp), # the deconvolved cell-type `beta` (celltypes x genes)
                                          m2 = t(as.matrix(mobProxyTheta2)), # the reference `beta` (celltypes x genes)
                                          type = "b", verbose = FALSE)
  rownames(corMtx_beta) <- paste0("decon_", seq(nrow(corMtx_beta)))
  # correlationPlot(mat = corMtx_beta)
  pairs_used <- STdeconvolve::lsatPairs(t(corMtx_beta))
  m <- t(corMtx_beta)[pairs_used$rowix, pairs_used$colsix]
  # correlationPlot(t(m))
  ### extract label
  ct_annotaion_idx <- sub("sc_label","",rownames(m))
  topic_annotation_idx <- sub("decon_", "", colnames(m))

  deconProp_final <- deconProp[,topic_annotation_idx]
  colnames(deconProp_final) <- ct_annotaion_idx
  deconProp_final <- deconProp_final[,cell_type]
  deconProp_final <- sweep(deconProp_final, 1, rowSums(deconProp_final), "/")
  return(deconProp_final)
}

####### apply DestVI #########
DestVI_run <- function(database, python_env = NULL,
                       n_top_genes = 2000, use_gpu = FALSE,
                       max_iter_sc = 400, max_iter_st = 3000){

  if(is.null(python_env)){
    cat("if you want to use the environment of your own ptython, please set the path of your python environment
        in yout computer. Or the code will install Miniconda and scvi-tools, scanpy, anndata and numpy to run DestVI.  We recommond
        the user to construct anaconda env with the requirments.yaml by run the code <conda env create -f requirments.yaml> we provided.")
  }else{
    reticulate::use_python(python_env, require = T)
    reticulate::py_config()
  }
  ### extract SC, CL, ST from database
  sc_exp <- database$sc_exp
  sc_label <- database$sc_label
  spot_exp <- database$spot_exp
  cell_type <- sort(unique(sc_label))

  if(nrow(sc_exp) != nrow(spot_exp))
    stop("The features of scRNA-seq and stRNA-seq need be equal!")

  scvi <- reticulate::import("scvi")
  np <- reticulate::import("numpy")
  anndata <- reticulate::import("anndata")
  pd <- reticulate::import("pandas")
  sc <- reticulate::import("scanpy")

  ############## training for scRNA-seq
  ## process of scRNA-seq data
  meta.data <- data.frame(cell_type = sc_label)
  rownames(meta.data) <- colnames(sc_exp)
  Seurat.sc <- Seurat::CreateSeuratObject(counts = sc_exp, meta.data = meta.data)
  Seurat.sc <- Seurat::NormalizeData(Seurat.sc, verbose = FALSE)
  Seurat.sc <- Seurat::FindVariableFeatures(Seurat.sc, nfeatures = n_top_genes, verbose = FALSE)

  sc_exp_norm <- as.matrix(Seurat.sc@assays$RNA@data)[Seurat.sc@assays$RNA@var.features,]
  sc_label_filter <- Seurat.sc$cell_type

  ## process of stRNA-seq data
  Seurat.st <- Seurat::CreateSeuratObject(counts = spot_exp)
  Seurat.st <- Seurat::NormalizeData(Seurat.st, verbose = FALSE)
  spot_exp_norm <- as.matrix(Seurat.st@assays$RNA@data)[rownames(sc_exp_norm), ]

  ## training scRNA-seq data
  adata_sc <- input_for_py(expr = sc_exp_norm, cell_type = sc_label_filter)
  scvi$model$CondSCVI$setup_anndata(adata_sc, labels_key = "cell_type")
  sc_model <- scvi$model$CondSCVI(adata_sc, weight_obs = FALSE)
  scvi$model$CondSCVI$train(sc_model, use_gpu = use_gpu, max_epochs = as.integer(max_iter_sc))

  ## training stRNA-seq data
  adata_st <- input_for_py(expr = spot_exp_norm)
  scvi$model$DestVI$setup_anndata(adata_st)
  st_model = scvi$model$DestVI$from_rna_model(adata_st, sc_model)
  scvi$model$DestVI$train(st_model, use_gpu = use_gpu, max_epochs = as.integer(max_iter_st))

  ## extract proportaion matrix
  weight_matrix <- as.matrix(scvi$model$DestVI$get_proportions(st_model))[,cell_type]
  weight_matrix_norm <- sweep(weight_matrix, 1, rowSums(weight_matrix), "/")
  return(weight_matrix_norm)
}




####### utils of interface with python using R scripts##########
input_for_py <- function(expr, cell_type = NULL){
  # pip install scanpy
  sc <- reticulate::import("scanpy")
  anndata <- reticulate::import("anndata")

  utils::write.csv(expr, "count.csv")
  adata <- sc$read_csv("count.csv",  first_column_names = TRUE, dtype='float32')$T
  if(! is.null(cell_type))
    adata$obs["cell_type"] <- cell_type

  unlink("count.csv", recursive=TRUE)
  return(adata)
}


data_frame_extracted <- function(adata){
  adata <- reticulate::py_to_r(adata)
  data_frame_total <- adata$varm['means_per_cluster_mu_fg']
  # slice_idx_tmp <- adata$uns['mod']['factor_names']$factor_names
  slice_idx_tmp <- adata$uns$mod$factor_names
  slice_idx <- paste0("means_per_cluster_mu_fg_", slice_idx_tmp)
  data_frame_extracted <- data_frame_total[slice_idx]
  names(data_frame_extracted) <- slice_idx_tmp
  results <- reticulate::r_to_py(data_frame_extracted)
  return(results)
}

################ run CARD ################
CARD_run <- function(database, minCountGene = 100, minCountSpot = 5){
  ### extract SC, CL, ST, ST.Loc from database
  sc_exp <- database$sc_exp
  sc_label <- database$sc_label
  spot_exp <- database$spot_exp
  spot_loc <- database$spot_loc
  cell_type <- sort(unique(sc_label))
  spot_names_ori <- colnames(spot_exp)
  ###
  spot_names <- paste0("spot", 1:nrow(spot_loc))
  colnames(spot_exp) <- spot_names
  rownames(spot_loc) <- spot_names
  colnames(spot_loc) <- c("x", "y")

  sampleInfo <- rep("sample1", ncol(sc_exp))
  names(sampleInfo) <- colnames(sc_exp)
  meta.data <- data.frame(cellID = colnames(sc_exp), cellType = sc_label, sampleInfo = sampleInfo)
  rownames(meta.data) <- colnames(sc_exp)

  CARD_obj <- CARD::createCARDObject(sc_count = sc_exp,
                                     sc_meta = meta.data,
                                     spatial_count = spot_exp,
                                     spatial_location = spot_loc,
                                     ct.varname = "cellType",
                                     ct.select = cell_type,
                                     sample.varname = "sampleInfo",
                                     minCountGene = minCountGene,
                                     minCountSpot = minCountSpot)
  CARD_obj <- CARD::CARD_deconvolution(CARD_object = CARD_obj)

  CARD_results <- as.matrix(CARD_obj@Proportion_CARD)
  # ccc <- setdiff(colnames(spot_exp),colnames(CARD_obj@spatial_countMat))
  #
  # cc <- spot_exp[,ccc]
  CARD_results <- CARD_results[colnames(spot_exp), cell_type]
  rownames(CARD_results) <- spot_names_ori
  return(CARD_results)
}


####### utilize of interface with python using R scripts ##########

#' Running each base deconvolution method individually to obtain the  base cell type deconvolution results on spatially resolved transcriptomics data.
#'
#' This function is implemented to perform individual deconvolution methods.  The current implementation of
#' EnDecon integrates twelve state-of-the-art methods:  CARD, cell2location, DeconRNASeq, DWLS, MuSiC (MuSiC weighted and MuSiC all gene), RCTD, SCDC, SpatialDWLS,
#' SPOTlight,Stereoscope, and SVR. These packages will be automatically installed along
#' with EnDecon. \cr
#'
#' @import pcaMethods
#' @import data.table
#' @import foreach
#' @importFrom Matrix colSums rowSums rowMeans
#' @importFrom methods new as
#' @importFrom Biobase ExpressionSet
#' @importFrom SCDC SCDC_prop
#' @importFrom spacexr Reference SpatialRNA create.RCTD run.RCTD
#' @importFrom MuSiC music_prop
#' @importFrom DeconRNASeq DeconRNASeq
#' @importFrom DWLS solveSVR
#' @importFrom DWLS solveDampenedWLS
#' @importFrom Seurat CreateSeuratObject SCTransform FindAllMarkers
#' @importFrom SPOTlight downsample_se_obj train_nmf topic_profile_per_cluster_nmf mixture_deconvolution_nmf
#' @importFrom Giotto createGiottoInstructions createGiottoObject normalizeGiotto findMarkers_one_vs_all calculateHVG fDataDT runPCA createNearestNetwork doLeidenCluster runDWLSDeconv
#' @importFrom reticulate use_python py_config import py_to_r r_to_py
#' @importFrom graphics hist par
#' @importFrom stats cor median
#' @importFrom utils head
#' @importFrom NMF basis
#' @importFrom NMF coef
#' @importFrom CARD createCARDObject CARD_deconvolution
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom STdeconvolve cleanCounts restrictCorpus fitLDA optimalModel getBetaTheta getCorrMtx lsatPairs
#'
#'
#' @param sc_exp scRNA-seq matrix, genes * cells. The format should be raw-counts. The matrix need include gene names and cell names.
#' @param sc_label cell type information. The cell need be divided into multiple categories.
#' @param spot_exp stRNA-seq matrix, genes * spots. The format should be raw-counts. The matrix need include gene names and spot names.
#' @param spot_loc coordinate matrix, spots * coordinates. The matrix need include spot names and coordinate name (x, y).
#' @param gene_det_in_min_cells_per a floor variable. minimum percent # of genes that need to be detected in a cell.
#' @param expression_threshold a floor variable. Threshold to consider a gene expressed.
#' @param nUMI a floor variable. 	minimum # of read count that need to be detected in a cell or spot.
#' @param verbose a logical variable that defines whether to print the processing flow of data process.
#' @param plot a logical variable that defines whether to plot the selected genes and selected cell expression.
#' @param python_env the path of python environment. We recommend user construct python environment by the .yml provided by ours.
#' @param use_gpu a logical variable whether to use GPU to train Stereoscope and cell2location.
#' @param saving_results a logical variable whether to save the results of individual deconvolution methods.
#' @param SCDC a logical variable whether to apply SCDC.
#' @param RCTD a logical variable whether to apply RCTD.
#' @param MuSiC a logical variable whether to apply MuSiC all gene and MuSiC weighted.
#' @param DeconRNASeq a logical variable whether to apply DeconRNASeq.
#' @param DestVI a logical variable whether to apply DestVI.
#' @param DWLS a logical variable whether to apply DWLS and SVR.
#' @param SPOTlight a logical variable whether to apply SPOTlight.
#' @param SpatialDWLS a logical variable whether to apply SpatialDWLS.
#' @param Stereoscope a logical variable whether to apply Stereoscope.
#' @param cell2location a logical variable whether to apply cell2location.
#' @param CARD a logical variable whether to apply CARD.
#' @param STdeconvolve a logical variable whether to apply STdeconvolve.
#' @param SCDC.iter.max a integer variable represents the maximum number of iteration in WNNLS of SCDC.
#' @param RCTD.CELL_MIN_INSTANCE a integer value represent the min cells in one cell type for reference scRAN-seq.
#' @param MuSiC.iter.max a integer variable represents maximum iteration number of MuSiC training.
#' @param MuSiC.nu a floor variable represents regulation parameter in MuSiC model.
#' @param MuSiC.eps a floor variable represents threshold of convergence of training model.
#' @param DeconRNASeq.perc a floor variable represents the values for filter cells.
#' @param DWLS.parallel a logical variable indicating whether to apply DWL with multiple CPU. Default setting is TRUE.
#' @param DWLS.is_select_DEGs a logical variable indicating whether to select genes for each cell type of scRNA-seq dataset. Default setting is TRUE.
#' @param SPOTlight.cl_n integer variable  indicating how many cells to keep from each cluster. If a cluster has n < cl_n then all cells will be selected, if it has more then cl_n will be sampled randomly. Default value is 100.
#' @param SPOTlight.hvg integer variable that represents number of highly variable genes to use on top of the marker genes. Default values is 3000.
#' @param SPOTlight.min_cont floor variable indicates the minimum contribution we expect from a cell in that spot. Default values is 0.001.
#' @param SpatialDWLS.findmarker_method a string vector indicating method to use to detect differentially expressed genes.
#' @param SpatialDWLS.ncp_spa a integer value indicating number of principal components to calculate. Default setting is 100.
#' @param SpatialDWLS.dimensions_to_use	 a integer value indicating number of dimensions to use as input for constructing KNN network. Default setting is 10.
#' @param SpatialDWLS.k a integer value indicating number of k neighbors to use for constructing KNN network. Default setting is 10.
#' @param SpatialDWLS.resolution  resolution in doLeidenCluster function in Giotto package. Default setting is 0.4.
#' @param SpatialDWLS.n_iterations number of interations to run the Leiden algorithm. If the number of iterations is negative, the Leiden algorithm is run until an iteration in which there was no improvement.
#' @param SpatialDWLS.n_cell number of cells per spot. Default setting is 50.
#' @param SpatialDWLS.is_select_DEGs a logical value whether to select genes before applying for the SpatialDWLS.
#' @param Stereoscope.sc_training_plot a logical variable whether to plot the training loss indicating whether to increase the number of maximum epoch for training for scRNA-seq dataset. Default setting is FALSE.
#' @param Stereoscope.sc_training_save_trained_model a logical variable whether to save the trained model for scRNA-seq dataset. Default setting is FALSE.
#' @param Stereoscope.sc_max_epochs an integer variable indicating the maximum epoches for training scRNA-seq. Default setting is 400.
#' @param Stereoscope.sc_lr  an integer variable indicating the learning rate for training scRNA-seq. Default setting is 0.01.
#' @param Stereoscope.select_HVG a logical variable whether to select highly variable genes for the scRNA-seq data. Default setting is TRUE.
#' @param Stereoscope.HVG_num number of selected highly variable genes if Stereoscope.select_HVG = TRUE. Default setting is 5000.
#' @param Stereoscope.st_training_plot a logical variable whether to plot the training loss indicating whether to increase the number of maximum epoch for training for stRNA-seq dataset. Default setting is FALSE.
#' @param Stereoscope.st_training_save_trained_model a logical variable whether to plot the training loss indicating whether to increase the number of maximum epoch for training for stRNA-seq dataset. Default setting is FALSE.
#' @param Stereoscope.st_max_epochs an integer variable indicating the maximum epoches for training sTRNA-seq. Default setting is 400.
#' @param Stereoscope.st_lr  an integer variable indicating the learning rate for training sTRNA-seq. Default setting is 0.01.
#' @param cell2location.sc_max_epoches  an integer variable indicating the maximum epoches for training scRNA-seq.
#' @param cell2location.sc_lr  an integer variable indicating the learning rate for training scRNA-seq.
#' @param cell2location.st_N_cells_per_location a integer variable indicating the number of cells in each spot.
#' @param cell2location.st_detection_alpha a floor variable indicating the super-parameter of regularization.
#' @param cell2location.st_max_epoches  an integer variable indicating the maximum epoches for training stRNA-seq.
#' @param CARD.minCountGene an integer variable indicating the minimum counts for each gene for the construct CARD object. Default setting is 100.
#' @param CARD.minCountSpot an integer variable indicating the minimum counts for each spatial location. Default setting is 5.
#' @param STdeconvolve.min.lib.size Minimum number of genes detected in a cell. Cells with fewer genes will be removed for the stRNA-seq dataset. Default setting is 1.
#' @param STdeconvolve.min.reads Minimum number of reads per gene. Genes with fewer reads will be removed. Default setting is 1.
#' @param STdeconvolve.nTopOD Number of top over-dispersed genes to use for the stRNA-seq data. Default setting is 1000.
#' @param STdeconvolve.betaScale Factor to scale the predicted cell-type gene expression profiles. Default setting is 1000.
#' @param DestVI.n_top_genes Number of selected HVGs by Seurat.V3. Default setting is 2000.
#' @param DestVI.max_iter_sc Maximum number of epoches for the training of scRNA-seq data. Default setting is 400.
#' @param DestVI.max_iter_st Maximum number of epoches for the training of stRNA-seq data. Default setting is 3000.
#'
#' @return a list contains all the results inferred by individual deconvolution methods and the times of running individual methods. The elements of list is a matrix, spots * cell-type and a time vector.
#'
#' @export
EnDecon_individual_methods <- function(sc_exp, sc_label, spot_exp, spot_loc,
                                       gene_det_in_min_cells_per = 0.01, expression_threshold = 1 ,
                                       nUMI = 100, verbose = FALSE, plot= FALSE, python_env = NULL, use_gpu = FALSE, saving_results = FALSE,
                                       SCDC = TRUE, RCTD = TRUE, MuSiC = TRUE, DeconRNASeq = TRUE, DestVI = TRUE, DWLS = TRUE, SPOTlight = TRUE,
                                       SpatialDWLS = TRUE, Stereoscope = TRUE, cell2location = TRUE, CARD = TRUE, STdeconvolve = TRUE,
                                       SCDC.iter.max = 1000, RCTD.CELL_MIN_INSTANCE = 10,MuSiC.iter.max = 1000, MuSiC.nu = 1e-04,
                                       MuSiC.eps = 0.01, DeconRNASeq.perc = 0.05,
                                       DWLS.parallel = TRUE, DWLS.is_select_DEGs = TRUE,
                                       SPOTlight.cl_n = 100, SPOTlight.hvg = 3000, SPOTlight.min_cont = 0.001,
                                       SpatialDWLS.findmarker_method = "gini", SpatialDWLS.ncp_spa = 100, SpatialDWLS.dimensions_to_use = 10, SpatialDWLS.k = 10,
                                       SpatialDWLS.resolution = 0.4, SpatialDWLS.n_iterations = 1000, SpatialDWLS.n_cell = 50, SpatialDWLS.is_select_DEGs = TRUE,
                                       Stereoscope.sc_training_plot = FALSE, Stereoscope.sc_training_save_trained_model = FALSE,
                                       Stereoscope.sc_max_epochs = 10000, Stereoscope.sc_lr = 0.01,
                                       Stereoscope.select_HVG = TRUE, Stereoscope.HVG_num = 5000,
                                       Stereoscope.st_training_plot = FALSE, Stereoscope.st_training_save_trained_model = FALSE,
                                       Stereoscope.st_max_epochs = 10000, Stereoscope.st_lr = 0.01,
                                       cell2location.sc_max_epoches = 1000, cell2location.sc_lr = 0.002, cell2location.st_N_cells_per_location = 30,
                                       cell2location.st_detection_alpha = 200.00, cell2location.st_max_epoches = 10000,
                                       CARD.minCountGene = 100, CARD.minCountSpot = 5,
                                       STdeconvolve.min.lib.size = 100, STdeconvolve.min.reads = 1, STdeconvolve.nTopOD = 1000,
                                       STdeconvolve.betaScale = 1000,
                                       DestVI.n_top_genes = 2000, DestVI.max_iter_sc = 400, DestVI.max_iter_st = 3000){
  if(is.null(python_env)){
    cat("if you want to use the environment of your own python, please set it manual.
         Otherwise, the code will create R_reticulate environment and install  scvi-tools and cell2location.
        We strongly recommend the user to construct anaconda env and install python packages by running the .yml file we provided.")
  }else{
    reticulate::use_python(python_env, require = T)
    reticulate::py_config()
  }

  Methods <- c("CARD", "cell2location", "DeconRNASeq", "DestVI", "DWLS", "MuSiC", "RCTD", "SCDC", "SpatialDWLS","SPOTlight", "Stereoscope", "STdeconvolve")

  Methods.idx <- c(CARD, cell2location, DeconRNASeq, DestVI, DWLS, MuSiC, RCTD, SCDC, SpatialDWLS, SPOTlight, Stereoscope, STdeconvolve)

  Methods.used <- Methods[Methods.idx]

  K <- length(Methods.used)

  database <- data_process(sc_exp, sc_label, spot_exp, spot_loc, gene_det_in_min_cells_per, expression_threshold, nUMI, verbose, plot)

  Results.Deconv <- list()

  k <- 1
  times_methods <- c()
  # apply CARD
  if(CARD){
    cat("Execute CARD analysis....\n")

    pred.CARD <- Sys.time()
    Results.Deconv$CARD <- CARD_run(database, minCountGene = CARD.minCountGene,
                                    minCountSpot = CARD.minCountSpot)
    end.CARD <- Sys.time()
    time.CARD <- difftime(end.CARD, pred.CARD, units = "mins")
    cat("Run time for CARD: ", time.CARD, "min","\n")
    times_methods <- c(times_methods, time.CARD)
    k <- k + 1
  }

  # apply cell2location
  if(cell2location){
    cat("Execute cell2location analysis....\n")

    pred.cell2location <- Sys.time()
    Results.Deconv$cell2location <- cell2location_run(database, python_env = python_env,
                                                      sc_max_epoches = cell2location.sc_max_epoches, sc_lr = cell2location.sc_lr,
                                                      sc_use_gpu = use_gpu, st_N_cells_per_location = cell2location.st_N_cells_per_location,
                                                      st_detection_alpha = cell2location.st_detection_alpha, st_max_epoches = cell2location.st_max_epoches)
    end.cell2location <- Sys.time()
    time.cell2location <- difftime(end.cell2location, pred.cell2location, units = "mins")
    cat("Run time for cell2location: ", time.cell2location, "min","\n")
    times_methods <- c(times_methods, time.cell2location)
    k <- k + 1
  }

  # apply DeconRNASeq
  if(DeconRNASeq){
    cat("Execute DeconRNASeq analysis....\n")

    pred.DeconRNASeq <- Sys.time()
    Results.Deconv$DeconRNASeq <- DeconRNASeq_run(database, perc = DeconRNASeq.perc)
    end.DeconRNASeq <- Sys.time()
    time.DeconRNASeq <- difftime(end.DeconRNASeq, pred.DeconRNASeq, units = "mins")
    cat("Run time for DeconRNASeq: ", time.DeconRNASeq, "min","\n")
    times_methods <- c(times_methods, time.DeconRNASeq)
    k <- k + 1
  }

  # apply DestVI
  if(DestVI){
    cat("Execute DestVI analysis....\n")

    pred.DestVI <- Sys.time()
    Results.Deconv$DestVI <- DestVI_run(database, python_env = python_env,
                                        n_top_genes = DestVI.n_top_genes, use_gpu = use_gpu,
                                        max_iter_sc = DestVI.max_iter_sc, max_iter_st = DestVI.max_iter_st)
    end.DestVI <- Sys.time()
    time.DestVI <- difftime(end.DestVI, pred.DestVI, units = "mins")
    cat("Run time for DestVI: ", time.DestVI, "min","\n")
    times_methods <- c(times_methods, time.DestVI)
    k <- k + 1
  }

  # apply DWLS
  if(DWLS){
    cat("Execute DWLS and SVR analysis, the DWLS is slow and non-efficient....\n")

    pred.DWLS <- Sys.time()
    temp <- DWLS_run(database, parallel = DWLS.parallel, is_select_DEGs = DWLS.is_select_DEGs, python_env = python_env)
    Results.Deconv$DWLS <- temp$DampenedWLS
    Results.Deconv$SVR <- temp$SVR
    end.DWLS <- Sys.time()
    time.DWLS <- difftime(end.DWLS, pred.DWLS, units = "mins")
    cat("Run time for DWLS: ", time.DWLS, "min","\n")
    times_methods <- c(times_methods, temp$time.DWLS, temp$time.SVR)
    k <- k + 2
  }

  # apply MuSiC
  if(MuSiC){
    cat("Execute MuSiC_weight and MuSiC_allgene analysis....\n")

    pred.MuSiC <- Sys.time()
    temp <- MuSiC_run(database, iter.max = MuSiC.iter.max, nu = MuSiC.nu, eps = MuSiC.eps)
    Results.Deconv$MuSiC_weight <- temp$Music_weight
    Results.Deconv$MuSiC_allgene <- temp$Music_allgene
    end.MuSiC <- Sys.time()
    time.MuSiC <- difftime(end.MuSiC, pred.MuSiC, units = "mins")
    cat("Run time for MuSiC: ", time.MuSiC, "min","\n")
    times_methods <- c(times_methods, time.MuSiC/2, time.MuSiC/2)
    k <- k + 2
  }

  # apply RCTD
  if(RCTD){
    cat("Execute RCTD analysis....\n")

    pred.RCTD <- Sys.time()
    Results.Deconv$RCTD <- RCTD_run(database, CELL_MIN_INSTANCE = RCTD.CELL_MIN_INSTANCE)
    end.RCTD <- Sys.time()
    time.RCTD <- difftime(end.RCTD, pred.RCTD, units = "mins")
    cat("Run time for RCTD: ", time.RCTD, "min","\n")
    times_methods <- c(times_methods, time.RCTD)
    k <- k + 1
  }

  # apply SCDC
  if(SCDC){
    cat("Execute SCDC analysis....\n")

    pred.SCDC <- Sys.time()
    Results.Deconv$SCDC <- suppressMessages(SCDC_run(database, iter.max = SCDC.iter.max))
    end.SCDC <- Sys.time()
    time.SCDC <- difftime(end.SCDC, pred.SCDC, units = "mins")
    cat("Run time for SCDC: ", time.SCDC, "min","\n")
    times_methods <- c(times_methods, time.SCDC)
    k <- k + 1
  }

  # apply SpatialDWLS
  if(SpatialDWLS){
    cat("Execute SpatialDWLS analysis....\n")

    pred.SpatialDWLS <- Sys.time()
    Results.Deconv$SpatialDWLS <- suppressWarnings(spatialDWLS_run(database, my_python_path = python_env, findmarker_method = SpatialDWLS.findmarker_method, is_select_DEGs = SpatialDWLS.is_select_DEGs,
                                                                   ncp_spa = SpatialDWLS.ncp_spa, dimensions_to_use = SpatialDWLS.dimensions_to_use, k = SpatialDWLS.k,
                                                                   resolution = SpatialDWLS.resolution, n_iterations = SpatialDWLS.n_iterations, n_cell = SpatialDWLS.n_cell))
    end.SpatialDWLS <- Sys.time()
    time.SpatialDWLS <- difftime(end.SpatialDWLS, pred.SpatialDWLS, units = "mins")
    cat("Run time for SpatialDWLS: ", time.SpatialDWLS, "min","\n")
    times_methods <- c(times_methods, time.SpatialDWLS)
    k <- k + 1
  }

  # apply SPOTlight
  if(SPOTlight){
    cat("Execute SPOTlight analysis....\n")

    pred.SPOTlight <- Sys.time()
    Results.Deconv$SPOTlight <- suppressWarnings(SPOTlight_run(database, cl_n = SPOTlight.cl_n, hvg = SPOTlight.hvg, min_cont = SPOTlight.min_cont))
    end.SPOTlight <- Sys.time()
    time.SPOTlight <- difftime(end.SPOTlight, pred.SPOTlight, units = "mins")
    cat("Run time for SPOTlight: ", time.SPOTlight, "min","\n")
    times_methods <- c(times_methods, time.SPOTlight)
    k <- k + 1
  }

  # apply Stereoscope
  if(Stereoscope){
    cat("Execute Stereoscope analysis....\n")

    pred.Stereoscope <- Sys.time()
    Results.Deconv$Stereoscope <- Stereoscope_run(database, python_env = python_env, use_gpu = use_gpu,
                                                  sc_training_plot = Stereoscope.sc_training_plot, sc_training_save_trained_model = Stereoscope.sc_training_save_trained_model,
                                                  sc_max_epochs = Stereoscope.sc_max_epochs, sc_lr = Stereoscope.sc_lr,
                                                  st_training_plot = Stereoscope.st_training_plot, st_training_save_trained_model = Stereoscope.st_training_save_trained_model,
                                                  st_max_epochs = Stereoscope.st_max_epochs, st_lr = Stereoscope.st_lr, select_HVG = Stereoscope.select_HVG,
                                                  HVG_num = Stereoscope.HVG_num)
    end.Stereoscope <- Sys.time()
    time.Stereoscope <- difftime(end.Stereoscope, pred.Stereoscope, units = "mins")
    cat("Run time for Stereoscope: ", time.Stereoscope, "min","\n")
    times_methods <- c(times_methods, time.Stereoscope)
    k <- k + 1
  }

  # apply STdeconvolve
  if(STdeconvolve){
    cat("Execute STdeconvolve analysis....\n")

    pred.STdeconvolve <- Sys.time()
    Results.Deconv$STdeconvolve <- STdeconvolve_run(database, min.lib.size = STdeconvolve.min.lib.size,
                                                    min.reads = STdeconvolve.min.reads, nTopOD = STdeconvolve.nTopOD,
                                                    betaScale = STdeconvolve.betaScale)
    end.STdeconvolve <- Sys.time()
    time.STdeconvolve <- difftime(end.STdeconvolve, pred.STdeconvolve, units = "mins")
    cat("Run time for STdeconvolve: ", time.STdeconvolve, "min","\n")
    times_methods <- c(times_methods, time.STdeconvolve)
    k <- k + 1
  }
  names(times_methods) <- names(Results.Deconv)
  if(saving_results)
  {
    cat("Saving individual Deconvolution analysis results...\n")

    save(Results.Deconv, file = "Results_Deconv.RData")
  }
  res <- list(Results.Deconv, times_methods)
  return(res)
}

#' Ensemble the results of individual deconvolution results.
#' This function uses the weighted median methods to integrate the results obtained by individual deconvolution methods.
#'
#' @importFrom spatstat.geom weighted.median
#' @importFrom abind abind
#' @param Results.Deconv a list contains all the results of individual deconvolution methods. The elements of list is a matrix, spots * cell-type.
#' @param lambda hyper-parameter constrain the weight of individual methods for ensemble. If the parameter is set to NULL, then, we will adopt the value in our algorithm.
#' @param prob.quantile numeric of probabilities with values in [0,1]. Default setting is 0.5.
#' @param niter a positive integer represents the maximum number of updating algorithm. Default setting is 100.
#' @param epsilon a parameter represents the stop criterion.
#'
#' @return a list contains a matrix of the ensemble deconvolution result and a vector of the weight assigned to individual methods.
#' @examples
#' data("breast.sc.ref")
#' data("breast.sc.cell.label")
#' data("breast.st")
#' data("breast.st.loc")
#' ##### path on ubuntu platform on our computer
#' python_env <- "~/.conda/envs/EnDecon_GPU/bin/python"
#' Results.dec.mouse <- EnDecon_individual_methods(sc_exp = breast.sc.ref,
#' sc_label = breast.sc.cell.label, spot_exp = breast.st,
#' spot_loc = breast.st.loc, python_env = python_env,
#' use_gpu = TRUE,gene_det_in_min_cells_per = 0.01,
#' RCTD.CELL_MIN_INSTANCE = 5, saving_results = FALSE)
#' ensemble.results <- solve_ensemble(Results.dec.mouse[[1]])
#'
#'@export
solve_ensemble <- function(Results.Deconv, lambda = NULL, prob.quantile = 0.5,
                           niter = 100, epsilon = 1e-5){
  # Results.Deconv <- Results.Deconv.all[[1]]
  num.methods <- length(Results.Deconv)
  num.spots <- nrow(Results.Deconv[[1]])
  num.cell.type <- ncol(Results.Deconv[[1]])

  ## initialization V by the mean of individual values
  w <- c(rep(1/num.methods, num.methods))
  H <-  Reduce("+", Map("*", Results.Deconv, w))

  if(is.null(lambda)){
    cat("We will adpote a value for lambda in our algorithm...", "\n")
  }

  k <- 1
  while (k <= niter) {
    if(k == 1){
      loss_all_temp <- 0
      temp2 <-  sapply(Results.Deconv, L1_norm, Y = H)
      lambda <- quantile(temp2, probs = prob.quantile)
    }else{
      loss_all_temp <- loss_all
    }
    ##### update w
    temp2 <-  sapply(Results.Deconv, L1_norm, Y = H)
    w <- exp(-temp2/lambda)/sum(exp(-temp2/lambda))
    ##### update V
    temp4 <- do.call(abind::abind, c(Results.Deconv, list(along = 3)))
    H <- apply(temp4, 1:2, median_weighted, w = w)

    loss_main <- sum(sapply(Results.Deconv, L1_norm, Y = H) * w)
    loss_entropy <- sum(w * log(w))
    loss_all <- loss_main + lambda * loss_entropy

    cat("iter: ", k, "loss_main: ", loss_main, "loss_entropy: ", loss_entropy,
        "loss_all: ", loss_all, "lambda: ", lambda, "\n")
    if(k == niter)
      cat("The method maybe not convergens, the algorithm need an larger max_epoches!", "\n")

    if(abs(loss_all - loss_all_temp) < epsilon | k >= niter)
      break
    k <- k + 1
  }

  colnames(H) <- colnames(Results.Deconv[[1]])
  H_norm <- sweep(H, 1, rowSums(H), "/")
  return(list(H_norm = H_norm, w = w))
}

L1_norm <- function(X, Y){
  return(sum(abs(X-Y)))
}

median_weighted <- function(x, w){
  index_non_zero <- which(x != 0, arr.ind = TRUE)
  x_use <- x[index_non_zero]
  w_use <- w[index_non_zero]
  results <- spatstat.geom::weighted.median(x,w)
  return(results)
}
