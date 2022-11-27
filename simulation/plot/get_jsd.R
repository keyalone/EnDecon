## Compute JSD
# tha package


get_jsd <- function(mtrx1, mtrx2) {
  ## spot * cell types
  ##### Get TRUE JSD between real-predicted proportions #####
  # Initialize matrix
  true_jsd_mtrx <- matrix(nrow = nrow(mtrx1), ncol = 1)
  
  # Loop over all the rows
  for (i in seq_len(nrow(mtrx1))) {
    # Create matrix to feed to JSD
    x <- rbind(mtrx1[i, ],
               mtrx2[i, ])
    # Calculate JSD and save it in true_JSD_mtrx
    if (sum(mtrx2[i, ]) > 0) {
      true_jsd_mtrx[i] <- suppressMessages(philentropy::JSD(x = x, unit = "log2",
                                                            est.prob = "empirical"))
    } else {
      true_jsd_mtrx[i] <- 1
    }
  }
  
  return(true_jsd_mtrx)
}
