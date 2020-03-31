#' Trait selection
#'
#' This function determines whether the selected traits exhibit or not a clustering/overdispersion
#' signal on the tested samples. For each trait, compares the observed Mean Pairwise Distance (MPD)
#' of each sample against a distribution of synthetic commmunities MPDs obtained by a randomization
#' test. Each synthetic community is build maintaining the original sample richness and randomly
#' selecting organisms form the global pool.
#'
#' @param table1 A data frame containing organisms names on the first column and its trait values on the
#' consecutive ones. It also has to contain two columns with the maximum and the minimum values of the tested
#' environmental variable where the organisms have been observed.
#'
#' @param table2 A presence-absence observations table with the organisms names on the first column and the
#' sample names as consecutive colnames.
#'
#' @param table3 A dataframe containing sample names on the first column and environmental parameters on the
#' consecutive ones.
#'
#' @param traits_columns Table 1 column numbers where different trait values appear.
#'
#' @param min_env_col Table 1 column number indicating the minimum value of the environmental variable were each
#' organism has been observed.
#'
#' @param max_env_col Table 1 column number indicating the maximum value of the environmental variable were each
#' organism has been observed.
#'
#' @param env_var_col Table 2 column number indicating the tested environmental variable.
#'
#' @param repetitions Number of simulated synthetic communities distributions.
#'
#' @return The function returns a dataframe with trait names as colnames and the p-value distribution of the different
#' traits.
#'
#' @examples
#'
#' \donttest{
#' data(group_information)
#' data(table_presence_absence)
#' data(metadata)
#' rtcc1(group_information, table_presence_absence, metadata, 2:11, 12, 12, 2, 100)
#' }
#'
#'
#' @export


rtcc1 <- function(table1, table2, table3, traits_columns, min_env_col, max_env_col, env_var_col, repetitions){

  dataset <- table1
  dataset <- as.data.frame(dataset)
  table_presence_absence <- table2
  metadata <- table3

  colnames(dataset)[min_env_col] <- "min"
  colnames(dataset)[max_env_col] <- "max"
  colnames(dataset)[1] <- "tax"

  colnames(metadata)[1] <- "sample_ID"
  colnames(metadata)[env_var_col] <- "env_variable"

  colnames(table_presence_absence)[1] <- "tax"

  pool <- as.vector(dataset$tax)
  l_pool <- length(pool)
  richnes <- vegan::specnumber(t(table_presence_absence[,-1]))

  metadata$sample_ID <- factor(metadata$sample_ID, levels = unique(metadata$sample_ID[order(metadata$env_variable)]))
  local_communities <- metadata$sample_ID

  y <- 1
  h_iteration <- 1
  results <- as.data.frame(matrix(nrow = length(local_communities), ncol = length(traits_columns), NA))
  colnames(results) <- colnames(dataset)[traits_columns]

  for(i in traits_columns){

    # Calculate MPDs for observed communities
    real_MPDs <- real_mpd(local_communities, table_presence_absence, dataset, i)

    trait <- as.numeric(dataset[,i])
    p <- 1
    species_presence <- rowSums(table_presence_absence[,colnames(table_presence_absence) %in% local_communities])
    p_values <- rep(NA, length(local_communities))
    species_presence[species_presence != 0] <- 1
    prob_vector <- species_presence
    for(m in local_communities){
      MPDs <- rep(NA, (repetitions))
      rich_m <- richnes[m]
      ncomps <- rich_m * (rich_m-1)
      MPDs <- sample_and_mpd(repetitions, h_iteration, trait, l_pool, rich_m, prob_vector, MPDs, ncomps)
      p_values[p] <- sum(MPDs < real_MPDs[p,2])/length(MPDs)
      p <- p + 1
    }
    results[,y] <- p_values
    y <- y + 1
  }
  return(results)
}
