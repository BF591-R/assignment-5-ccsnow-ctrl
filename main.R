library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')
library('biomaRt')
library('testthat')
library('fgsea')

#' Function to generate a SummarizedExperiment object with counts and coldata
#' to use in DESeq2
#'
#' @param csv_path (str): path to the file verse_counts.tsv
#' @param metafile (str): path to the metadata sample_metadata.csv
#' @param selected_times (list): list of sample timepoints to use
#' 
#'   
#' @return SummarizedExperiment object with subsetted counts matrix
#'   and sample data. Ensure that the timepoints column used as input 
#'   to the model design has 'vP0' set as the reference factor level. Your 
#'   colData dataframe should have columns named samplename and timepoint.
#' @export
#' Define a function that reads from these two files to construct a
 # counts matrix and matching sample dataframe containing only samples belonging to
# the vP0 and vAd timepoints. This should correspond to a total of 4 samples
#(vP0_1, vP0_2, vAd_1, vAd_2). Store this counts matrix and the sample dataframe
#in a SummarizedExperiments object. Ensure that the name of your counts matrix
#stored as an assay in the SummarizedExperiment object is `counts`.
#' @examples se <- make_se('verse_counts.tsv', 'sample_metadata.csv', c('vP0', 'vAd'))
make_se <- function(counts_csv, metafile_csv, selected_times) {
  metafile <- read_csv(metafile_csv)
  counts <- read_tsv(counts_csv)
  
  metafile <- metafile %>% 
    filter(timepoint %in% selected_times)
  
  metafile$timepoint <- factor(metafile$timepoint, levels = selected_times)
  
  counts <- counts %>% 
    column_to_rownames("gene") %>% 
    dplyr::select(all_of(metafile$samplename)) 
  
  se <- SummarizedExperiment(
    assays = list(counts = as.matrix(counts)),
    colData = DataFrame(
      samplename = metafile$samplename,
      timepoint = metafile$timepoint,
      row.names = metafile$samplename
    )
  )
  
  return(se)
}

#' Function that runs DESeq2 and returns a named list containing the DESeq2
#' results as a dataframe and the dds object returned by DESeq2
#'
#' @param se (obj): SummarizedExperiment object containing counts matrix and
#' coldata
#' @param design: the design formula to be used in DESeq2
#'
#' @return list with DESeqDataSet object after running DESeq2 and results from
#'   DESeq2 as a dataframe
#' @export
#'
#' @examples results <- return_deseq_res(se, ~ timepoint)
return_deseq_res <- function(se, design) {
  dds <- DESeqDataSet(se, design = design)
  dds <- DESeq(dds)
  res <- as.data.frame(results(dds))
  return(list(dds = dds, res = res))
}


#' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#' adds a column to denote plotting status in volcano plot. Column should denote
#' whether gene is either 1. Significant at padj < .10 and has a positive log
#' fold change, 2. Significant at padj < .10 and has a negative log fold change,
#' 3. Not significant at padj < .10. Have the values for these labels be UP,
#' DOWN, NS, respectively. The column should be named `volc_plot_status`. Ensure
#' that the column name for your rownames is called "genes". 
#'
#' @param deseq2_res (df): results from DESeq2 
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return Tibble with all columns from DESeq2 results and one additional column
#'   labeling genes by significant and up-regulated, significant and
#'   downregulated, and not significant at padj < .10.
#'   
#' @export
#'
#' @examples labeled_results <- label_res(res, .10)
#' 
label_res <- function(deseq2_res, padj_threshold) {
  res <- deseq2_res %>%
    rownames_to_column("genes") %>%
    as_tibble()
  res <- res %>%
    mutate(volc_plot_status = case_when(
      padj < padj_threshold & log2FoldChange > 0 ~ "UP",
      padj < padj_threshold & log2FoldChange < 0 ~ "DOWN",
      TRUE ~ "NS"
    ))
  
  return(res)
}

#' Function to plot the unadjusted p-values as a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#'
#' @return ggplot: a histogram of the raw p-values from the DESeq2 results
#' @export
#'
#' @examples pval_plot <- plot_pvals(labeled_results)
plot_pvals <- function(labeled_results) { 
  plot <- ggplot(labeled_results, aes(x = pvalue)) +
    geom_histogram(binwidth = 0.01) + 
    labs( title = "Pvalues Plotted", x = "pvalue", y = "Counts") +
    theme_bw() +
    theme(plot.title= element_text(face="bold"))
  return(plot) }


#' Function to plot the log2foldchange from DESeq2 results in a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return ggplot: a histogram of log2FC values from genes significant at padj 
#' threshold of 0.1
#' @export
#'
#' @examples log2fc_plot <- plot_log2fc(labeled_results, .10)
plot_log2fc <- function(labeled_results, padj_threshold) {
  sig_results <- labeled_results %>%
  filter(padj < padj_threshold)
  plot <- ggplot(sig_results, aes(x = log2FoldChange)) + 
    geom_histogram(binwidth = 0.01) +
    labs( title = "Log2 Fold Change Plotted", x = "log2fc", y = "Counts") +   
     
    theme(plot.title= element_text(face="bold")) 
  return(plot)
}
#' Function to make scatter plot of normalized counts for top ten genes ranked
#' by ascending padj
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param dds_obj (obj): The object returned by running DESeq (dds) containing
#' the updated DESeqDataSet object with test results
#' @param num_genes (int): Number of genes to plot
#'
#' @return ggplot: a scatter plot with the normalized counts for each sample for
#' each of the top ten genes ranked by ascending padj
#' @export
#'
#' @examples norm_counts_plot <- scatter_norm_counts(labeled_results, dds, 10)
scatter_norm_counts <- function(labeled_results, dds_obj, num_genes){
  top_genes <- labeled_results %>% 
    arrange(padj) %>% 
    head(num_genes) 
  
  norm_counts <- counts(dds_obj, normalized = TRUE)
  
  top_counts <- norm_counts[top_genes$genes, ] 
  
  top_long_counts <- as.data.frame(top_counts) %>% 
    rownames_to_column("gene") %>%
    pivot_longer(
      cols = -gene, 
      names_to = "samplename", 
      values_to = "counts") 
  plot <- ggplot(top_long_counts, aes(x = samplename, y = counts, color = gene))+
           geom_jitter() +
          theme_bw() +
           labs(
             x = "Sample",
             y ="Counts",
             color = "Gene",
             title = "Normalized counts of top ten genes"
           )
  return(plot)

                 
}

#' Function to generate volcano plot from DESeq2 results
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#'
#' @return ggplot: a scatterplot (volcano plot) that displays log2foldchange vs
#'   -log10(padj) and labeled by status
#' @export
#'
#' @examples volcano_plot <- plot_volcano(labeled_results)
#' 
plot_volcano <- function(labeled_results) {
  plot <- ggplot(labeled_results, aes(x = log2FoldChange, y = -log10(padj), color = volc_plot_status)) +
    geom_point() +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_bw()+
    labs( 
      x = "Log2 Fold Change", y = "Padj", color= "Plot Status", title = "DESeq2 Results")
  return(plot)
}
#' Function to generate a named vector ranked by log2FC descending
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param id2gene_path (str): Path to the file containing the mapping of
#' ensembl IDs to MGI symbols
#'
#' @return Named vector with gene symbols as names, and log2FoldChange as values
#' ranked in descending order
#' @export
#'
#' @examples rnk_list <- make_ranked_log2fc(labeled_results, 'data/id2gene.txt')

make_ranked_log2fc <- function(labeled_results, id2gene_path) {
  id2gene <- read_tsv(id2gene_path, 
                        col_names = c("genes", "symbol"))
  
  labeled_results <- labeled_results %>%
    left_join(id2gene, by = "genes")  %>%
    distinct(symbol, .keep_all = TRUE) %>%
    filter(!is.na(symbol))
  v <- setNames(labeled_results$log2FoldChange, labeled_results$symbol)
  v <- sort(v, decreasing = TRUE)
  return(v)
  
  
}


#' Function to run fgsea with arguments for min and max gene set size
#'
#' @param gmt_file_path (str): Path to the gene sets of interest in GMT format
#' @param rnk_list (named vector): Named vector generated previously with gene 
#' symbols and log2Fold Change values in descending order
#' @param min_size (int): Minimum number of genes in gene sets to be allowed
#' @param max_size (int): Maximum number of genes in gene sets to be allowed
#'
#' @return Tibble of results from running fgsea
#' @export
#'
#' @examples fgsea_results <- run_fgsea('data/m2.cp.v2023.1.Mm.symbols.gmt', rnk_list, 15, 500)
run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
  pathways <- gmtPathways(gmt_file_path) 
  resulting<- fgsea(pathways = pathways, 
        stats = rnk_list,
        minSize = min_size,
        maxSize = max_size)
  resulting<- as_tibble(resulting)
  return(resulting)
}

#' Function to plot top ten positive NES and top ten negative NES pathways
#' in a barchart
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param num_paths (int): the number of pathways for each direction (top or
#'   down) to include in the plot. Set this at 10.
#'
#' @return ggplot with a barchart showing the top twenty pathways ranked by positive
#' and negative NES
#' @export
#'
#' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
top_pathways <- function(fgsea_results, num_paths){
  
 
  top_pos <- fgsea_results %>% 
    arrange(desc(NES)) %>% 
    head(num_paths) %>%
    mutate(direction = "Positive") 
  
  top_neg <- fgsea_results %>% 
    arrange(NES) %>% 
    head(num_paths) %>%
    mutate(direction = "Negative") 
  combined <- dplyr::bind_rows(top_pos, top_neg)
  combined$pathway <- gsub("_", " ", combined$pathway)
  plot <- ggplot(combined, aes(x = NES, y = reorder(stringr::str_wrap(pathway, width = 30), NES), fill = direction)) +
    geom_col() +
    theme_bw() +
    labs(
      x = "Normalized Enrichment Score (NES)", 
      y = "Pathway",
      fill = "Direction"
    )
  return(plot)
}

