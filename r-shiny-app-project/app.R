library(shiny)
library(tidyverse)
library(DT)
library(pheatmap)
library(DESeq2)
library('fgsea')
library(biomaRt)
library('RColorBrewer')


# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("BF591 Final Project"),
  
  tabsetPanel(
    
    # panel 1 : samples
    tabPanel(
      "Samples",
      sidebarLayout(
        sidebarPanel(
          # load sample information file in CSV format
          fileInput("sample_info", "Sample information file in CSV format:", accept = c(".csv")),
        ),
        mainPanel(
          tabsetPanel(
            # tab with summary of table that includes a summary of the type and values in each column
            tabPanel('Summary',
                     tableOutput("sample_info_summary")),
            # tab with data table displaying sample information, with sortable columns
            tabPanel('Table',
                     DT::dataTableOutput("sample_info_table")),
            # tab with histograms, density plots, or violin plots of continuous variables
            tabPanel('Plots',
                     # chose x variable to plot
                     sidebarPanel(
                       uiOutput("variable_ui"),
                       # chose binwidth
                       selectInput(
                         inputId = "binwidth", 
                         label = "Binwidth:", 
                         choices = c(0.1,1,5,10),
                         selected = 1
                       )
                     ),
                     # show histogram in main panel
                     mainPanel(
                       plotOutput("sample_info_histogram")
                       )
                     )
          )
        )
      )
    ),
    
    # panel 2 : counts
    tabPanel(
      "Counts",
      sidebarLayout(
        sidebarPanel(
          # load normalized count matrix in CSV format
          fileInput("norm_counts", "Normalized count matrix in CSV format:", accept = c(".csv")),
          # slider to include genes with at least X percentile of variance
          sliderInput("var_percentile", "Percentile of Variance:", min = 0, max = 100, value = 25),
          # slider to include genes with at least X samples that are non-zero
          sliderInput("nonzero_samples", "Minimum Non-Zero Samples:", min = 0, max = 69, value = 10)
        ),
        mainPanel(
          tabsetPanel(
            # tab with text or a table summarizing the effect of the filtering, including:
              # number of samples, total number of genes, number and % of genes passing current filter, and those not
            tabPanel('Filtering Effect',
                     tableOutput("norm_counts_summary")),
            # tab with diagnostic scatter plots, where genes passing filters are marked in a darker color, and genes filtered
              # out are lighter
            tabPanel('Diagnostic Plots',
                     sidebarPanel(
                       # choose to scale variance by log10 in diagnostic plot
                       radioButtons(
                         inputId = "scale_var", 
                         label = "Scale median count and variance by log10:", 
                         choices = c("TRUE", "FALSE"),
                         selected = "TRUE"
                       )
                     ),
                     # show histogram in main panel
                     mainPanel(
                       plotOutput("variance_plot"),
                       plotOutput("nonzero_plot")
                     )
            ),
            # tab with clustered heatmap of counts remaining after filtering
            tabPanel('Heatmap',
                     plotOutput("counts_heatmap")),
            # tab with a scatter plot of PCA projections
            tabPanel('PCA',
                     # choose PC for x axis
                     sidebarPanel(
                       selectInput(
                         inputId = "x_pc", 
                         label = "Principal Component for the X axis:", 
                         choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"),
                         selected = "PC1"
                       ),
                       # choose PC for y axis
                       selectInput(
                         inputId = "y_pc", 
                         label = "Principal Component for the Y axis", 
                         choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"),
                         selected = "PC2"
                       )
                     ),
                     # show PCA scatterplot in main panel
                     mainPanel(
                       plotOutput("counts_pca")
                     )
            ),
          )
        )
      )
    ),
    
    # panel 3: DE
    tabPanel(
      "DE",
      sidebarLayout(
        sidebarPanel(
          # load raw count matrix in CSV format
          fileInput("raw_counts", "Raw count matrix in CSV format:", accept = c(".csv")),
          # load metadata / sample info in CSV format
          fileInput("metadata", "Metadata, or sample information, file in CSV format:", accept = c(".csv"))
        ),
        mainPanel(
          tabsetPanel(
            # tab with sortable table displaying differential expression results
            tabPanel('Table',
                     DT::dataTableOutput("de_table")),
            # tab with volcano plot comparing the log2FoldChange with the adj.p-value
            tabPanel('Volcano Plot',
                     sidebarLayout(
                       sidebarPanel(
                         # significant p-adjusted value
                         sliderInput("volcano_p", "Significant p-adjusted threshold:", min = 0, max = 0.1, value = 0.05),
                         # significant log2FC values
                         sliderInput("volcano_log2FC", "Significant log2FoldChange threshold:", min = 0, max = 6, value = 0)
                       ),
                       mainPanel(plotOutput("volcano_plot"))
                     )
            )
          )
        )
      )
    ),
    
    # panel 4: GSEA
    tabPanel(
      "GSEA",
      sidebarLayout(
        sidebarPanel(
          # notify user to run DE before GSEA
          textOutput('GSEA_explanation'),
          # input gene set file
          fileInput("gene_set", "File of gene set collection:", accept = c(".gmt")),
          # choose ranking metric from DESeq2 results
          radioButtons(
            inputId = "metric", 
            label = "Ranking metric for differential expression results:", 
            choices = c("log2FoldChange","-log10(p-value)","-log10(p-adj.)"),
            selected = "log2FoldChange"
          ),
          # button to download FGSEA results
          downloadButton("download_res", "Download FGSEA results"),
          # input fgsea results
          fileInput("fgsea_results", "File of FGSEA results:"),
          # file upload button
          actionButton("upload", "Upload FGSEA results")
        ),
        mainPanel(
          tabsetPanel(
            # tab with barplot of fgsea NES for top pathways
            tabPanel('Barplot',
                     sidebarPanel(
                       # choose number of top pathways to plot by adjusted p-value
                       sliderInput("n_pathways", "Number of top pathways to plot (determined by adjusted p-value):", min = 1, max = 50, value = 10)
                     ),
                     mainPanel(
                       plotOutput("fgsea_barplot")
                     )
            ),
            # tab with sortable table displaying differential expression results
            tabPanel('Table',
                     sidebarPanel(
                       # filter table by adjusted p-value
                       sliderInput("padj_filter_table", "Filter table by adjusted p-value:", min = 0, max = 1, value = 0.05),
                       # radio buttons to select all, positive, or negative NES pathways
                       radioButtons(
                         inputId = "path_filter", 
                         label = "Filter table by normalized enrichment score (NES):", 
                         choices = c("all", "positive","negative"),
                         selected = "all"
                       ),
                       # button to download FGSEA results
                       downloadButton("download_table", "Download filtered table")
                     ),
                     mainPanel(
                       DT::dataTableOutput("fgsea_table")
                     )
            ),
            tabPanel('Scatter plot',
                     sidebarPanel(
                       # filter scatter plot by adjusted p-value
                       sliderInput("padj_filter_scatter", "Filter plot by adjusted p-value:", min = 0, max = 1, value = 0.05)
                     ),
                     mainPanel(
                       plotOutput("fgsea_scatter")
                     )
            )
          )
        )
      )
    )
  )
)

server <- function(input, output) {
  # increasing size of file upload
  options(shiny.maxRequestSize=30*1024^2)
  
  # tab 1 : sample info
  # loading in sample info
  load_sample_info <- reactive({
    # requires file
    req(input$sample_info)
    # reads in csv file setting first row as header
    samples <- read.csv(input$sample_info$datapath, header=TRUE, stringsAsFactors = TRUE)
    # if data is tab delimited instead of comma-separated, it will load in as one big column
    if (ncol(samples) == 1) {
      print('Error: sample information data has 1 column. Ensure data is comma-separated with a header.')
    }
    return(samples)
  })
  # sample info summary table
  output$sample_info_summary <- renderTable({
    # load sample info data
    sample_info <- load_sample_info()
    # create data frame to show in table
    data_type = sapply(sample_info, class)
    # determine values of each column
    values <- sapply(seq_along(sample_info), function(i) {
      column <- sample_info[[i]]
      # if col is numeric type, calculate mean and sd
      if (is.numeric(column)) {
        paste0(round(mean(column, na.rm=TRUE), 2), " (+/- ", round(sd(column, na.rm=TRUE),2), ")")
      } else if (is.factor(column)) {
        levels <- levels(column)
        if (length(levels) > 3) {
          # if more than 3 levels, find common prefix and assign to be levels
          common_pre <- substr(levels, 1, 7)
          # use * to denote prefix
          paste(paste(common_pre[1], '*'))
        } else {
          # if levels <= 3, just display levels
          paste(paste(levels, collapse = ", "))
        }
      } else {
        NA
      }
    })
    # create data frame to show in table
    df <- data.frame(
      Column_Name = names(data_type),
      Type = as.character(data_type),
      Mean_or_Distinct_Values = values
    )
    # change formatting of col names
    colnames(df) <- c('Column Name', 'Type', 'Mean (Sd) or Distinct Values')
    return(df)
  })
  # use DT package to show sortable and searchable data table
  output$sample_info_table <- DT::renderDataTable({
    datatable(as.data.frame(load_sample_info()), 
              class = 'cell-border stripe',
              filter = 'top'
    )
  })
  # generate selectInput dynamically
  output$variable_ui <- renderUI({
    # load sample info
    sample_info <- load_sample_info()
    # generate numeric variable column names
    numeric_variables <- names(sample_info)[sapply(sample_info, is.numeric)] 
    # if column has a sd of 0 (like AvgSpotLen and version), it should be treated as
    # categorical / not continuous -> doesn't need histogram
    numeric_variables <- numeric_variables[sapply(sample_info[numeric_variables], sd, na.rm = TRUE) != 0]
    # generate select input
    selectInput(
      inputId = "variable", 
      label = "Choose a Variable:", 
      choices = numeric_variables,
      selected = if (length(numeric_variables) > 0) numeric_variables[1] else NULL
    )
  })
  # display histograms of continuous variables
  output$sample_info_histogram <- renderPlot({
    # get numeric variables from sample info
    sample_info <- load_sample_info()
    # assign x variable and axis label
    x_variable <- input$variable
    x_label <- x_variable
    # if variable is one with very large values, plot by billions or millions
    if (x_variable %in% c('Bases','Bytes')) {
      sample_info[[x_variable]] = sample_info[[x_variable]] / 1000000000
      # update x label
      x_label <- paste(x_variable, 'in Billions')
    }
    if (x_variable == 'mrna.seq_reads') {
      sample_info[[x_variable]] = sample_info[[x_variable]] / 1000000
      # update x label
      x_label <- paste(x_variable, 'in Millions')
    }
    # generate histogram plot
    ggplot(sample_info, aes_string(x = x_variable)) +
      geom_histogram(binwidth = as.numeric(input$binwidth), fill = "skyblue", color = "black") + 
      labs(title = paste("Histogram of", x_variable), x = x_label, y = "Frequency") +
      theme_minimal()
  })
  
  # tab 2 : counts
   load_norm_counts <- reactive({
     # requires file 
     req(input$norm_counts)
     # reads in csv file setting first col and row as gene and sample names
     counts <- read.csv(input$norm_counts$datapath, row.names=1, header=TRUE)
     if (ncol(counts) == 0) {
       print('Error: normalized count data has 0 columns. Ensure data is comma-separated with a header (sample names) and row names (genes).')
     }
     return(counts)
   })
   filter_norm_counts <- reactive({
     # load norm count data
     data <- load_norm_counts()
     # define filter-passing thresholds
     percentile <- input$var_percentile
     min_nonzeros <- input$nonzero_samples
     # calculate variances, median count, and number of nonzero samples
     variances <- apply(data, 1, var)
     median_counts <- apply(data, 1, median)
     variance_threshold <- quantile(variances, probs = percentile / 100)
     nonzero_counts <- rowSums(data > 0)
     num_zeros = rowSums(data == 0)
     # determine if genes pass filter
     passes_filter <- (variances >= variance_threshold) & (nonzero_counts >= min_nonzeros)
     # add to data frame
     filtered_data <- cbind(
       data,
       pass_filter = passes_filter,
       variance = variances,
       median_count = median_counts,
       nonzeros = nonzero_counts,
       zeros = num_zeros
     )
     return(as.data.frame(filtered_data))
   })
   # summary table output
   output$norm_counts_summary <- renderTable({
     # load all data and just filtered data
     counts <- filter_norm_counts()
     filtered <- counts[counts$pass_filter, ]
     # remove non-count columns from data
     counts <- counts %>% dplyr::select(-pass_filter, -variance, -median_count, -nonzeros, -zeros)
     filtered <- filtered %>% dplyr::select(-pass_filter, -variance, -median_count, -nonzeros, -zeros)
     # calculate percentages
     percent_pass = round((nrow(filtered) / nrow(counts))*100, 2)
     pass = paste('(', percent_pass, '%)', sep='')
     percent_fail = round(((nrow(counts)-nrow(filtered)) / nrow(counts))*100,2)
     fail = paste('(', percent_fail, '%)', sep='')
     # create data frame to show in table
     data.frame(
       Metric = c("Total Samples", "Total Genes", "Genes Passing Filter", "Genes Failing Filter"),
       Value = c(
         ncol(counts), 
         nrow(counts),
         paste(nrow(filtered), pass, sep=' '),
         paste((nrow(counts)-nrow(filtered)), fail, sep=' ')
       )
     )
   })
   # diagnostic scatter plots:
   # median vs variance (consider log scale for plot)
   output$variance_plot <- renderPlot({
     # data <- filter_stats()
     data <- filter_norm_counts()
     scatter <- ggplot(data, aes(x=median_count, y=variance, color=pass_filter)) +
       geom_point() + theme_minimal() +
       scale_color_manual(values = c('gold','blue'),labels=c("Does not pass filter", "Passes filter")) +
       labs(title='Scatterplot of Median Count vs Variance', x='Median Count', y='Variance')
     # scale variance and median with log10 if specified by input
     if (input$scale_var) {
       scatter <- scatter + scale_y_log10() + scale_x_log10() +
         labs(y = "log10(Variance)", x = 'log10(Median Count)')
     }
     return(scatter)
   })
   # median count vs number of zeros
   output$nonzero_plot <- renderPlot({
     data <- filter_norm_counts()
     scatter <- ggplot(data, aes(x=median_count, y=zeros, color=pass_filter)) +
       geom_point() + theme_minimal() +
       scale_color_manual(values = c('gold','blue'),labels=c("Does not pass filter", "Passes filter")) +
       labs(title='Scatterplot of Median Count vs Number of Zeros', x='Median Count', y='Number of Zeros')
     # scale median with log10 if specified by input
     if (input$scale_var) {
       scatter <- scatter + scale_x_log10() +
         labs(x = 'log10(Median Count)')
     }
     return(scatter)
   })
   # clustered heatmap of counts remaining after filtering w/ log transformed counts
   output$counts_heatmap <- renderPlot({
     # load filtered normalized count data
     counts <- filter_norm_counts()
     counts <- counts[counts$pass_filter, ]
     counts <- counts %>% 
       dplyr::select(-pass_filter, -variance, -median_count, -nonzeros, -zeros)
     # log transform counts
     log_counts <- log10(counts + 0.1)
     # convert to matrix for pheatmap()
     log_counts_matrix <- as.matrix(log_counts)
     # generating palette from RColorBrewer
     palette <- brewer.pal(n = 11, name = 'RdBu')
     # creating heatmap for signficantly differentially expressed probes
     pheatmap(log_counts_matrix,
              cluster_rows = FALSE,
              cluster_cols = TRUE,
              scale = 'none',
              legend = TRUE,
              show_colnames = TRUE,
              show_rownames = FALSE,
              main = 'Heatmap of Log-Transformed Normalized Counts After Filtering',
     )
   })
   # scatter plot of PCA projections
   output$counts_pca <- renderPlot({
     # load filtered normalized count data (samples as rows)
     counts <- filter_norm_counts() %>% 
       dplyr::select(-pass_filter, -variance, -median_count, -nonzeros, -zeros)
     # convert to matrix to run pca
     counts_matrix <- as.matrix(counts)
     pca <- prcomp(counts_matrix, scale. = TRUE)
     # create data frame of first 10 principal components
     pca_data <- as.data.frame(pca$x[, 1:10])
     colnames(pca_data) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
     # calculate indexes for input PCs by extracting number
     x_index <- as.numeric(sub("PC", "", input$x_pc)) # Extract number from "PC3" -> 3
     y_index <- as.numeric(sub("PC", "", input$y_pc)) # Extract number from "PC5" -> 5
     # calculate percent variance explained for input PCs
     percentVar <- round(100 * (pca$sdev[c(x_index, y_index)]^2 / sum(pca$sdev^2)), 2)
     # plot the PCA using ggplot2
     ggplot(pca_data, aes_string(x = input$x_pc, y = input$y_pc)) +
       geom_point(color="blue") + theme_minimal() +
       labs(title = 'PCA of Normalized Counts',
            x = paste0("PC1: ", percentVar[1], "% variance"),
            y = paste0("PC2: ", percentVar[2], "% variance"))
   })
   
   # tab 3 : differential expression (DE)
   load_raw_counts <- reactive({
     # requires file 
     req(input$raw_counts)
     # reads in raw counts csv file and sets row names and header
     counts <- read.csv(input$raw_counts$datapath, row.names=1, header=TRUE)
     if (ncol(counts) == 0) {
       print('Error: raw count data has 0 columns. Ensure data is comma-separated with a header (sample names) and row names (genes).')
     }
     # convert counts to matrix format for DESeq2
     counts <- as.matrix(counts)
     return(counts)
   })
   load_metadata <- reactive({
     # requires file 
     req(input$metadata)
     # reads in sample info csv file and sets header and factors
     metadata <- read.csv(input$metadata$datapath, header = TRUE, stringsAsFactors = TRUE)
     if (ncol(metadata) == 0) {
       print('Error: metadata has 0 columns. Ensure data is comma-separated with a header and row names.')
     }
     # only keep metadata needed for DESeq2: samples and diagnosis
     metadata <- metadata %>% dplyr::select(Sample.Name, diagnosis)
     # set 'diagnosis' to be a factor
     #metadata$diagnosis <- as.factor(metadata$diagnosis)
     return(metadata)
   })
   # creating reactive values to store DE results to use in GSEA
   DE_results <- reactiveVal(NULL)
   run_DESeq2 <- reactive({
     # load counts matrix and coldata
     cts <- load_raw_counts()
     coldata <- load_metadata()
     # DESeq2 suggests 10 for count filter
     count_filter = 10
     # pre-filtering based on count_filter
     keep <- rowSums(cts) >= count_filter
     cts <- cts[keep, ]
     # setting up dds for DESeq
     dds <- DESeqDataSetFromMatrix(countData = cts,
                                   colData = coldata,
                                   design = ~ diagnosis)
     # run DESeq on dds
     dds <- DESeq(dds)
     res <- as.data.frame(results(dds))
     # match entrez gene ids to hgnc symbols
     # using the `hsapiens_gene_ensembl` dataset in biomart
     res$gene <- row.names(res)
     ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
     ID_tbl <- getBM(
       # fetching attributes "entrezgene_id" and "hgnc_symbol"
       attributes = c("entrezgene_id","hgnc_symbol"),
       filters = "entrezgene_id", 
       values = res$gene,
       mart = ensembl
     )
     # add gene symbols to res
     ID_tbl$entrezgene_id <- as.character(ID_tbl$entrezgene_id)
     res <- res %>%
       left_join(ID_tbl, by = c("gene" = "entrezgene_id"))
     # reorder tibble and save as df
     res <- res %>% 
       dplyr::select(gene, hgnc_symbol, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
     # store res df as reactive result
     DE_results(as.data.frame(res))

   })
   # use DT package to show sortable and searchable data table
   output$de_table <- DT::renderDataTable({
     req(run_DESeq2())
     datatable(DE_results(), 
               class = 'cell-border stripe',
               filter = 'top'
     )
   })
   output$volcano_plot <- renderPlot({
     res <- DE_results()
     # determining expression of genes
     log2fc_threshold <- input$volcano_log2FC
     padj_threshold <- input$volcano_p
     res$expression <- ifelse(res$padj <= padj_threshold & res$log2FoldChange >= log2fc_threshold, "UP",
                              ifelse(res$padj <= padj_threshold & res$log2FoldChange <= 0-log2fc_threshold, "DOWN", "NA"))
     # volcano plot
     ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = expression)) +
       geom_point() + theme_minimal() +
       labs(title = "Volcano Plot of Gene Expression",
            x = "log2FoldChange", y = "-log10(adjusted p-value)") +
       scale_color_manual(values = c("UP" = "blue", "DOWN" = 'gold')) +
       theme(legend.position = 'bottom')
   })
   
   # tab 4 : gene set enrichment analysis (GSEA)
   output$GSEA_explanation <- renderText({
     "Since GSEA uses DE results, run DE in the previous tab before continuing with GSEA."
   })
   run_fgsea <- reactive({
     req(input$gene_set)
     # first, rank genes in DE results
     res <- DE_results()
     # determine ranking metric from input
     rank_metric <- input$metric
     if (rank_metric == "-log10(p-value)") {
       res$rank_metric = -log10(res$pvalue)
     } else if (rank_metric == "-log10(p-adj.)") {
       res$rank_metric = -log10(res$padj)
     } else {
       res$rank_metric = res$log2FoldChange
     }
     # arrange and descending order and form vector
     res <- res %>% arrange(desc(rank_metric))
     gene_symbol <- c(res$hgnc_symbol)
     rank_metric <- c(res$rank_metric)
     ranked_genes <- setNames(rank_metric, gene_symbol)
     # next, run FGSEA comparing ranked genes to gene sets
     gene_sets <- fgsea::gmtPathways(input$gene_set$datapath)
     gene_sets_upper <- lapply(gene_sets, toupper)
     fgsea_results <- fgsea::fgsea(
       pathways = gene_sets_upper,
       stats = ranked_genes,
       minSize = 15,
       maxSize = 500
     )
     fgsea_results_df <- as.data.frame(fgsea_results) %>% arrange(padj)
     return(fgsea_results_df)
   })
   output$download_res <- downloadHandler(
     filename = function() {
       paste0("fgsea_results.csv")
     },
     content = function(file) {
       fgsea_results_df <- run_fgsea()
       # write out results to file after collapsing leading edge column into string
       fgsea_results_df$leadingEdge <- sapply(fgsea_results_df$leadingEdge, function(x) paste(x, collapse = ";"))
       write.csv(fgsea_results_df, file, row.names = FALSE)
       showNotification("FGSEA results saved to fgsea_results.csv", type = "message")
     }
   )
   load_fgsea_results <- reactive({
     # requires file 
     req(input$fgsea_results)
     req(input$upload)
     # reads in sample info csv file and sets header and factors
     res <- read.csv(input$fgsea_results$datapath, header = TRUE, stringsAsFactors = TRUE)
     if (ncol(res) == 1) {
       print('Error: GSEA data has 1 column. Ensure data is comma-separated with a header.')
     }
     return(res)
   })
   fgsea_table_filter <- reactive({
     data <- load_fgsea_results()
     # filter by significant padj
     padj_threshold <- input$padj_filter_table
     data <- data %>% filter(padj < padj_threshold)
     # change data in dat abased on path filter input
     if (input$path_filter == 'positive') {
       data <- data %>% filter(NES > 0)
     } else if (input$path_filter == 'negative') {
       data <- data %>% filter(NES < 0)
     }
     return(data)
   })
   output$fgsea_table <- DT::renderDataTable({
     datatable(fgsea_table_filter(), 
               class = 'cell-border stripe',
               filter = 'top'
     )
   })
   # download button for fgsea table
   output$download_table <- downloadHandler(
     filename = function() {
       paste0("filtered_fgsea_table.csv")
     },
     content = function(file) {
       fgsea_results_df <- fgsea_table_filter()
       # write out results to file after collapsing leading edge column into string
       fgsea_results_df$leadingEdge <- sapply(fgsea_results_df$leadingEdge, function(x) paste(x, collapse = ";"))
       write.csv(fgsea_results_df, file, row.names = FALSE)
       showNotification("Filtered table saved to filtered_fgsea_table.csv", type = "message")
     }
   )
   # determining top n pathways in fgsea results using padj
   top_pathways <- reactive({
     fgsea_results <- load_fgsea_results()
     # select top positive and negative NES pathways
     top_pos <- fgsea_results %>% filter(NES > 0) %>%
       arrange(padj) %>% slice_head(n = input$n_pathways)
     top_neg <- fgsea_results %>% filter(NES < 0) %>% 
       arrange(padj) %>% slice_head(n = input$n_pathways)
     # bind the selected pathways
     top_paths <- bind_rows(top_pos, top_neg) %>%
       mutate(pathway = str_replace_all(pathway, "_", " ")) %>%
       mutate(pathway = sapply(pathway, function(x) paste(strwrap(x, width = 75), collapse = "\n")))
   })
   # barplot of fgsea NES for top pathways
   output$fgsea_barplot <- renderPlot({
     top_paths <- top_pathways()
     ggplot(top_paths, aes(x = reorder(pathway, NES), y = NES, fill = NES > 0)) +
       geom_bar(stat = "identity") + scale_fill_manual(values = c("TRUE" = "gold", "FALSE" = "blue")) +
       theme_minimal() + coord_flip() + 
       labs(y = "Normalized Enrichment Score (NES)", x = '',
            title = "FGSEA Results") + 
       theme(legend.position = 'none', 
             plot.title = element_text(size = 14, hjust = 0),
             axis.text.y = element_text(size = 6))
   })
   output$fgsea_scatter <- renderPlot({
     res <- load_fgsea_results()
     # determining expression of genes
     padj_threshold <- input$padj_filter_scatter
     res$enrichment <- ifelse(res$padj <= padj_threshold & res$NES > 0, "POSITIVE",
                              ifelse(res$padj <= padj_threshold & res$NES < 0, "NEGATIVE", "NA"))
     # volcano plot
     ggplot(res, aes(x = NES, y = -log10(padj), color = enrichment)) +
       geom_point() + theme_minimal() +
       labs(title = "Scatterplot of NES vs Adjusted p-Value",
            x = "Normalized Enrichment Score (NES)", y = "-log10(Adjusted P-Value)") +
       scale_color_manual(values = c("POSITIVE" = "gold", "NEGATIVE" = 'blue')) +
       theme(legend.position = 'bottom')
   })
}

# Run the application 
shinyApp(ui = ui, server = server)