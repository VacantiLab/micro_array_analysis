plot_distribution <- function(micro_array_data_file)
# this function plots the distribution of log2 transformed probe intensities
#     the input is a micro-array tab-delimited table
#        the probe IDs are down the rows and the sample names across the columns
#        the raw table is usually processed by a separate python script

{

  library(ggplot2)

  output_directory <- get_output_directory()

  #read the data from the file into a data frame using a user-defined function
  DATA <- read_txt_to_df(micro_array_data_file)

  #record the number of samples
  n_samples <- length(colnames(DATA))

  #record the number of probes present before filtering
  n_probes <- length(rownames(DATA))

  #log2 transform the data
  DATA <- log2(DATA)

  #filter the probes
  DATA_filtered <- DATA
  filter_threshold <- 5
  filter_criteria <- 0.2
  counter <- 1
  #probes_to_delete <- NULL
  probes_to_delete <- rep(FALSE,n_probes)
  for (probe in rownames(DATA))
  {
      array_of_interest <- DATA_filtered[probe,]
      samples_meeting_criteria <- array_of_interest >= filter_threshold
      n_meeting_threshold <- sum(samples_meeting_criteria)
      fraction_meeting_threshold <- n_meeting_threshold/n_samples
      criteria_met <- fraction_meeting_threshold >= filter_criteria
      if (!criteria_met)
      {
          #probes_to_delete <- c(probes_to_delete, probe)
          probes_to_delete[counter] <- TRUE
      }

      if (counter %% 100 == 0)
      {
          print(paste('probe being filtered: ',counter))
      }
      counter = counter + 1
  }

  #n_probes_to_delete <- length(probes_to_delete)
  n_probes_to_delete <- sum(probes_to_delete)
  n_probes_filtered <- n_probes - n_probes_to_delete

  #create the filtered data frame
  #rows_to_delete <- rownames(DATA_filtered) %in% probes_to_delete
  #rows_to_keep <- !rows_to_delete
  rows_to_keep <- !probes_to_delete
  DATA_filtered <- DATA_filtered[rows_to_keep,]


  #put the data in long format
  #    columns correspond to probe ID, intensity, and sample ID
  print('putting the data in long data frame format')
  DATA_long <- FatToLongDF(DATA)
  DATA_filtered_long <- FatToLongDF(DATA_filtered)

  print('plotting the intensity distributions')
  #plot the distribution of the probe density in the data
  XLabel <- 'log2(Intensity)'
  YLabel <- 'Number of Probes'
  ggplot(DATA_long, aes_string(x='value')) +
      #geom_histogram(binwidth=1) +
      labs(x = XLabel) +
      labs(y = YLabel) +
      geom_density(alpha=0.2,fill='#FF6666')

  pdf_width <- 4
  pdf_height <- 4
  ggsave(paste(output_directory,'data_probe_distribution.pdf',sep=''), width = pdf_width, height = pdf_height, dpi = 300, limitsize=FALSE)

  #plot the distribution of the probe density in the filtered data
  XLabel <- 'log2(Intensity)'
  YLabel <- 'Number of Probes'
  ggplot(DATA_filtered_long, aes_string(x='value')) +
      #geom_histogram(binwidth=1) +
      labs(x = XLabel) +
      labs(y = YLabel) +
      geom_density(alpha=0.2,fill='#FF6666')

  pdf_width <- 4
  pdf_height <- 4
  ggsave(paste(output_directory,'data_filtered_probe_distribution.pdf',sep=''), width = pdf_width, height = pdf_height, dpi = 300, limitsize=FALSE)

  #print('collapsing the filtered probes to gene symbols')
  #DATA_gene_symbols <- probes_to_genes(data_set)

  #transform data back to untransformed
  DATA <- 2^DATA
  DATA_filtered <- 2^DATA_filtered

  #to_return <- list(DATA,DATA_filtered,DATA_gene_symbols)
  to_return <- list(DATA,DATA_filtered)
  return(to_return)

}

###############################################################################

read_txt_to_df <- function(txt_directory)
#Function to read a tab delimited text file into a data frame
{
    DF <- read.table(file=txt_directory,head=TRUE,check.names=FALSE,sep='\t') #check.names=FALSE prevents an 'X' from being added to the numeric column names
    #Name the rows of the data frame as the genes given in the first column of the data frame
    RowNames <- as.character(DF[,1])
    rownames(DF) <- RowNames
    DF[,1] <- NULL #remove the first column of the data frame as it is no longer needed
    return(DF)
}

###############################################################################

get_output_directory <- function()
{
    #Input where the heatmap is stored
    working_directory <- getwd()
    working_directory_up1 <- gsub('/[^/]*$','/',working_directory) #matches '/' followed by 0 or more characters other than '/' followed by the end of the string, and replaces with '/'
    output_directory <- paste(working_directory_up1,'output/',sep='')
    if (!file.exists(output_directory)){dir.create(output_directory)} #creates the directory for the output
    return(output_directory)
}

###############################################################################

FatToLongDF <- function(DATA)
{
    n_row_DATA <- nrow(DATA)
    n_col_DATA <- ncol(DATA)
    DATA_long_len <- nrow(DATA)*ncol(DATA)
    DATA_long <- data.frame(matrix(nrow=DATA_long_len,ncol=3))
    colnames(DATA_long) <- c('probe','sample','value')
    probe_names <- rownames(DATA)
    sample_names <- colnames(DATA)

    for (i in 1:n_col_DATA)
    {
        current_sample <- sample_names[i]
        start_iteration <- (i-1)*n_row_DATA + 1
        end_iteration <- i*n_row_DATA
        DATA_long[start_iteration:end_iteration,'probe'] <- rownames(DATA)
        DATA_long[start_iteration:end_iteration,'sample'] <- rep(current_sample,n_row_DATA)
        DATA_long[start_iteration:end_iteration,'value'] <- DATA[,current_sample]
    }

    return(DATA_long)
}
