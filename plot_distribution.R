plot_distribution <- function(micro_array_data_file)
{

  library(ggplot2)

  output_directory <- get_output_directory()

  #read the data from the file into a data frame using a user-defined function
  DATA <- read_txt_to_df(micro_array_data_file)

  DATA <- log2(DATA)

  DATA_long <- FatToLongDF(DATA)


  XLabel <- 'log2(Intensity)'
  YLabel <- 'Number of Probes'
  ggplot(DATA_long, aes_string(x='value')) +
      geom_histogram(binwidth=1) +
      labs(x = XLabel) +
      labs(y = YLabel)

  pdf_width <- 4
  pdf_height <- 4
  ggsave(paste(output_directory,'probe_distribution.pdf',sep=''), width = pdf_width, height = pdf_height, dpi = 300, limitsize=FALSE)

  return()

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
