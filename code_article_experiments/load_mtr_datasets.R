
load_mtr_datasets <- function(dataset){

  # list of datasets
  datasets_list <- c("enb", "jura", "wq", "scm20d", "rf1", "oes97",
                     "air", "births2", "wage")
  
  # check dataset exists
  if (!dataset %in% datasets_list){
    stop('Unknown dataset.')
  }
  
  if (dataset %in% c("air", "births2", "wage")){
    
    if (dataset == "air"){
      load("mtr_datasets/air_data_benchmark.Rdata")
    }
    if (dataset == "births2"){
      load("mtr_datasets/births_benchmark2.Rdata")
    }
    if (dataset == "wage"){
      load("mtr_datasets/wage_benchmark.Rdata")
    }
    p <- ncol(X)
    
  }else{
  
    # define number of inputs
    input_number <- read.csv('mtr_datasets/input_number.csv', header = F)
    p <- input_number[input_number[, 1] == dataset, 2]
  
    # load dataset
    data <- read.arff(paste0('mtr_datasets/', dataset, '.arff'))
    data <- data[!apply(data, 1, anyNA),]
    X <- as.matrix(data[1:p])
    Y <- as.matrix(data[(p+1):ncol(data)])
    
  }

  return(list(X = X, Y = Y, p = p))
  
}
