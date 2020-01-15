# Splitting the data into training and test sets
# Keaton Wilson
# keatonwilson@me.com
# 2019-11-25

# packages
require(tidyverse)


#' train_test_split
#'A small function that separates a data set of occurence and environmental data
#'into training and test sets
#' @param extra_prepped_data input data generated from the 
#' \code{\link{prep_data_2}} function.
#' @param blocked_obj a blockCV object generated from the 
#' \code{\link{run_block_cv}} function
#'
#' @return a list two items long, the first being the training data and the 
#' second being the test data
#' @export
#'
#' @examples
train_test_split = function(extra_prepped_data, blocked_obj){
  
  extract_index = function(list_of_folds = NULL) {
    for(k in 1:length(list_of_folds)){
      train_index <- unlist(list_of_folds[[k]][1]) # extract the training set indices
      test_index <- unlist(list_of_folds[[k]][2])# extract the test set indices
    }
    mini_list = list(train_index, test_index)
    return(mini_list)
  }
  
  indices = extract_index(blocked_obj$folds)
  print(length(indices[[1]]))
  print(length(indices[[2]]))
  
  #applying indexes and splitting data
  training_data = extra_prepped_data[indices[[1]],] %>%
    drop_na()
  test_data = extra_prepped_data[-indices[[2]],] %>%
    drop_na()
  
  return(list(training_data = training_data, 
              test_data = test_data))
}

