install.packages('rjson', dependencies = TRUE)

library("rjson")

data <- fromJSON(file = "C:\\Users\\Daniel\\OneDrive - Erasmus University Rotterdam\\Documents\\Master\\Computer Science\\TVs-all-merged.json")

keys <- names(data)
n <- length(keys)

# Function to make a qGram from the title
qGram <- function(title, q) {
  qGram <- c()
  for(i in 1: (nchar(title)-q+1)){
    qGram <- c(qGram,substring(title,i,q+i-1))
  }
  return(qGram)
}

# Product representation
q <- 5
#Create data list with the observations not sorted by their ModelID
dataList2 <- list(list(title=data[[1]][[1]]$title, 
                       modelID=data[[1]][[1]]$modelID,
                       shop = data[[1]][[1]]$shop,
                      # Remove all no alphanumeric tokens from the title to construct modelwords and qGram product representations
                       title2 = gsub("[^[:alnum:]|^.]", "",gsub("[\"]", "inch", tolower(data[[1]][[1]]$title))), 
                       modelwords=c(unlist(strsplit(tolower(gsub("[^[:alnum:]|^.]", " ",gsub("[\"]", "inch", data[[1]][[1]]$title))), split= " "))),
                       qGram = qGram(gsub("[^[:alnum:]|^.]", "",gsub("[\"]", "inch", tolower(data[[1]][[1]]$title))), q)
                       ))
for(i in 2:length(keys)){
  m = length(data[[i]])
  for(j in 1:m){
    title_ij <- data[[i]][[j]]$title
    shop_ij <- data[[i]][[j]]$shop
    modelID_ij <- data[[i]][[j]]$modelID
    obs2 <- list(list(title=title_ij, 
                     modelID=modelID_ij,
                     shop = shop_ij,
                     title2=gsub("[^[:alnum:]|^.]", "",gsub("[\"]", "inch", tolower(data[[i]][[j]]$title))),
                     modelwords=c(unlist(strsplit(tolower(gsub("[^[:alnum:]|^.]", " ",gsub("[\"]", "inch", title_ij))), split= " "))),
                     qGram = qGram(gsub("[^[:alnum:]|^.]", "",gsub("[\"]", "inch", tolower(data[[i]][[j]]$title))), q)
                      ))
    dataList2 <- append(dataList2, obs2)
  }
}

# Create a set of all model words
modelwords2 <- dataList2[[1]]$modelwords
for (i in 1:length(dataList2)){
  modelwords2 <- c(modelwords2, dataList2[[i]]$modelwords)
}

modelwords2 <- unique(modelwords2[nchar(modelwords2)>=1])
# There are 1653 unique model words 

# Similarity
# Create a sparse matrix of binary product representations
# Every binary indicator gives whether a model word is in the title of a product
sparse <- matrix(0,nrow=length(modelwords2), ncol=length(dataList2))
for (i in 1:length(modelwords2)) {
  for (j in 1:length(dataList2)) {
    if (modelwords2[i] %in% dataList2[[j]]$modelwords){
      sparse[i,j]=1
    }
  }
}

# Minhashing
# Convert the matrix with binary product representations into a signature matrix
minhash <- function(matrix, n_permutations) {
  set.seed(0)
  signature <- matrix(nrow=n_permutations, ncol=ncol(matrix))
  for (j in 1:n_permutations) {
    permutation <-matrix[sample(nrow(matrix)),]
    row <- c()
    for (i in 1:ncol(matrix)) {
      index <- which(permutation[,i]!=0)[1]
      row <- c(row,index)
    }
    signature[j,]=row
  }
  return(signature)
}

# Reduce the matrix of binary vector representations to a signature matrix 
# Such that it is half the dimension of the binary vector matrix
n <- nrow(sparse)
signature4 <- minhash(sparse,floor(0.5*n)) 
# reduce the signature to 826 x 1624 matrix


# LHS # CHECK
# The factors of 1624 are:
factors <- c(1, 2, 7, 14, 59, 118, 413, 812)
threshold <- matrix(nrow = length(factors)+1, ncol = length(factors)+1)
threshold[1,]<-c(0,factors)
threshold[,1]<-c(0,factors)
for (i in 1:length(factors)){
    threshold[(length(factors)+2-i),(i+1)]=1/tail(factors,i)[1]^(1/factors[i])
}
threshold

# Performs locality sensitivity hashing on a signature matrix
lsh_new <- function(m, r, b){
  n <- ncol(m)
  potential <- diag(0,nrow=n)
  for (k in 1:(n-1)){
    for (l in (k+1):n) {
      for (j in 1:b) {
        if (potential[k,l]!=1){
          if(identical(m[(1+r*(j-1)):(r*(j-1)+r),k],m[(1+r*(j-1)):(r*(j-1)+r),l])){
            potential[k,l]=1
            potential[l,k]=1
          }
        }
      }
      if(dataList2[[k]]$shop==dataList2[[l]]$shop){
        potential[k,l]=0
        potential[l,k]=0
      }
    }
  }
  return(potential)
}

c1 <- lsh_new(signature4,59,14)
sum(c1)/2 # 0 comparisons, more bands required

c2 <- lsh_new(signature4,14,59)
sum(c2)/2 # 24, more bands required
t_c2 <- (1/59)^(1/14)

c3 <- lsh_new(signature4,7,118)
sum(c3)/2 #  4030, lot of comparisons
t_c3 <- (1/118)^(1/7)

# To test scalability use both c2 and c3

# High jump for b going from 29 bands to 58, run for both to show trade-off between pair quality and pair completeness

# Take the two potential comparisons matrix and calculate the dissimilarity matrices
# With the dissimilarity matrices perform the bootstrap to find the optimal value of the threshold
# Use the optimal value of the threshold on the test set
# Calculate the scores of the test set

# Calculate the similarity between two products using the overlapping q-grams
qGramDistance <- function(q1, q2){
  n1 <- length(q1)
  n2 <- length(q2)
  counter <- 0
  for (i in 1:n1){
    if (!(q1[i] %in% q2)){
      counter <- counter + 1
    }
  }
  for (j in 1:n2){
    if (!(q2[j] %in% q1)){
      counter <- counter + 1
    }
  }
  return((n1+n2-counter)/(n1+n2))
}

# Construct a dissimilarity matrix for only the comparisons provided by LSH
dissimilarityMatrix <- function(data, comparisons){
  N <- length(data)
  dissimilarity <- 1-comparisons
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      if(dissimilarity[i,j]==0) {
        dissimilarity[i,j]=1-qGramDistance(data[[i]]$qGram,data[[j]]$qGram)
        dissimilarity[j,i]=dissimilarity[i,j]
      }
    }
  }
  return(dissimilarity)
}

# Construct a matrix indicating which products are identical
N <- length(dataList2)
identical <- diag(0, nrow=N)
for (i in 1:(N-1)){
  for (j in (i+1):N){
    if (dataList2[[i]]$modelID==dataList2[[j]]$modelID) {
      identical[i,j]=1
      identical[j,i]=1
    }
  }
}
sum(identical)/2 # 399 duplicates

# EVALUATION

# Run a bootstrap for different threshold levels
runBootstrap <- function(b, d,threshold_values){
  # Set seed to ensure same samples
  set.seed(1)
  n <- nrow(d)
  # Construct matrix to store results
  results <- matrix(nrow=b, ncol=5)
  for (i in 1:b) {
    results_b <- data.frame()
    # Generate bootstrap samples
    train_idx <- sample(1:n,replace = TRUE)
    test_idx <- c(1:n)[-train_idx]
    best_score <- 0
    best_threshold <- threshold_values[1]
    # Construct dissimilarity matrix and identical matrix for training set
    dissimilarity_train <- d[train_idx,train_idx]
    identical_train <- identical[train_idx,train_idx]
    # Use training set the optimal parameters by using a grid search
    for (epsilon in threshold_values) {
      # Construct matrix of pairs which are labelled as duplicates 
      duplicates_train <- dissimilarity_train<=epsilon
      TP <- 0
      TN <- 0
      FP <- 0
      FN <- 0
      for (k in 1:(n-1)){
        for (j in (k+1):n){
          if(identical_train[k,j]==1 & duplicates_train[k,j]==1){
            TP <- TP + 1
          }
          if(identical_train[k,j]==1 & duplicates_train[k,j]==0){
            FN <- FN + 1
          }
          if(identical_train[k,j]==0 & duplicates_train[k,j]==1){
            FP <- FP + 1
          }
          if(identical_train[k,j]==0 & duplicates_train[k,j]==0){
            TN <- TN + 1
          }
        }
      }
      percision <- TP/(TP+FP)
      recall <- TP/(TP+FN)
      score <- (2*percision*recall)/(percision+recall)
      if(is.nan(score)){
        score <- 0
      }
      # Update score and threshold level if score is improved 
      if(score>best_score){
        best_score <- score
        best_threshold <- epsilon
      }
    }
    # Construct dissimilarity, identical and duplicate matrix for test set
    dissimilarity_test <- d[test_idx,test_idx]
    identical_test <- identical[test_idx,test_idx]
    # Within bootstrap choose optimal threshold value with the highest F1 score
    duplicates_test <- (dissimilarity_test <= best_threshold)
    TP <- 0
    TN <- 0
    FP <- 0
    FN <- 0
    n <- length(test_idx)
    for (l in 1:(n-1)){
      for (m in (l+1):n){
        if(identical_test[l,m]==1 & duplicates_test[l,m]==1){
          TP <- TP + 1
        }
        if(identical_test[l,m]==1 & duplicates_test[l,m]==0){
          FN <- FN + 1
        }
        if(identical_test[l,m]==0 & duplicates_test[l,m]==1){
          FP <- FP + 1
        }
        if(identical_test[l,m]==0 & duplicates_test[l,m]==0){
          TN <- TN + 1
        }
      }
    }
    # Calculate results for the test set
    percision_test <- TP/(TP+FP)
    if(is.nan(percision_test)){
      percision_test <- 0
    }
    recall_test <- TP/(TP+FN)
    if(is.nan(recall_test)){
      recall_test <- 0
    }
    score_test <- (2*percision_test*recall_test)/(percision_test+recall_test)
    if(is.nan(score_test)){
      score_test <- 0
    }
    num_comparisons <- sum(dissimilarity_test[dissimilarity_test!=1])/2
    num_duplicates <- sum(identical_test)/2
    pair_quality <- TP/num_comparisons
    if(is.nan(pair_quality)){
      pair_quality <-0
    }
    pair_completeness <- TP/num_duplicates
    if(is.nan(pair_completeness)){
      pair_completeness <- 0
    }
    f1_star <- (2*pair_quality*pair_completeness)/(pair_quality+pair_completeness)
    if(is.nan(f1_star)){
      f1_star <-0
    }
    results_b <- c(PC = pair_completeness, PQ = pair_quality, F_1_star = f1_star, F_1 = score_test, threshold = best_threshold)
    results[i,] <- results_b
  }
  return(results)
}

# Run over the 5 bootstraps
b <- 5 
# Set a grid search for the threshold value
threshold_values <- c(1:9/10)
# First run for dissimilarities with few comparisons
d1 <- dissimilarityMatrix(dataList2,c2)
results <- runBootstrap(b,d1, threshold_values)
avg_results <- colMeans(results)
avg_results
# Second run with a broader dissimilarity matrix 
d2 <- dissimilarityMatrix(dataList2,c3)
results2 <- runBootstrap(b,d2, threshold_values)
avg_results2 <- colMeans(results2)
avg_results2

# Maybe extra comparisons
c4 <- lsh_new(signature4,2,413)
t_c4 <- (1/413)^(1/2)
sum(c4)/2 # 628780
d3 <- dissimilarityMatrix(dataList2,c4)
results3 <- runBootstrap(b,d3,threshold_values)
avg_results3 <- colMeans(results3)