#' SRS
#' 
#' Normalize species counts by scaling with ranked sub-sampling
#' 
#' @param data an OTU table with samples in columns and taxa in rows; must be supplied as a data frame rather than a matrix
#' @param Cmin the desired counts per sample; by default, the minimum sample count in the submitted OTU table
#'
#' @return - an OTU table with equal sample counts
#' 
#' @author  L. Beule and P.Karlovsky
#' 
#' @references Beule L, Karlovsky P. Improved normalization of species count data in ecology by scaling with ranked sub-sampling (SRS): application to microbial communities. PeerJ. 2020;8:e9593.
#' 
SRS <- function(data, Cmin){
  if(Cmin > min(colSums(data))){
    print("ERROR: Cmin > minimum library size. Please select a Cmin that is <= the minimum library size of the dataset.")
  } else {
    if(Cmin < 0){
      print("ERROR: Cmin < 0. Please select a Cmin >= 0.")
    } else {
      if(Cmin %% 1 > 0){
        print("ERROR: Please select a Cmin without decimal places.")
      } else {
        counter = 1 #counting the loops
        for(i in seq(1, ncol(data), 1)){
          if (i == 1) {
            fixed_factor <- (data[,i]/(sum(data[,i])/Cmin))
            assign(paste(names(data)[i],sep=""),fixed_factor)
            fixed_factor_1 <- data.frame(get(names(data)[i]))
            names(fixed_factor_1)[i] <- names(data)[i]
          } else {
            fixed_factor <- (data[,i]/(sum(data[,i])/Cmin))
            assign(paste(names(data)[i],sep=""),fixed_factor)
            fixed_factor_1 <- cbind(fixed_factor_1, fixed_factor)
            names(fixed_factor_1)[{
              counter = counter + 1
            }] <- names(data)[i]
          }
        }
        
        fixed_factor_1
        
        revtrunc_fixed_factor_1 <- floor(fixed_factor_1) # floor (e.g. 1.9 will become 1)
        revtrunc_fixed_factor_1
        
        diff_counts <- Cmin-colSums(revtrunc_fixed_factor_1) #how many counts differences to the selected library size?
        diff_counts
        
        revtrunc <- function(x) { sign(x) * (x - floor(x)) } 
        revtrunc_fixed_factor <- (round(revtrunc(fixed_factor_1),10000000))
        revtrunc_fixed_factor
        
        x <- as.data.frame(revtrunc_fixed_factor)
        counter = 1
        for(i in seq(1, ncol(x), 1)){
          if (i == 1) {
            if(diff_counts[i] == 0){
              fixed_factor <- revtrunc_fixed_factor_1[,i]
              assign(paste(names(data)[i],sep=""),fixed_factor)
              fixed_factor_1 <- data.frame(get(names(data)[i]))
              colnames(fixed_factor_1)[i] <- names(data)[i]
            } #if the sum of the counts in the library = Cmin
            else {
              maxN <- function(x, N=diff_counts[i]){
                len <- length(x)
                if(N>len){
                  warning('N greater than length(x).  Setting N=length(x)')
                  N <- length(x)
                }
                sort(x,partial=len-N+1)[len-N+1]
              }
              max <- which(x[,i] == maxN(x[,i]), arr.ind = TRUE)
              max
              sum(x[,i] > unique(x[,i][max]))
              normalization_value <- diff_counts[i] - sum(x[,i] > unique(x[,i][max]))
              normalization_value
              
              lowest_level_choise <- as.data.frame(which(x[,i] == unique(maxN(x[,i]))))
              lowest_level_choise 
              length(t(lowest_level_choise)) #how many counts have to be sampled?
              
              if(sum(revtrunc_fixed_factor_1[,i][lowest_level_choise[,1]]) == 0){
                lowest_level <- as.numeric(as.vector(sample(as.factor(lowest_level_choise[,1]), normalization_value, replace = F)))
                y <- as.vector(rep(0,length(x[,1]))) #create an empty vector
                y[lowest_level] = 1 #set the randomly selected counts to 1
                y
                
              } #if all of the integer values of the lowest rank = 0, do random subsmpling
              else {
                sub_int <- subset(lowest_level_choise, (revtrunc_fixed_factor_1[,i][lowest_level_choise[,1]] >= 1) == TRUE)
                sub_int_bind <- as.data.frame(cbind(sub_int,revtrunc_fixed_factor_1[,i][sub_int[,1]]))
                names(sub_int_bind)[1] <- "V1"
                names(sub_int_bind)[2] <- "V2"
                
                sub_int_bind_ordered <- sub_int_bind[order(sub_int_bind$V2, decreasing = TRUE),]  
                sub_int_bind_ordered
                
                sub_int_bind_ordered_V1 <- sub_int_bind_ordered$V1
                sub_int_bind_ordered_V2 <- sub_int_bind_ordered$V2
                
                if((length(unique(sub_int_bind_ordered_V2)) == 1 & length(sub_int_bind_ordered_V2)>as.vector(normalization_value))){
                  lowest_level <- as.numeric(as.vector(sample(as.factor(sub_int_bind_ordered_V1), normalization_value, replace = F)))
                  y <- as.vector(rep(0,length(x[,1]))) #create an empty vector
                  y[lowest_level] = 1 #set the randomly selected counts to 1
                  y
                } #if all of the integer values of the lowest rank are equal, do random subsampling
                else {
                  if(length(sub_int_bind_ordered_V1)>normalization_value){
                    
                    maxN_1 <- function(x, N=normalization_value){
                      len <- length(x)
                      if(N>len){
                        warning('N greater than length(x).  Setting N=length(x)')
                        N <- length(x)
                      }
                      sort(x,partial=len-N+1)[len-N+1]
                    }
                    max_1 <- which(as.data.frame(sub_int_bind_ordered_V2)[,1] == maxN_1(as.data.frame(sub_int_bind_ordered_V2)[,1]), arr.ind = TRUE)
                    max_1
                    sum(as.data.frame(sub_int_bind_ordered_V2)[,1] > unique(as.data.frame(sub_int_bind_ordered_V2)[,1][max_1]))
                    
                    normalization_value_1 <- normalization_value - sum(as.data.frame(sub_int_bind_ordered_V2)[,1] > unique(as.data.frame(sub_int_bind_ordered_V2)[,1][max_1])) # how may values are above the maxima that is the one to be normalized?
                    normalization_value_1
                    lowest_level_choise_1 <- as.data.frame(which(as.data.frame(sub_int_bind_ordered_V2)[,1] == unique(maxN_1(as.data.frame(sub_int_bind_ordered_V2)[,1]))))
                    lowest_level_choise_1
                    lowest_level <- as.numeric(as.vector(sample(as.factor(lowest_level_choise_1[,1]), normalization_value_1, replace = F)))
                    lowest_level <- sub_int_bind_ordered_V1[lowest_level]
                    lowest_level
                    
                    lowest_level_1 <- sub_int_bind_ordered_V1[(as.data.frame(sub_int_bind_ordered_V2)[,1] > unique(as.data.frame(sub_int_bind_ordered_V2)[,1][max_1]))]
                    
                    lowest_level <-c(lowest_level_1, lowest_level)
                    y <- as.vector(rep(0,length(x[,1]))) #create an empty vector
                    y[lowest_level] = 1 #set the randomly selected counts to 1
                    y
                    
                  } #if integer ranks are > normalization_value, do ranked subsampling
                  else {
                    if(length(sub_int_bind_ordered_V1)<normalization_value){
                      sub_int_zeros <- subset(lowest_level_choise, (revtrunc_fixed_factor_1[,i][lowest_level_choise[,1]] < 1) == TRUE)
                      length(t(sub_int_zeros))
                      
                      lowest_level_2 <- as.numeric(as.vector(sample(as.factor(sub_int_zeros[,1]), (normalization_value-length(sub_int_bind_ordered_V1)), replace = F)))
                      lowest_level_2
                      lowest_level_3 <- c(sub_int_bind_ordered_V1,lowest_level_2)
                      
                      y <- as.vector(rep(0,length(x[,1]))) #create an empty vector
                      y[lowest_level_3] = 1 #set the randomly selected counts to 1
                      y
                      
                    } #if integer ranks are < normalization_value, do ranked subsampling of the zero integers
                    else {
                      y <- as.vector(rep(0,length(x[,1]))) #create an empty vector
                      y[sub_int_bind_ordered_V1] = 1 #set the randomly selected counts to 1
                      y
                    } # if integer ranks are = normalization_value, sample all
                  }
                  
                  
                }
                
              } 
              
              SRS <- revtrunc_fixed_factor_1[,i] + ceiling(x[,i] > unique(x[,i][max])) + y #sum it all u
              SRS
              sum(SRS) #verification
              assign(paste(colnames(data)[i],sep=""),SRS)
              fixed_factor_1 <- data.frame(get(names(data)[i]))
              colnames(fixed_factor_1)[i] <- names(data)[i]
            } #if the sum of the counts in the library > Cmin
          } #for the first libraray
          else {
            if(diff_counts[i] == 0){
              fixed_factor <- revtrunc_fixed_factor_1[,i]
              assign(paste(colnames(data)[i],sep=""),fixed_factor)
              fixed_factor_1 <- cbind(fixed_factor_1, fixed_factor)
              colnames(fixed_factor_1)[{
                counter = counter + 1
              }] <- names(data)[i]
            } #if the sum of the counts in the library = Cmin
            else { 
              maxN <- function(x, N=diff_counts[i]){
                len <- length(x)
                if(N>len){
                  warning('N greater than length(x).  Setting N=length(x)')
                  N <- length(x)
                }
                sort(x,partial=len-N+1)[len-N+1]
              }
              max <- which(x[,i] == maxN(x[,i]), arr.ind = TRUE)
              max
              sum(x[,i] > unique(x[,i][max]))
              normalization_value <- diff_counts[i] - sum(x[,i] > unique(x[,i][max]))
              normalization_value
              
              lowest_level_choise <- as.data.frame(which(x[,i] == unique(maxN(x[,i]))))
              lowest_level_choise 
              length(t(lowest_level_choise)) #how many counts have to be sampled?
              
              if(sum(revtrunc_fixed_factor_1[,i][lowest_level_choise[,1]]) == 0){
                lowest_level <- as.numeric(as.vector(sample(as.factor(lowest_level_choise[,1]), normalization_value, replace = F)))
                y <- as.vector(rep(0,length(x[,1]))) #create an empty vector
                y[lowest_level] = 1 #set the randomly selected counts to 1
                y
                
              } #if all of the integer values of the lowest rank = 0, do random subsmpling
              else {
                sub_int <- subset(lowest_level_choise, (revtrunc_fixed_factor_1[,i][lowest_level_choise[,1]] >= 1) == TRUE)
                sub_int_bind <- as.data.frame(cbind(sub_int,revtrunc_fixed_factor_1[,i][sub_int[,1]]))
                names(sub_int_bind)[1] <- "V1"
                names(sub_int_bind)[2] <- "V2"
                
                sub_int_bind_ordered <- sub_int_bind[order(sub_int_bind$V2, decreasing = TRUE),]  
                sub_int_bind_ordered
                
                sub_int_bind_ordered_V1 <- sub_int_bind_ordered$V1
                sub_int_bind_ordered_V2 <- sub_int_bind_ordered$V2
                
                if((length(unique(sub_int_bind_ordered_V2)) == 1 & length(sub_int_bind_ordered_V2)>as.vector(normalization_value))){
                  lowest_level <- as.numeric(as.vector(sample(as.factor(sub_int_bind_ordered_V1), normalization_value, replace = F)))
                  y <- as.vector(rep(0,length(x[,1]))) #create an empty vector
                  y[lowest_level] = 1 #set the randomly selected counts to 1
                  y
                } #if all of the integer values of the lowest rank are equal, do random subsampling
                else {
                  if(length(sub_int_bind_ordered_V1)>normalization_value){
                    
                    maxN_1 <- function(x, N=normalization_value){
                      len <- length(x)
                      if(N>len){
                        warning('N greater than length(x).  Setting N=length(x)')
                        N <- length(x)
                      }
                      sort(x,partial=len-N+1)[len-N+1]
                    }
                    max_1 <- which(as.data.frame(sub_int_bind_ordered_V2)[,1] == maxN_1(as.data.frame(sub_int_bind_ordered_V2)[,1]), arr.ind = TRUE)
                    max_1
                    sum(as.data.frame(sub_int_bind_ordered_V2)[,1] > unique(as.data.frame(sub_int_bind_ordered_V2)[,1][max_1]))
                    
                    normalization_value_1 <- normalization_value - sum(as.data.frame(sub_int_bind_ordered_V2)[,1] > unique(as.data.frame(sub_int_bind_ordered_V2)[,1][max_1])) # how may values are above the maxima that is the one to be normalized?
                    normalization_value_1
                    lowest_level_choise_1 <- as.data.frame(which(as.data.frame(sub_int_bind_ordered_V2)[,1] == unique(maxN_1(as.data.frame(sub_int_bind_ordered_V2)[,1]))))
                    lowest_level_choise_1
                    lowest_level <- as.numeric(as.vector(sample(as.factor(lowest_level_choise_1[,1]), normalization_value_1, replace = F)))
                    lowest_level <- sub_int_bind_ordered_V1[lowest_level]
                    lowest_level
                    
                    lowest_level_1 <- sub_int_bind_ordered_V1[(as.data.frame(sub_int_bind_ordered_V2)[,1] > unique(as.data.frame(sub_int_bind_ordered_V2)[,1][max_1]))]
                    
                    lowest_level <-c(lowest_level_1, lowest_level)
                    y <- as.vector(rep(0,length(x[,1]))) #create an empty vector
                    y[lowest_level] = 1 #set the randomly selected counts to 1
                    y
                    
                  } #if integer ranks are > normalization_value, do ranked subsampling
                  else {
                    if(length(sub_int_bind_ordered_V1)<normalization_value){
                      sub_int_zeros <- subset(lowest_level_choise, (revtrunc_fixed_factor_1[,i][lowest_level_choise[,1]] < 1) == TRUE)
                      length(t(sub_int_zeros))
                      
                      lowest_level_2 <- as.numeric(as.vector(sample(as.factor(sub_int_zeros[,1]), (normalization_value-length(sub_int_bind_ordered_V1)), replace = F)))
                      lowest_level_2
                      lowest_level_3 <- c(sub_int_bind_ordered_V1,lowest_level_2)
                      
                      y <- as.vector(rep(0,length(x[,1]))) #create an empty vector
                      y[lowest_level_3] = 1 #set the randomly selected counts to 1
                      y
                      
                    } #if integer ranks are < normalization_value, do ranked subsampling of the zero integers
                    else {
                      y <- as.vector(rep(0,length(x[,1]))) #create an empty vector
                      y[sub_int_bind_ordered_V1] = 1 #set the randomly selected counts to 1
                      y
                    } # if integer ranks are = normalization_value, sample all
                  }
                  
                  
                }
                
              }    
              
              SRS <- revtrunc_fixed_factor_1[,i] + ceiling(x[,i] > unique(x[,i][max])) + y #sum it all up
              SRS
              sum(SRS) #verification
              assign(paste(names(data)[i],sep=""),SRS)
              fixed_factor_1 <- cbind(fixed_factor_1, SRS)
              colnames(fixed_factor_1)[{
                counter = counter + 1
              }] <- names(data)[i]
            } #if the sum of the counts in the library > Cmin
          } #for all other libaries 
        }
        
        SRS_output <- fixed_factor_1
        SRS_output
      }
    }
  }
}