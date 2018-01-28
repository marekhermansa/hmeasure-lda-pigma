# Define a function that calculates the hmeasure
my.HMeasure <- function(true.class)
{
# input check

  true.class <- as.factor(as.character(true.class))
  input.true.class <- levels(true.class)
  
  levels(true.class) <- c('0', '1')
  true.class <- as.factor(true.class)
  true.class <- as.numeric(true.class)-1 # turn into numeric array of 0s and 1s
  
    
  complete.rows <- complete.cases(scores)
  scores <- subset(scores,subset=complete.rows)
  
#===
  
  name.now <- colnames(scores)[1]
  s <- scores[,1]
  n <- length(s)

  y <- subset(true.class,subset=complete.rows)
      
  n1 <- sum(y) 
  n0 <- n-n1
  pi0 <- n0/n
  pi1 <- n1/n
    
  # Set severity ratio to default
  severity.ratio <- pi1/pi0

# Calculate ROC curve
# ROC starts at F0[1]=0, F1[1]=0, and ends at F0[K1]=1, F1[K1]=1.

  # Count the instances of each unique score, and ranks them by score
  s1 <- unname(tapply(y, s, sum))/n1
  s1 <- c(0,s1,1-sum(s1)) # add the points 0,0 and 1,1
  s0 <- unname(tapply(1-y, s, sum))/n0
  s0 <- c(0,s0,1-sum(s0)) # add the points 0,0 and 1,1
    
  # Number of unique scores
  S <- length(s1)
  
  F1 <- cumsum(s1)
  F0 <- cumsum(s0)
  
  
  # Consider ROC above diagonal only
  chull.points <- chull(1-F0,pmax(1-F1,1-F0))

  hc <- length(chull.points)
  
  # Extract shape
  if (severity.ratio > 0){
    shape1 <- 2
    shape2 <- 1+(shape1-1)*1/severity.ratio
  }
  if (severity.ratio < 0){
    shape1 <- pi1+1
    shape2 <- pi0+1
  }

# Calculate h-measure
  
  b00 <- beta(shape1,shape2)
  b10 <- beta(1+shape1,shape2)
  b01 <- beta(shape1,1+shape2)
  
  # lshape
  b0 <- c(1:hc+1)
  b1 <- c(1:hc+1)
  cost <- c(1:(hc+1))
  
  cost[1] <- 0
  cost[hc+1] <- 1
  
  b0[1] <- pbeta(cost[1], shape1=(1+shape1), shape2=shape2)*b10/b00
  b0[hc+1] <- pbeta(cost[hc+1], shape1=(1+shape1), shape2=shape2)*b10/b00
  b1[1] <- pbeta(cost[1], shape1=shape1, shape2=(1+shape2))*b01/b00
  b1[hc+1] <- pbeta(cost[hc+1], shape1=shape1, shape2=(1+shape2))*b01/b00
    
  # lshape
  G0 <- 1-F0[chull.points]
  G1 <- 1-F1[chull.points] 
  
  for (i in 2:hc){
    cost[i] <- pi1*(G1[i]-G1[i-1]) / (pi0*(G0[i]-G0[i-1]) + pi1*(G1[i]-G1[i-1]))
      
    b0[i] <- pbeta(cost[i], shape1=(1+shape1), shape2=shape2)*b10/b00
      
    b1[i] <- pbeta(cost[i], shape1=shape1, shape2=(1+shape2))*b01/b00
  }
    
  LHshape1 <- 0
  for (i in 1:hc){
  LHshape1 <- LHshape1 + pi0*(1-G0[i])*(b0[(i+1)]-b0[i]) + pi1*G1[i]*(b1[(i+1)]-b1[i])
  }
    
  B0 <- pbeta(pi1, shape1=(1+shape1), shape2=shape2)*b10/b00
    
  B1 <- pbeta(1, shape1=shape1, shape2=(1+shape2))*b01/b00 -
    pbeta(pi1, shape1=shape1, shape2=(1+shape2))*b01/b00
    
  H <- 1 - LHshape1/(pi0*B0 + pi1*B1)
    
  return(H)
} 

#===============================================

library(MASS)
library(class)
data(Pima.te)

# Split data set into training and test subsets
n <- dim(Pima.te)[1]
ntrain <- floor(2*n/3)
ntest <- n-ntrain
pima.train <- Pima.te[seq(1,n,3),]
pima.test <- Pima.te[-seq(1,n,3),]
true.class<-pima.test[,8]

# Train the learner
pima.lda <- lda(formula=type~., data=pima.train)
out.lda <- predict(pima.lda,newdata=pima.test)

# Make predictions
class.lda <- out.lda$class
scores.lda <- out.lda$posterior[,2]

# Evaluate the implemented learner
scores <- data.frame(LDA=scores.lda)
my.results <- my.HMeasure(true.class)
my.results

# Evaluate the learner from hmeasure package
scores <- data.frame(LDA=scores.lda)
results <- HMeasure(true.class,scores)
results$metrics[c('H')]