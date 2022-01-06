### 
# KATZ'S BACKOFF MODEL
###
#' Values known:
#' n-1 sequence = the first n-1 words of the ngram
#'  
#' Values needed:
#' k = a user set constant determining the threshold below which to treat frequecy 
#'     counts as zero
#'     
#' d = discount ratio of Good-Turing estimate to the raw frequency count (C*/C)
#' 
#' delta = sum of discounted maximum likelihood (n-1)gram probabilities for the
#'         nth word where the nth word is unseen in the ngram, but seen in the
#'         (n-1)gram
#' 
#' beta = the left over probability mass for the (n-1)gram. Said another way,
#'        One minus the sum of all discounted maximum likelihood probabilites of 
#'        the nth word for the ngram that have been seen.
#' 
#' alpha = the ratio of left over probabilty (beta) to the sum of all discounted 
#'         maximum likelihood (n-1)gram probabilites for unseen ngram words (delta)
#'         e.g. beta/delta
#' 


##### Packages #####
library(data.table)
library(ngram)
library()


 
##### Functions #####


#' check.existence() checks for the existence of the n-1gram in the given ngram table
#' and returns TRUE/FALSE

#' NOTES: Here, and throughout the code, 'context words' is used to refer to an ngram without
#'        its last word because the last word is the word being predicted.

#'       Additionally, 'context' is the variable name in the ngramtable containing 
#'       the context words for each ngram.
check.existence <- function(contextwords, ngramtable){
  ngramtable <- as.data.table(ngramtable)
  if (nrow(ngramtable[contextwords == context]) > 0){
    return(TRUE)
  } else {
    return(FALSE)
  }
}


#' rm.lastword() removes the last word from a given ngram and returns the ngram's 
#' context words.
rm.lastword <- function(ngram){
  index.spaces <- gregexpr(' ', ngram)[[1]]
  index.lastspace <- length(index.spaces)
  contextwords <- substr(ngram, 
                     start = 1 , 
                     stop = index.spaces[index.lastspace-1])
  return(contextwords)
}


#' get.lastword() returns the last word of an ngram.
get.lastword <- function(ngram){
  index.spaces <- gregexpr(' ', ngram)[[1]]
  index.lastspace <- length(index.spaces)
  lastword <- substr(ngram, 
                     start = index.spaces[index.lastspace-1] + 1, 
                     stop = index.spaces[index.lastspace] - 1)
  return(lastword)
}


#' get.context() returns the last n-1 words from a given text string (txt), where 
#' n is the second argument.
get.context <- function(txt, n){
  index.spaces <- gregexpr(' ', txt)[[1]]
  txt.len <- nchar(txt)
  if (substr(txt, start = txt.len, stop = txt.len) == ' '){
    n.spaces <- n
    contextwords <- substr(txt, 
                           start = index.spaces[length(index.spaces) - n.spaces] + 1,
                           stop = txt.len - 1)
  } else {
    n.spaces <- n-1
    contextwords <- substr(txt, 
                           start = index.spaces[length(index.spaces) - n.spaces] + 1,
                            stop = txt.len)
  }
  return(contextwords)
}


#' get.discount() calculates the Katz's Back-Off Models' value 'd' by dividing 
#' Good-Turing smoothed frequency counts by the observed frequency counts.
#' 
#' NOTE: k is an integer for which all frequencies <= k are treated as zero
get.discount <- function(ngramtable, k){
  allNs <- ngramtable[, .N, by = freq]
  zero <- data.table('freq' = 0, 'N' = 0)
  allNs <- rbind(allNs, zero)
  
  totalN <- allNs[,sum(freq*N)]
  maxN <- max(allNs$freq) 
  N1 <- allNs[freq == 1, N]
  
  allNs[, q := shift(freq, type = 'lead')]
  allNs[, t := shift(freq, type = 'lag')]
  
  allNs[freq == 0, N := N1]
  allNs[freq == 1, Z := N/(.5*(t-0))]
  allNs[freq > 1 & freq < maxN, Z := N/(.5*(t-q))]
  allNs[freq == maxN, Z := N/(.5*(((2*freq - q) - q)))]
  
  allNs <- allNs[order(freq)]
  
  fit <- lm(log(Z) ~ log(freq), data = allNs[2:nrow(allNs)])
  
  allNs[freq <= k & freq != 0, gtcounts := round(exp(predict(fit, newdata = data.table('freq' = 1:k))),5)
      ][freq > k | freq == 0, gtcounts := N
      ][, gtprob := gtcounts/totalN
      ][, d := gtcounts/N]
  
  ngt.discounted <- ngramtable[allNs[, .(freq,N, gtcounts, d)], on = .(freq), nomatch = 0]
  ngt.discounted[, d := gtcounts/N]
  return(ngt.discounted)
}













