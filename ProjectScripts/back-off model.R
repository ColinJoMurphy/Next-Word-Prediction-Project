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



 
##### Helper Functions #####


#' check.existence() checks for the existence of the n-1gram in the given ngram table
#' and returns TRUE/FALSE.

#' NOTES: Here, and throughout the code, 'context words' is used to refer to an ngram without
#'        its last word because the last word is the word being predicted.

#'       Additionally, 'context' is the variable name in the ngramtable containing 
#'       the context words for each ngram.
#'       
#'contextwords: the first n-1 words of an ngram
#'ngramtable: a table containing ngrams from a given vocabulary and information 
#'            about each, i.e., frequency, frequency of frequency, etc.

check.existence <- function(contextwords, ngramtable){
  ngramtable <- as.data.table(ngramtable)
  return(nrow(ngramtable[contextwords == context]) > 0)
}


#' rm.lastword() removes the last word from a given ngram and returns the ngram's 
#' context words.
#' 
#' ngram: a set of n words separated by spaces

rm.lastword <- function(ngram){
  index.spaces <- gregexpr(' ', ngram)
  index.lastspace <- lengths(index.spaces)
  txt.len <- nchar(ngram)
  print(sum(index.lastspace!=2))
  
  index.secondtolast <- index.lastspace - 1 
  # index.secondtolast <- fifelse(index.lastspace == 1, 
  #                              {
  #                               warning("There is an ngram here with only one space")
  #                               1
  #                              }, 
  #                              index.lastspace - 1)
  
  contextwords <- fifelse(substr(ngram, txt.len, txt.len) == ' ',
                          substr(ngram, 
                                 start = 1 , 
                                 stop = mapply('[', index.spaces, index.secondtolast) - 1),
                          substr(ngram,
                                 start = 1,
                                 stop = mapply('[', index.spaces, index.lastspace) - 1)
  )
         
  
  return(contextwords)
}


#' get.lastword() returns the last word of an ngram.
#' 
#' ngram: a set of n words separated by spaces

get.lastword <- function(ngram){
  index.spaces <- gregexpr(' ', ngram)
  index.lastspace <- lengths(index.spaces)
  txt.len <- nchar(ngram)
  
  index.secondtolast <- index.lastspace - 1 
  lastword <- fifelse(substr(ngram, txt.len, txt.len) == ' ',
                      substr(ngram, 
                             start = mapply('[', index.spaces, index.secondtolast) + 1, 
                             stop = mapply('[', index.spaces, index.lastspace) - 1),
             
                      substr(ngram, 
                             start = mapply('[', index.spaces, index.lastspace) + 1, 
                             stop = txt.len)
  )
  return(lastword)
}


#' get.context() returns the last n-1 words from a given text string.
#' 
#' txt: text to use for next word prediction 
#' n: integer refering what size ngram should be used (usually 1,2, or 3)

get.context <- function(txt, n){
  index.spaces <- gregexpr(' ', txt)[[1]]
  txt.len <- nchar(txt)
  if (substr(txt, start = txt.len, stop = txt.len) == ' '){
    n.spaces <- n-1
    contextwords <- substr(txt, 
                           start = index.spaces[length(index.spaces) - n.spaces] + 1,
                           stop = txt.len - 1)
  } else {
    n.spaces <- n-2
    contextwords <- substr(txt, 
                           start = index.spaces[length(index.spaces) - n.spaces] + 1,
                            stop = txt.len)
  }
  return(contextwords)
}


#' get.discount() calculates the Katz's Back-Off Models' value 'd' by dividing 
#' Good-Turing smoothed frequency counts by the observed frequency counts.
#' 
#' ngramtable: a table containing ngrams from a given vocabulary and information 
#'             about each, i.e., frequency, frequency of frequency, etc.
#'k: an integer such that frequency counts <= k are discounted using a Good Turing 
#'   smoothing

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
  
  ngt.discounted <- ngramtable[allNs[, .(freq, N, gtcounts, d)], on = .(freq), nomatch = 0]
  ngt.discounted[, d := gtcounts/N]
  return(ngt.discounted)
}





##### Main Model Function #####


#' predict.nextword() predicts the next possible words in a given phrase using Katz's
#' Back-Off Model
#' 
#' txt: the phrase for which the next word needs predicting
#' k: an integer such that frequency counts <= k are discounted using a Good Turing 
#'   smoothing
#' nsize: an integer indicating the larges ngram to use in the model
#' ngramtablelist: a list containing pre-made ngramtables such that 
#'                 ngramtablelist[[1]] is the unigram table, ngramtablelist[[2]]
#'                 is the bigram table, etc.
#' predictnum: an integer indicating how many next word candidates should be returned

predict.nextword <- function(txt,
                             k,
                             nsize,
                             ngramtablelist,
                             predictnum){
  ### Check if the context words are in the ngramtable for the largest n
  n <- nsize
  contextwords <- get.context(txt = txt, n = n)
  table <- ngramtablelist[[n]]
  table[, context := rm.lastword(ngrams)]
  
  print(table)
  
  exist <- check.existence(contextwords = contextwords, ngramtablelist[[n]])
  if(isTRUE(exist)){
    table <- ngramtablelist[[n]]
    table.sd <- get.discount(table, k = k)
    table.sd[, probability :=  d*(nrow(.SD[context == contextwords])/sum(freq))]
    nextwords <- table.sd[context == contextwords][order(probability)]
    
    print(nextwords)
    if (nrow(nextwords) >= predictnum) { 
            nextwords[1:predictnum, 'nextword' := get.lastword(ngrams)]
      return(nextwords[1:predictnum, .(nextword, probability)])
    } else {  
            nextwords[,'nextword' := get.lastword(ngrams)]
      message("'predictnum' argument is larger than number of candidates from vocabulary\n\n returning all candidates")
      return(nextwords[, .(nextword, probability)])
    }
    
  }#end of if loop for first check true
  
  
  
}#end of main function

nw <- predict.nextword('school', 5, 2, tablelist, predictnum = 5)
nw





























