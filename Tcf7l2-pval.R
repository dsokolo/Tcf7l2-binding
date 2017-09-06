#setwd("path-to-tcf2la-binding.RData")

#geneList <- sample(x = counts,size = 1000)
Tcf2la.Pval <- function(geneList, represent = "over", permutations = 3000, seed = 0) {
# this function takes a gene lists and calculates if the Tcf7l2 binding is over-underrepresented
# by computing a permutation p-value on the total number of binding sites in your gene lists vs the total number of 
# binding sites in an unrelated gene list. The null hypothesis is that there is no difference in binding between
# your list and the unrelated list. each permutation compares your binding to a random list that is differentially selected each time
# if your list is over or under represented, then you have a sucess (depending on what you set)
# p = (# permutations) - (# sucesses) / (# permutations). the more permutations the more accurate and
# the longer the running time.
#inputs
    # geneList = a list of mouse gene symbols, MGI/Ensembl/Entrez no accepted
    # representation = 'over' or 'under' depending on what you want
    # permutaitons = how many permutations you want to run this for
    # seed = what seed we're selecting our random populations from

load("tcf2la-binding.RData") #load named vector of gene lists with Chip-seq binding scores
bindingVector <- thecount.ol[names(thecount.ol) %in% geneList] # get vector of binding for our gene
#hist(bindingVector, main = "distribution of Tcf7l2 binding in list")
emperical.score <- sum(bindingVector) # get total binding sites in our gene list

success <- 0
for(i in 1:permutations){ 
  random.score <- sum(sample(thecount.ol, length(geneList))) #get total binding sites of random list of genes within our population
  if(represent == "over") { 
    if(emperical.score > random.score) {
      success <- success + 1
    }
  }
  if(represent == "under") {
    if(emperical.score < random.score) {
      success <- success + 1
    }
  }
  if(!(represent %in% c("over","under"))) {
    return("Error: You did not specify representation (or you made a spelling mistake), it is 'over' or 'under'")
  }
  
    
}
p.value <- (permutations - success) / permutations

ans <- list(pVal = p.value, Representation = paste0("this was a test for ",represent,"-representation of Tcf7l2 binding"),
            bindingScore = emperical.score)
return(ans)
}
#The code below here tests a random selection of genes against itself before
#calculating the p-value. The proccess is repreated 2000 times (1000 overrepresented, 1000 underrepresented)
#and then qqplots are generated to insure that the p-values are normally distributed.
#currently this code is commented out because it takes about 30 minutes to run


#pValsOver <- c()
#pValsUnder <- c()
#for(i in 1:1000) {
#  test <- Tcf2la.Pval(geneList)
#  testUnder <- Tcf2la.Pval(geneList, represent = "under")
#  pValsOver[i] <- test$pVal
#  pValsUnder[i] <- testUnder$pVal

#  print(i)
#}
#qqnorm(pValsOver, main = "qqplot of over-representation test with random genes selected")
#qqnorm(pValsUnder, main = "qqplot of under-representation test with random genes selected")

