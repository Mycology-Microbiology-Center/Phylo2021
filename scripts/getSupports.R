## The getSupports function compares all the subclades in the primary tree (the
## one used for visualization) to all the subclades in the secondary tree
## (needed only for its supports values), and vice versa, and identifies all
## those clades that have identical content tip-wise. It produces a data frame
## supportsTable with node IDs of primaryTree and their respective supports from
## secondaryTree, for further use in visualizations. Primary and secondary trees
## can be, e.g., Bayesian tree and Maximum Likelihood bootstrap tree, or vice
## versa.

## Requires geiger and ape (which is a dependency for geiger).
require(geiger)

getSupports <- function(primaryTree, secondaryTree)
{
  ## The getAllSubTrees subfunction below atomizes a tree into each individual
  ## subclade, so then we can compare and match these subclades between our 2 trees.
  getAllSubtrees <- function(phy, minSize = 2)
  {
    res <- list()
    count = 1
    ntip <- length(phy$tip.label)
    for (i in 1:phy$Nnode)
    {
      l <- tips(phy, ntip + i)
      bt <- match(phy$tip.label, l)
      if (sum(is.na(bt)) == 0)
      {
        st <- phy
      }
      else
        st <- drop.tip(phy, phy$tip.label[is.na(bt)])
      if (length(st$tip.label) >= minSize)
      {
        res[[count]] <- st
        count <- count + 1
      }
    }
    res
  }
  getAllSubtrees(primaryTree) -> primarySub
  getAllSubtrees(secondaryTree) -> secondarySub
  supportsList <- matrix('-', Nnode(primaryTree), 1)
  for (i in 1:Nnode(primaryTree))
  {
    for (j in 1:Nnode(secondaryTree))
    {
      match(primarySub[[i]]$tip.label[order(primarySub[[i]]$tip.label)], secondarySub[[j]]$tip.label[order(secondarySub[[j]]$tip.label)]) -> shared
      match(secondarySub[[j]]$tip.label[order(secondarySub[[j]]$tip.label)], primarySub[[i]]$tip.label[order(primarySub[[i]]$tip.label)]) -> shared2
      if (sum(is.na(c(shared, shared2))) == 0)
      {
        secondaryTree$node.label[j] -> supportsList[i]
        supportsTable <-
          data.frame(min(as.data.frame(primaryTree[["edge"]])$V1):max(as.data.frame(primaryTree[["edge"]])$V1),
                     supportsList)
        rownames(supportsTable) <- NULL
        colnames(supportsTable) <- c("node", "bootstrap")
      }
    }
  }
  return(supportsTable)
}
