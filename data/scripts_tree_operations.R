## OPERATIONS ON PHYLOGENETIC TREES WITH GGTREE
## For documentation on ggtree and treeio see: 
## https://yulab-smu.top/treedata-book/index.html
## Original code for plotting ML and bootstrap supports: 
## https://www.protocols.io/view/plotting-phylogenetic-tree-with-branch-supports-fr-n9fdh3n?step=1

## To install ggtree / treeio
# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# BiocManager::install(c("ggtree", "treeio"))

library(treeio)
library(ggtree)
library(ape)
library(ggplot2)
library(geiger)

## TREE IMPORT

## Read trees: ML tree and Bayesian tree. This works for MrBayes output built
## with conformat=simple. If you used default conformat=FigTree in MrBayes block,
## import like specified the detour below.
bootTree <- treeio::read.newick ('data/tree_raxml.nwk')
bayesTree <- ape::read.nexus ('data/tree_mrbayes.nex')

## This Bayesian tree file includes 2 trees: first with  posterior supports, and
## second without them. We will operate on the first.
bayesTree <- bayesTree[[1]]

######## DETOUR TO TRANSFORM DEFAULT MRBAYES FORMATTING

## For downstream plotting we must have posteriors as node labels. Default
## MrBayes output (conformat=FigTree) contains posteriors in comments that, when
## imported, go to 'data' part of treedata object, and not in 'phylo$node.label'
## where we need them. It can be fixed as follows.

# bayesTree <- read.mrbayes ('data/tree_mrbayes_figtree.nex')

## Sort nodes in data to match node order in phylo. (Below is MrBayes
## version, but with few adjustments it also works for BEAST)
# bayesTree@data <- bayesTree@data[order(as.numeric(bayesTree@data$node)), ]

## Write posteriors as node labels. 
## If tree is already rooted, do not append '' to vector. 
# bayesTree@phylo$node.label <-
#   append ('', as.vector(subset(
#     bayesTree@data$prob,
#     as.numeric(bayesTree@data$node) > length(bayesTree@phylo$tip.label)
#   )))

## Translate treedata into phylo.
# bayesTree <- as.phylo(bayesTree)

######## END OF DETOUR

## Visual check of the trees.
ggtree(bootTree) +
  geom_tiplab() +
  geom_text2(
    aes(label = label, subset = !isTip),
    hjust = 1.3,
    vjust = -0.4,
    size = 2
  ) +
  xlim(0, 0.3)

ggtree(bayesTree) +
  geom_tiplab() +
  geom_text2(
    aes(label = round(as.numeric(label), 2), subset = !isTip),
    hjust = 1.3,
    vjust = -0.4,
    size = 2
  ) +
  xlim(0, 0.3)
## Note: Warnings about NAs introduced by coercion can be safely ignored.

## REROOTING

## For the next step we need to root both trees. To root the tree, you need to
## know either 1) labels of outgroup taxa, or 2) an ID of a node joining all the
## taxa belonging to the outgroup.

## First rerooting option, by specifying outgroup labels.
bayesTree_rooted <-
  ape::root(bayesTree,
            outgroup = c('MM35985', 'MYX328'),
            edgelabel = TRUE)
## note: edgelabel = T is essential for correct label placement, otherwise
## labels can slip to wrong branches. For details see Czech et al. 2017 "A
## Critical Review on the Use of Support Values [...]"
## https://doi.org/10.1093/molbev/msx055

## Same for ML tree
bootTree_rooted <-
  ape::root(bootTree,
            outgroup = c('MM35985', 'MYX328'),
            edgelabel = TRUE)

## If outgroup is large, the option with node ID may be preferable. The easiest
## way to locate needed ID is to plot all of them (on example of bayesian tree).
# ggtree(bootTree) +
#   geom_tiplab() +
#   geom_text(
#     aes(label = node),
#     hjust = 1.3,
#     vjust = 1.3,
#     size = 2,
#     color = 'red'
#   ) +
#   xlim(0, 0.3)

## If your tree is large and crowded, it may be difficult to spot the desired
## node. Another way is to make a table from the tree, find a taxon from the
## outgroup and trace node-to-node connections down to the common node of all
## outgroup taxa. 
# bayesTree_tib <- as_tibble(bayesTree)

## In case of Bayesian tree it's node 89.
# bayesTree_rooted <- ape::root(bayesTree, node = 89, edgelabel = TRUE)

## check the result
ggtree(bootTree_rooted) +
  geom_tiplab() +
  geom_text2(
    aes(label = round(as.numeric(label), 2), subset = !isTip),
    hjust = 1.3,
    vjust = -0.4,
    size = 2
  ) +
  xlim(0, 0.3)

ggtree(bayesTree_rooted) +
  geom_tiplab() +
  geom_text2(
    aes(label = round(as.numeric(label), 2), subset = !isTip),
    hjust = 1.3,
    vjust = -0.4,
    size = 2
  ) +
  xlim(0, 0.3)

## COMBINE BOOTSTRAPS AND POSTERIORS

## The getSupports function below compares all the subclades in the Bayes tree
## to all the subclades in the bootstrap tree, and vice versa, and identifies
## all those clades that are identical. It produces a data frame supportsTable
## with node IDs and their respective bootstraps for further use in ggtree.

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

## Apply function to our trees and proceed to tree annotation.
supportsTable <- getSupports(primaryTree = bayesTree_rooted, secondaryTree = bootTree_rooted)

## RELABELING

## Read substitution table for tip labels. The first column must contain the
## same values as tips in the tree, the other column(s) - your pretty labels.
labels <- read.csv('data/labels_simple.csv')

## Add bootstrap values obtained with a getSupports function above.
tree_data <- ggtree(bayesTree_rooted) %<+% supportsTable
## Note: in the same way, the table can contain any other data you may want to
## plot, e.g. values of morphological traits. %<+% is a ggtree-specific operator, 
## basically left join.

## Combine tree and new labels.
tree_data <- tree_data %<+% labels

## Write the result.
tree_data_lab <-
tree_data +
  geom_tiplab(aes(label = label_pretty), offset = 0.001) +
  geom_label2(
    aes(
      label = paste(
        ifelse(is.na(round(as.numeric(label), 2)),
        '-', round(as.numeric(label), 2)),
        ifelse(is.na(bootstrap),
               '-', bootstrap),
        sep = '/'
      ),
      subset = !isTip & label != 'Root'
    ),
    hjust = 1.1,
    vjust = -0.3,
    alpha = 0.8,
    label.size = NA,
    label.padding = unit(0.2, 'mm'),
    size = 2
  ) +
  xlim(0, 0.5)
## Note: xlim is needed to fix truncation of long labels. You may need to change
## the second parameter or even drop xlim altogether, depending on your tree.

## Check the result.
tree_data_lab

## You can change some default aesthetics of the tree post mortem.
## Change all edges from round to miter.
tree_data_lab[["layers"]][[1]][["geom_params"]][["lineend"]] <- 'square'
tree_data_lab[["layers"]][[2]][["geom_params"]][["lineend"]] <- 'square'
## Change thickness of tree lines.
tree_data_lab[["layers"]][[1]][["geom"]][["default_aes"]][["size"]] <- 0.7

## Check the result.
tree_data_lab

## Save the result.
dir.create('output')
ggsave(
  filename = 'output/tree_data_lab.pdf',
  device = 'pdf',
  width = 210,
  height = 297,
  units = 'mm'
)
## Note: This saves the last plotted graph. If you need to save a particular
## object, specify plot = object_name

## BONUS RENAMING

## Read table with more information for labels.
labels_extended <- read.csv('data/labels_extended.csv', na.strings = '')

## Add more data.
tree_data <- tree_data %<+% labels_extended

## Plot the tree with ridiculously detailed labels.
tree_data +
  geom_tiplab(
    aes(label = paste(
      label,
      ifelse(is.na(Organism),
             '', paste(',', Organism)),
      ifelse(is.na(Country),
             '', paste(',', Country)),
      ifelse(is.na(Altitude),
             '', paste(', altitude:', Altitude, 'm')),
      sep = ''
    )),
    geom = 'label',
    label.size = NA,
    label.padding = unit(0, 'mm'),
    offset = 0.001
  ) +
  geom_label2(
    aes(
      label = paste(
        ifelse(is.na(round(as.numeric(label), 2)),
        '-', round(as.numeric(label), 2)),
        ifelse(is.na(bootstrap),
               '-', bootstrap),
        sep = '/'
      ),
      subset = !isTip & label != 'Root'
    ),
    hjust = 1.1,
    vjust = -0.3,
    alpha = 0.8,
    label.size = NA,
    label.padding = unit(0.2, 'mm'),
    size = 2
  ) +
  xlim(0, 0.5)

## Save the result.
ggsave(
  filename = 'output/tree_data_lab_ext.pdf',
  device = 'pdf',
  width = 420,
  height = 297,
  units = 'mm'
)

## BONUS HIGHLIGHTING CLADES

## For this purpose ggtree provides several geoms, e.g. geom_hilight (yep),
## geom_balance (see again https://yulab-smu.top/treedata-book/chapter5.html).
## With geom_cladelabel we additionally can mark clades with labeled vertical
## bars. We will highlight 2 clades and add 1 labeled bar by adding these
## geometries.

tree_data +
  geom_label2(
    aes(
      label = paste(
        ifelse(is.na(round(as.numeric(label), 2)),
        '-', round(as.numeric(label), 2)),
        ifelse(is.na(bootstrap),
               '-', bootstrap),
        sep = '/'
      ),
      subset = !isTip & label != 'Root'
    ),
    hjust = 1.1,
    vjust = -0.3,
    alpha = 0.8,
    label.size = NA,
    label.padding = unit(0.2, 'mm'),
    size = 2
  ) +
  geom_hilight(
    node = 65,
    fill = "#1B9E77",
    alpha = 0.3,
    extend = 0.22
  ) +
  geom_hilight(
    node = 56,
    fill = "#F54748",
    alpha = 0.3,
    extend = 0.20
  ) +
  geom_cladelabel(
    56,
    "SP. NOV.",
    offset = 0.12,
    barsize = 2,
    align = TRUE,
    angle = 90,
    offset.text = 0.008,
    extend = 0.5,
    hjust = 0.5,
    fontsize = 5
  ) +
  geom_tiplab(aes(label = label_pretty), offset = 0.001) +
  xlim(0, 0.5)

## Save the result.
ggsave(
  filename = 'output/tree_data_highlighted.pdf',
  device = 'pdf',
  width = 210,
  height = 297,
  units = 'mm'
)
