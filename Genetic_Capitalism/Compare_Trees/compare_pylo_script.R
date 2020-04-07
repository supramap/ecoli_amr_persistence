rm(list=ls())
# A script to compare phylogenetic trees
library(ape)
# Read in the newick files here
github_tree <- read.tree(file = 'PDG000000004.1024.reference_target.tree.newick.txt')
ncbi_tree <- read.tree(file = 'PDG000000004.1024.newick.txt')
# For the following, plot=TRUE may be set, but the tree topology is very large
# For unrooted
unroot_compare <- comparePhylo(github_tree, ncbi_tree, plot=FALSE)
# Forces the comparison as rooted
root_compare <- comparePhylo(github_tree, ncbi_tree, plot=FALSE, force.rooted = TRUE)
# Another method, might have difficulties with unrooted topology per https://rdrr.io/cran/ape/man/all.equal.phylo.html
final_equal <- all.equal(github_tree, ncbi_tree)
final_identical <- identical(github_tree, ncbi_tree)
final_equal_list <- all.equal.list(github_tree, ncbi_tree)

# Printing these results demonstrates that both newick files are identical