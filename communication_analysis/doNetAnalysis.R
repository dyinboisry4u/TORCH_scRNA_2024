# options(future.rng.onMisue = "ignore")
suppressPackageStartupMessages(library(CellChat, lib.loc = "/share/home/Grape/software/install_pkg/miniconda3/envs/mainenv/lib/R/library"))
future::plan("multisession",workers=3)
options(future.globals.maxSize=80000*1024^2)

# fix future bug, ref
# https://github.com/sqjin/CellChat/pull/703
netAnalysis_computeCentrality <- function(object = NULL, slot.name = "netP", net = NULL, net.name = NULL, thresh = 0.05) {
  if (is.null(net)) {
    prob <- methods::slot(object, slot.name)$prob
    pval <- methods::slot(object, slot.name)$pval
    pval[prob == 0] <- 1
    prob[pval >= thresh] <- 0
    net = prob
  }
  if (is.null(net.name)) {
    net.name <- dimnames(net)[[3]]
  }
  if (length(dim(net)) == 3) {
    nrun <- dim(net)[3]
    p <- progressr::progressor(nrun)
    centr.all <- future.apply::future_sapply(
      X = 1:nrun,
      FUN = function(x) {
        Sys.sleep(1/nrun)
        p(sprintf("%g of %g", x, nrun)) # Use with_progress() to see progress bar in client-side
        net0 <- net[ , , x]
        computeCentralityLocal(net0)
      },
      future.seed = TRUE,
      simplify = FALSE
    )
  } else {
    centr.all <- as.list(computeCentralityLocal(net))
  }
  names(centr.all) <- net.name
  if (is.null(object)) {
    return(centr.all)
  } else {
    slot(object, slot.name)$centr <- centr.all
    return(object)
  }
}

computeCentralityLocal <- function(net) {
  centr <- vector("list")
  G <- igraph::graph_from_adjacency_matrix(net, mode = "directed", weighted = T)
  centr$outdeg_unweighted <- rowSums(net > 0)
  centr$indeg_unweighted <- colSums(net > 0)
  centr$outdeg <- igraph::strength(G, mode="out")
  centr$indeg <- igraph::strength(G, mode="in")
  centr$hub <- igraph::hub_score(G)$vector
  centr$authority <- igraph::authority_score(G)$vector # A node has high authority when it is linked by many other nodes that are linking many other nodes.
  centr$eigen <- igraph::eigen_centrality(G)$vector # A measure of influence in the network that takes into account second-order connections
  centr$page_rank <- igraph::page_rank(G)$vector
  igraph::E(G)$weight <- 1/igraph::E(G)$weight
  centr$betweenness <- igraph::betweenness(G)
  #centr$flowbet <- try(sna::flowbet(net)) # a measure of its role as a gatekeeper for the flow of communication between any two cells; the total maximum flow (aggregated across all pairs of third parties) mediated by v.
  #centr$info <- try(sna::infocent(net)) # actors with higher information centrality are predicted to have greater control over the flow of information within a network; highly information-central individuals tend to have a large number of short paths to many others within the social structure.
  centr$flowbet <- tryCatch({
    sna::flowbet(net)
  }, error = function(e) {
    as.vector(matrix(0, nrow = nrow(net), ncol = 1))
  })
  centr$info <- tryCatch({
    sna::infocent(net, diag = T, rescale = T, cmode = "lower")
    # sna::infocent(net, diag = T, rescale = T, cmode = "weak")
  }, error = function(e) {
    as.vector(matrix(0, nrow = nrow(net), ncol = 1))
  })
  return(centr)
}

args=commandArgs(T)
cellchat.obj=args[1]
input.dir=args[2]
out.dir=args[3]

print(paste0("Get ", input.dir, cellchat.obj, "_cellchat_obj.rds"))
cellchat<-readRDS(paste0(input.dir, cellchat.obj, "_cellchat_obj.rds"))
cellchat<-netAnalysis_computeCentrality(cellchat, slot.name="netP")
saveRDS(cellchat,paste0(out.dir, cellchat.obj, "_withCentrality.rds"))
print(paste0("Finish ", out.dir, cellchat.obj, "_withCentrality.rds"))
