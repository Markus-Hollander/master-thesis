library("igraph")

args = commandArgs(trailingOnly=TRUE)

df = read.table(args[1])
network = graph.data.frame(df, directed=TRUE)

n = as.integer(args[3])

rewired = rewire(network, with = keeping_degseq(loops=FALSE, niter=n * 2))

write_graph(rewired, args[2], format="ncol")
