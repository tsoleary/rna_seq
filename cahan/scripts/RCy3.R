# creating a cytoscape for the differentially expressed genes ------------------
# this information from the vignettes for each GENIE3 and RCy3
# https://bioconductor.org/packages/release/bioc/vignettes/GENIE3/inst/doc/GENIE3.html
# https://bioconductor.org/packages/release/bioc/vignettes/RCy3/inst/doc/Overview-of-RCy3.html

require(GENIE3)
require(RCy3)
require(RColorBrewer)

# from the GENIE3 vingnette

exprMatr <- matrix(sample(1:10, 100, replace = TRUE), nrow = 20)
rownames(exprMatr) <- paste("Gene", 1:20, sep = "")
colnames(exprMatr) <- paste("Sample", 1:5, sep = "")
head(exprMatr)

set.seed(123) # For reproducibility of results

# now call the GENIE3() function to get the weighted matrix
weightMat <- GENIE3(exprMatr)


# so now we need to define a subset of genes as canidate regulators
# Genes that are used as candidate regulators
regulators <- c(2, 4, 7)
# Or alternatively:
regulators <- c("Gene2", "Gene4", "Gene7")
weightMat <- GENIE3(exprMatr, regulators = regulators)

weightMat

# genie3 is based on regression trees
# they can be either Random Forests or Extra-Trees methods 

# Use Extra-Trees (ET) method
# 7 randomly chosen candidate regulators at each node of a tree
# 5 trees per ensemble
weightMat <- GENIE3(exprMatr, treeMethod = "ET", K = 7, nTrees = 50)

linkList <- getLinkList(weightMat)
dim(linkList)

# seth said to set a threshold for the linkList at 0.1
linkList <- getLinkList(weightMat, threshold = 0.1)

# RCy3 -------------------------------------------------------------------------

# first you need to lauch Cytoscape the desktop and then you can run these two
# lines
cytoscapePing()
cytoscapeVersionInfo()


# my first network
nodes <- data.frame(id = c("node 0", "node 1", "node 2", "node 3"),
                    group = c("A", "A", "B", "B"), # categorical strings
                    score = as.integer(c(20, 10, 15, 5)), # integers
                    stringsAsFactors = FALSE)
edges <- data.frame(source = c("node 0", "node 0", "node 0", "node 2"),
                    target = c("node 1", "node 2", "node 3", "node 3"),
                    interaction = c("inhibits", "interacts", "activates", 
                                  "interacts"),  # optional
                    weight = c(5.1, 3.0, 5.2, 9.9), # numeric
                    stringsAsFactors = FALSE)

createNetworkFromDataFrames(nodes, edges, 
                            title = "my first network", 
                            collection = "DataFrame Example")

# there are other visual styles that you can try and you can also make a new 
# visual style
setVisualStyle("Marquee")
# to get a list of the different styles available
RCy3::getVisualStyleNames()
setVisualStyle("Minimal")

# imporing data ----------------------------------------------------------------
# now this is from the RCy3 vignette on how to import the data
# https://bioconductor.org/packages/release/bioc/vignettes/RCy3/inst/doc/Importing-data.html

csv <- system.file("extdata", "galExpData.csv", package = "RCy3")
data <- read.csv(csv, stringsAsFactors = FALSE)
loadTableData(data, data.key.column = "name") # for some reason this throws an error
# "Failed to load data: Provided key columns do not contain any matches"


# network functions and visualization ------------------------------------------
lesmis <- system.file("extdata","lesmis.txt", package="RCy3")
dataSet <- read.table(lesmis, header = FALSE, sep = "\t")

gD <- igraph::simplify(igraph::graph.data.frame(dataSet, directed = FALSE))

# verify the number of nodes and edges
igraph::vcount(gD)
igraph::ecount(gD)

# common igraph functions ----
# calculate the degree for all nodes
degAll <- igraph::degree(gD, v = igraph::V(gD), mode = "all")

# calculate the betweenness for all nodes
betAll <- igraph::betweenness(gD, v = igraph::V(gD), directed = FALSE) / (((igraph::vcount(gD) - 1) * (igraph::vcount(gD)-2)) / 2)
betAll.norm <- (betAll - min(betAll))/(max(betAll) - min(betAll))
rm(betAll)

# calculate dice similarities between all nodes
dsAll <- igraph::similarity.dice(gD, vids = igraph::V(gD), mode = "all")

# add new attributes to the network based on these calculated values

gD <- igraph::set.vertex.attribute(gD, "degree", index = igraph::V(gD), value = degAll)
gD <- igraph::set.vertex.attribute(gD, "betweenness", 
                                   index = igraph::V(gD), value = betAll.norm)

summary(gD)

# and now for the edge attributes

F1 <- function(x) {data.frame(V4 = dsAll[which(igraph::V(gD)$name == as.character(x$V1)), which(igraph::V(gD)$name == as.character(x$V2))])}
dataSet.ext <- plyr::ddply(dataSet, .variables=c("V1", "V2", "V3"), function(x) data.frame(F1(x)))

gD <- igraph::set.edge.attribute(gD, "weight", index = igraph::E(gD), value = 0)
gD <- igraph::set.edge.attribute(gD, "similarity", index = igraph::E(gD), value = 0)

for (i in 1:nrow(dataSet.ext)){
  igraph::E(gD)[as.character(dataSet.ext$V1) %--% 
                  as.character(dataSet.ext$V2)]$weight <- as.numeric(dataSet.ext$V3)
  igraph::E(gD)[as.character(dataSet.ext$V1) %--% 
                  as.character(dataSet.ext$V2)]$similarity <- as.numeric(dataSet.ext$V4)
}
rm(dataSet,dsAll, i, F1)
# you should see weight and similarity added to the summary
summary(gD)

# lets put this shit back into cytoscape
createNetworkFromIgraph(gD, title = "Les Miserables", collection = "Books")
# you can choose different layouts and get a list from this function
getLayoutNames()

getLayoutPropertyNames("fruchterman-rheingold") 

layoutNetwork('force-directed defaultSpringLength=70 defaultSpringCoefficient=0.000003')


setNodeColorMapping('degree', 
                    c(min(degAll), 
                      mean(degAll), 
                      max(degAll)), 
                    c('#F5EDDD', '#F59777', '#F55333'))
lockNodeDimensions(TRUE)
setNodeSizeMapping('betweenness', 
                   c(min(betAll.norm),
                     mean(betAll.norm), 
                     max(betAll.norm)), 
                   c(30, 60, 100))




setEdgeLineWidthMapping('weight', 
                        c(min(as.numeric(dataSet.ext$V3)), 
                          mean(as.numeric(dataSet.ext$V3)), 
                          max(as.numeric(dataSet.ext$V3))), 
                        c(1,3,5))
setEdgeColorMapping('weight', 
                    c(min(as.numeric(dataSet.ext$V3)), 
                      mean(as.numeric(dataSet.ext$V3)), 
                      max(as.numeric(dataSet.ext$V3))), 
                    c('#BBEE00', '#77AA00', '#558800'))



setBackgroundColorDefault('#D3D3D3')
setNodeBorderColorDefault('#000000')
setNodeBorderWidthDefault(3)
setNodeShapeDefault('ellipse')
setNodeFontSizeDefault(20)
setNodeLabelColorDefault('#000000')


# okay cool now let's try to incorporate gene networks from string

# https://bioconductor.org/packages/release/bioc/vignettes/RCy3/inst/doc/Cancer-networks-and-data.html

# example from the vignette
string.cmd = 'string disease query disease="breast cancer" cutoff=0.9 species="Homo sapiens" limit=150'
commandsRun(string.cmd)

# get another disease
string.cmd = 'string disease query disease="ovarian cancer" cutoff=0.9 species="Homo sapiens" limit=150'
commandsRun(string.cmd)


# so now let's see what we can do with Cytoscape when we are interacting with R
# to see what networks are out there
getNetworkList()

layoutNetwork(layout.name='circular') 

getLayoutNames()

getLayoutPropertyNames(layout.name='force-directed')
layoutNetwork('force-directed defaultSpringCoefficient=0.0000008 defaultSpringLength=70')

getTableColumnNames('node')

disease.score.table <- getTableColumns('node','stringdb::disease score')
disease.score.table

par(mar=c(1,1,1,1))
plot(factor(row.names(disease.score.table)),disease.score.table[,1], ylab=colnames(disease.score.table)[1])
summary(disease.score.table)

# from the top quartile of the diease score
top.quart <- quantile(disease.score.table[,1], 0.75)
top.nodes <- row.names(disease.score.table)[which(disease.score.table[,1]>top.quart)]
createSubnetwork(top.nodes,subnetwork.name ='top disease quartile')
#returns a Cytoscape network SUID


# subset even more to create a new subnetwork with connected nodes only
createSubnetwork(edges='all',subnetwork.name='top disease quartile connected') 

# from using the main network ----

# set the current network back to the main network
setCurrentNetwork(network="String Network - ovarian cancer")
top.nodes <- row.names(disease.score.table)[tail(order(disease.score.table[,1]),3)]
selectNodes(nodes=top.nodes)
selectFirstNeighbors()
createSubnetwork('selected', subnetwork.name='top disease neighbors') # selected nodes, all connecting edges (default)



# trying to query the drosophila database
string.cmd = 'string disease query disease="breast cancer" cutoff=0.9 species="Homo sapiens" limit=150'
commandsRun(string.cmd)


# visualizing data on networks

load(system.file("extdata","tutorial-ovc-expr-mean-dataset.robj", package="RCy3"))
load(system.file("extdata","tutorial-ovc-mut-dataset.robj", package="RCy3"))
load(system.file("extdata","tutorial-brc-expr-mean-dataset.robj", package="RCy3"))
load(system.file("extdata","tutorial-brc-mut-dataset.robj", package="RCy3"))

str(brc.expr)  # gene names in row.names of data.frame
str(brc.mut)  # gene names in column named 'Hugo_Symbol'

setCurrentNetwork(network="String Network - breast cancer")
layoutNetwork('force-directed') #uses same settings as previously set

?loadTableData
loadTableData(brc.expr,table.key.column = "display name")  #default data.frame key is row.names
loadTableData(brc.mut,'Hugo_Symbol',table.key.column = "display name")  #specify column name if not default


style.name = "dataStyle"
createVisualStyle(style.name)
setVisualStyle(style.name)

setNodeShapeDefault("ellipse", style.name) #remember to specify your style.name!
setNodeSizeDefault(60, style.name)
setNodeColorDefault("#AAAAAA", style.name)
setEdgeLineWidthDefault(2, style.name)
setNodeLabelMapping('display name', style.name)

brc.expr.network = getTableColumns('node','expr.mean')  
min.brc.expr = min(brc.expr.network[,1],na.rm=TRUE)
max.brc.expr = max(brc.expr.network[,1],na.rm=TRUE)
data.values = c(min.brc.expr,0,max.brc.expr)


display.brewer.all(length(data.values), colorblindFriendly=TRUE, type="div") # div,qual,seq,all
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping('expr.mean', data.values, node.colors, style.name=style.name)



# trying to query the drosophila database
string.cmd = 'string disease query disease="breast cancer" cutoff=0.9 species="Homo sapiens" limit=150'
commandsRun(string.cmd)



################################################################################
# how are my top genes related
# https://bioconductor.github.io/BiocWorkshops/cytoscape-automation-in-r-using-rcy3.html#use-case-1---how-are-my-top-genes-related

require(RCy3)
cytoscapePing()

RNASeq_gene_scores <- read.table( 
  file.path(getwd(),"230_Isserlin_RCy3_intro","data","TCGA_OV_RNAseq_All_edgeR_scores.txt"),
  header = TRUE, sep = "\t", quote="\"", stringsAsFactors = FALSE)

# heat and cold shock data -----------------------------------------------------

cold_deg <- res_cold %>%
  dplyr::filter(padj < 0.01)


commandsHelp("help string")
commandsHelp("help string protein query")

cold_deg_string_interaction_cmd <- paste('string gene query taxonID=7227 limit=150 cutoff=0.9 query="', paste(cold_deg$gene, collapse = ","), '"', sep = "")
commandsGET(cold_deg_string_interaction_cmd)

mesen_string_interaction_cmd <- paste('string protein query taxonID=7227 limit=150 cutoff=0.9 query="',paste(cold_deg$gene, collapse=","),'"',sep="")
commandsGET(mesen_string_interaction_cmd)
