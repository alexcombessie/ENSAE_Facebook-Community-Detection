################################INTRODUCTION#####################################################

#installing missing libraries
install.packages('MCL')
install.packages('ape')
install.packages('stringr')
install.packages("combinat")
install.packages("tables")
install.packages("xtables")
install.packages("partitions")

#Loading libraries
library(igraph)
library(MASS)
library(MCL)
library(ape)
library(stringr)
library(combinat)
library(tables)
library(xtables)
library(partitions)
library(moments)
library(reldist)
library(stargazer)
library(HistogramTools)
library(RCurl)

#!!!!!!!!!!!!CHANGE THIS FOR EACH USER!!!!!!!!
workingdir<-"D:/Users/acombess/Dropbox/ENSAE/Semestre 1/Statistical Analysis of Network Data/Project SAND/4-Code"
#workingdir<-"C:/Users/dduong/OneDrive/ProjetSAND"
#workingdir<-"C:/Users/wymeka/Documents/ENSAE/Network data"

setwd(workingdir)


################################LOADING DATA#####################################################
#See script Data Loading

#Generating the list of egonet files
directory <- "./raw-data/egonets/"
egonets_list<-sapply(list.files(directory), function(x) paste(directory, x, sep="/"), USE.NAMES = FALSE)
tr_directory <- "./raw-data/training/"
training_list<-sapply(list.files(tr_directory), function(x) paste(tr_directory, x, sep="/"), USE.NAMES = FALSE)

#Generating the list of all egonetworks objects
n_egonets<-length(egonets_list)
graph_list<-vector("list", length=n_egonets)
for (i in 1:n_egonets) {
    graph_list[[i]]<-CreateGraph(egonets_list[i])
}

#Saving all adjacency matrix as csv in a dedicated folder
for (x in egonets_list) {
  CreateAdjacency(x,"AdjacencyCSV")
}

#Apply the function to get all egonetwork objects from adjacency CSVs
adjancyFiles = list.files("./AdjacencyCSV/")
facebook_graphs = vector("list", length=length(adjancyFiles))
truecircle_list<-vector("list", length=length(adjancyFiles))
i=1;j=1
for(file in adjancyFiles){
  tmp_fn = paste("./raw-data/training//",regmatches(file,gregexpr("[0-9]+", file)) ,".circles", sep="")
  facebook_graphs[[i]] = FromCSVtoiGraph(str_replace_all(paste("./AdjacencyCSV/",file)," ",""))
  if (tmp_fn %in% training_list){
    trueCircles = read.csv(file = tmp_fn, header = FALSE, stringsAsFactors=FALSE, sep=' ')
    trueCircles = trueCircles[,-1]
    truecircle_list[[i]] <- trueCircles
  }else{
    truecircle_list[[i]] <- NULL
  }
  i=i+1;
}
#true circles
head(truecircle_list)
head(adjancyFiles)
# Load true circles 
#trueCircles = loadCircles("./raw-data/training/")


#One network for testing purposes
cc <- FromCSVtoiGraph("./AdjacencyCSV/Adjacency1839.csv")

################################ADDING FEATURES#####################################################

#Gather all possible attributes for vertices
allAttributes= list()
for(i in seq(1,length(facebook_graphs))){
  allAttributes = union(allAttributes,list.vertex.attributes(facebook_graphs[[i]]))
}

#Initialize all attributes vertices to NA , except for the "id" attribute
for (g in seq(1,length(facebook_graphs))) {
  for(attribute in allAttributes){
    if(attribute != "id"){
      vertex_attr(facebook_graphs[[g]],attribute)=rep(NA,length(V(facebook_graphs[[g]])))
    }
  }
}

features <- read.csv("./features/features.csv", sep=',',stringsAsFactors = FALSE)

userAttributes = apply(features,1,splitLine)

# Transform features to a key value list
keyValueAttributes = apply(as.array(userAttributes),1,parseFeatures)

# Set each vertex attributes
for (g in seq(1,length(facebook_graphs))) {
  for (vertexId in V(facebook_graphs[[g]])$id) {
    userId = as.integer(vertexId)
    k = userId +1
    index = which(as.array(V(facebook_graphs[[g]])$id == vertexId))
    keysvalues = keyValueAttributes[[k]]
    for (j in seq(3,ncol(keysvalues))) {
      oldAttributeValue = get.vertex.attribute(facebook_graphs[[g]],keysvalues[1,j],index =
                                                 index)
      newAttributeValue = keysvalues[2,j]
      if (length(oldAttributeValue)!=0) {
        if (!is.na(oldAttributeValue)) {
          if(keysvalues[1,j]!="id"){
            newAttributeValue = paste(oldAttributeValue,newAttributeValue,sep = ",")
          }
        }
      }
      # !!! Bottleneck !!!
      # each keyvalues first line stands for the attribute name and the second one for its value
      facebook_graphs[[g]] = set.vertex.attribute(facebook_graphs[[g]], keysvalues[1,j],index =
                                                    index,value = newAttributeValue)
    }
  }
  print(g)
}

################################WEIGHTING#####################################################

facebook_graphs_weights<-vector("list", length=n_egonets)

# THE ITERATIONS WERE VERY VERY LONG SO ONLY 1:91 WERE DONE
for(i in 1:n_egonets){
  facebook_graphs_weights[[i]]<-rep(NA,ecount(facebook_graphs[[i]]))
  pb <- txtProgressBar(min = 1, max = ecount(facebook_graphs[[i]]), style = 3)
  for(j in 1:ecount(facebook_graphs[[i]])){
    facebook_graphs_weights[[i]][j]<-similarity(facebook_graphs[[i]],E(facebook_graphs[[i]])[[j]])
    setTxtProgressBar(pb, j)
  }
}

#In this section we consider only networks with the ego removed

#Remove the ego from all the networks
facebook_graphs_noego<-vector("list", length=n_egonets)
for (i in 1:n_egonets)  {
  ego_vertex<-V(facebook_graphs[[i]])[[1]]$name
  facebook_graphs_noego[[i]]<-delete_vertices(facebook_graphs[[i]],ego_vertex)}

##Removing ego edges and adding +1 to all weights (standardization)
facebook_graphs_weights_noego<-vector("list", length=n_egonets)
for(i in 1:91){
  facebook_graphs_weights_noego[[i]]<-1+tail(facebook_graphs_weights[[i]],ecount(facebook_graphs_noego[[i]])-ecount(facebook_graphs[[i]]))
}

#Adding weights as attributes to the graphs
facebook_graphs_noego_withweights<-facebook_graphs_noego
for(i in 1:91){
  E(facebook_graphs_noego_withweights[[i]])$weight<-facebook_graphs_weights_noego[[i]]
}

################################DESCRIPTIVE STATISTICS#####################################################

#Visualizing a few networks, with the ego this time
par(mfrow=c(2,2),mar=c(1,1,1,1))
graph_list_sample4=sample(facebook_graphs,4)
for (x in graph_list_sample4){
  SimpleGraphPlot_nolabel(x)}

#Size of the graph (Number of edges)
sizegraph_list<-rep(NA,n_egonets)
for (i in 1:n_egonets){
  sizegraph_list[i]<-gsize(facebook_graphs_noego[[i]])}
summary(sizegraph_list)

#Order of the graph (Number of vertices)
ordergraph_list<-rep(NA,n_egonets)
for (i in 1:n_egonets){
  ordergraph_list[i]<-gorder(facebook_graphs_noego[[i]])}
summary(ordergraph_list)

#Plot both order and size distribution
par(mfrow=c(1,2),mar=c(5.1,4.1,4.1,2.1))
hist(ordergraph_list,col="grey84", main="", 
     xlab="Network Order",ylab="Number of networks")
hist(sizegraph_list,col="grey42", main="",xlab="Network Size",ylab="Number of networks")

#Output for latex
xtable(cbind(summary(ordergraph_list),summary(sizegraph_list)),digits=0)
skewness(ordergraph_list)
gini(ordergraph_list)
skewness(sizegraph_list)
gini(sizegraph_list)

#Degree distribution
#Strength has no meaning here as the graphs are not weighted (yet)
degree_all<-c()
for (i in 1:n_egonets){
  degree_all<-c(degree_all,degree(facebook_graphs_noego[[i]],mode="all"))}
par(mfrow=c(1,1))
hist(degree_all,xlab="Degree", ylab="Frequency",
     main="",col="grey")
xtables(summary(data.frame(values=degree_all,row.names = NULL)),digits=0)
skewness(degree_all)
gini(degree_all)

#Fit a power law to the degree distribution
par(mfrow=c(1,1),mar=c(4.5,4.1,1,1))
log_degree_powerlaw(graph_list_noego)
#output table regression for latex
stargazer(powerlaw_reg(graph_list_noego),no.space=TRUE,font.size ="footnotesize")

#Distribution of closeness, betweeness and eigenvalue
closeness_all<-c()
betweenness_all<-c()
eigen_centrality_all<-c()
for (i in 1:n_egonets){
  closeness_all<-c(closeness_all,closeness(facebook_graphs_noego[[i]]))
  betweenness_all<-c(betweenness_all,betweenness(facebook_graphs_noego[[i]]))
  eigen_centrality_all<-c(eigen_centrality_all,eigen_centrality(facebook_graphs_noego[[i]])$value)}
par(mfrow=c(3,1))
hist(log(closeness_all),xlab="Log-scale",ylab="Frequency",
     main="Closeness Centrality",col="grey30",breaks=20)
hist(log(betweenness_all),xlab="Log-scale", ylab="Frequency",
     main="Betweenness Centrality",col="grey50",breaks=20)
hist(log(eigen_centrality_all),xlab="Log-scale",ylab="Frequency",
     main="Eigenvector Centrality",col="grey70",breaks=20)
#Latex table
tableginiskewness<-data.frame(
  Mean=c(mean(closeness_all),mean(betweenness_all),mean(eigen_centrality_all)),
  Skewness=c(skewness(closeness_all),skewness(betweenness_all), skewness(eigen_centrality_all)),
  Gini=c(gini(closeness_all),gini(betweenness_all),gini(eigen_centrality_all)),
  row.names = c("Closeness","Betweenness","Eigenvector"))
xtable(tableginiskewness)


#Distribution of the density of the graphs
density_list<-rep(NA,n_egonets)
for (i in 1:n_egonets){
  density_list[i]<-graph.density(facebook_graphs_noego[[i]])}
par(mfrow=c(1,1))
hist(density_list,xlab="",main="",col="grey42",breaks="Scott")
summary(density_list)
xtable(summary(data.frame(values=density_list,row.names = NULL)))

#Component decomposition
component_list<-vector("list", length=n_egonets)
ratio_component_list<-rep(NA,n_egonets)
components_size<-c()
pb <- txtProgressBar(min = 1, max = n_egonets, style = 3)
for (i in 1:n_egonets){
  component_list[[i]]<-sapply(decompose.graph(facebook_graphs_noego[[i]]), vcount)
  components_size<-c(components_size,component_list[[i]])
  ratio_component_list[i]<-component_list[[i]][1]/sum(component_list[[i]])
  setTxtProgressBar(pb, i)}
par(mfrow=c(1,1))
hist(ratio_component_list,xlab="",main="",col="grey70",breaks="Scott")
par(mfrow=c(1,1))
plot(ecdf(components_size),main="",xlab="Order of the component",
     ylab="Empirical cumulative distribution",col="red",ylim=c(0.65,1.05))


#Clique decomposition
order_clique<-rep(NA,n_egonets)
pb <- txtProgressBar(min = 1, max = n_egonets, style = 3)
for (i in 1:n_egonets){
  order_clique[i]<-clique_num(facebook_graphs_noego[[i]])
  setTxtProgressBar(pb, i)}
par(mfrow=c(1,1),mar=c(2.5,5,1,1))
hist(order_clique,xlab="",main="",col="grey60",breaks="Scott")
summary(order_clique)


#Feature completion
length(list.vertex.attributes(facebook_graphs[[12]]))
completion_bynetwork<-rep(NA,n_egonets)
feature_notna<-vector("list", length=n_egonets)
number_vertexfeatures<-vector("list", length=n_egonets)
pb <- txtProgressBar(min = 1, max = n_egonets, style = 3)
for (i in 1:n_egonets){
  feature_notna[[i]]<-lapply(vertex_attr(facebook_graphs[[i]]),function(x) sum(!is.na(x)))
  completion_bynetwork[i]<-sum(!is.na(unlist(vertex_attr(facebook_graphs[[i]]))))/length(unlist(vertex_attr(facebook_graphs[[i]])))
  setTxtProgressBar(pb, i)}
feature_notna_complete<-tail(unlist(lapply(do.call(Map, c(c, feature_notna)),
                                           function(x) sum(x)/(sum(ordergraph_list)+n_egonets))),-1)
feature_notna_complete_simple<-feature_notna_complete[names(feature_notna_complete)[!grepl("id",names(feature_notna_complete))]]
par(mfrow=c(1,1),mar=c(2.5,12,1,2))
barplot(sort(feature_notna_complete_simple,decreasing = F),
        main="",horiz=T,las=1,cex.names=0.85,xaxt="n")
axis(1,at=c(0,0.2,0.4,0.6,0.8,1),labels=paste(pretty(c(0,20,40,60,80,100)), " %"))
grid(ny=NA)
summary(completion_bynetwork)
par(mfrow=c(1,1),mar=c(2.5,5,1,1))
hist(completion_bynetwork,xlab="",main="",col="grey88",breaks="Scott",xaxt="n")
axis(1,at=c(0.2,0.25,0.3,0.35,0.4),labels=paste(pretty(c(20,25,30,35,40)), " %"))


#Social circles analysis
par(mfrow=c(1,1),mar=c(5,5,1,1))
hist(unlist(circlecount("./raw-data/training/")),main="",xlab="Number of circles",col="grey30")
hist(circle_order("./raw-data/training/"),main="",xlab="Number of users per circle",col="grey50")


barplot(sapply(overlap("./raw-data/training/")[1:7],function(x) x/sum(overlap("./raw-data/training/")[1:7])),names.arg=1:7,
        xlab="Number of circles per user",ylab="% of total users",col="grey70",yaxt="n")
axis(2,at=c(0,0.1,0.2,0.3,0.4,0.5,0.6),labels=paste(pretty(c(0,0.1,0.2,0.3,0.4,0.5,0.6))*100, " %"))

hist(unlist(circlecompletion("./raw-data/training/","./raw-data/egonets/")),
     xlab="Percentage of users labeled in a circle",main="",col="grey90",breaks=5,xaxt="n")
axis(1,at=c(0,0.2,0.4,0.6,0.8,1),labels=paste(pretty(c(0,20,40,60,80,100)), " %"))

summary(unlist(circlecount("./raw-data/training/")))
summary(circle_order("./raw-data/training/"))
sqrt(var(circle_order("./raw-data/training/")))
gini(circle_order("./raw-data/training/"))
sum(overlap("./raw-data/training/")[2:10])/sum(overlap("./raw-data/training/")[1:10])
summary(unlist(circlecompletion("./raw-data/training/","./raw-data/egonets/")))


################################COMMUNITY DETECTION#####################################################

#Run quick classifier from igraph 3 clustering)
fstcom <- leading.eigenvector.community(cc)
plot(fstcom, cc, vertex.label=NA)
fstcom <- walktrap.community(cc)
plot(fstcom, cc, vertex.label=NA)
fstcom <- fastgreedy.community(cc)
plot(fstcom, cc, vertex.label=NA)
dendPlot(fstcom, mode="phylo")

#sizes(fstcom)
#membership(fstcom)
#end

#Unsupervised MCL - Monte Carlo Clustering
A <- as.matrix(get.adjacency(cc))  #transform as matrix to run MCL
A<- A[-1,-1] #Delete the first colum andfirst row otherwise the MC will fail

pp = mcl(A, addLoops = TRUE, expansion = 2, inflation = 2, allow1 = TRUE,max.iter = 100, ESM = FALSE) #Run the MCL
#Info generated by the MCL
  pp$Cluster
  pp$n.iterations
  pp$K
#end of info

#Take the MCL resultst and inject to the group property of each vertex
circleMCL = pp$Cluster
circleMCL = append(circleMCL,0,0)
V(cc)$circle = circleMCL

SimpleGraphPlot_nolabel(cc) #plot without label for visibility

#V(cc)$circle
#vcount(cc)
#length(circleMCL)

#latent model network from the book. To be tested and build. Not working yet.
#install.packages('eigenmodel')
#library(eigenmodel)
#sizes(cc)

#A <- get.adjacency(cc, sparse=FALSE)
#cc.leig.fit1 <- eigenmodel_mcmc(A, R=2, S=11000,burn=10000)
#lat.sp.1 <- eigen(cc.leig.fit1$ULU_postmean)$vec[, 1:2]

#plot(cc, layout=lat.sp.1)

################################EVALUATING THE MODEL#################################################
# Testing BER function with cluster_leading_eigen clustering

#max (1 -ber) for each row (= each predicted community)
#each row index stands for the predicted community and each row value stands for the true community
# if length(predictedCommunities) < length(trueCommunities)
predicted <- cluster_leading_eigen(facebook_graphs_noego[[2]])

# Using method 2 = combinations (when length(predictedCommunities) > length(trueCommunities))
communitiesWithBER <- associatedCommunitiesWithBER(truecircle_list[[2]],predicted,V(facebook_graphs_noego[[2]]),2)

# Associated communities after fusion
communitiesWithBER$communities

# (1-BER) matrix after fusion 
communitiesWithBER$oneminusber

# Max BER for the associated communities
communitiesWithBER$maxber

# Using method 1 = no combinations (when length(predictedCommunities) > length(trueCommunities))
communitiesWithBER <- associatedCommunitiesWithBER(truecircle_list[[2]],predicted,V(facebook_graphs_noego[[2]]),1)


##################################
##### Performance evaluation #####
##################################

#### Leading eigen ##########
oneminusberi = c()
j=1
for(i in seq(1,110)[c(-6,-7,-28,-33,-55,-74,-86,-100,-103)]){
  print(i)
  if(!is.null(truecircle_list[[i]] )){
    predicted = cluster_leading_eigen(facebook_graphs_noego[[i]],8)

    oneminusberii = associatedCommunitiesWithBER(truecircle_list[[i]],predicted,V(facebook_graphs_noego[[i]]),2)$maxber
    oneminusberi[j] = mean(oneminusberii)
    j=j+1
  }
}
accuracyBER.leading = mean(oneminusberi)


####### MCL ###############
oneminusberi = c()
j=1
for(i in seq(1,110)[c(-6,-7,-33,-55,-65,-74,-85,-86,-95,-103,-106)]){
  print(i)
  if(!is.null(truecircle_list[[i]] )){
    predicted = facebookMCL(i)
    
    oneminusberii = associatedCommunitiesWithBER(truecircle_list[[i]],predicted,V(facebook_graphs_noego[[i]]),2)$maxber
    oneminusberi[j] = mean(oneminusberii)
    j=j+1
  }
}
accuracyBER.mcl = mean(oneminusberi)



######## WALKTRAP ############
samplegraph<-c(4, 8, 15, 28, 35, 43, 50, 56, 57, 74, 75, 76, 78, 81, 86, 87, 89, 90)


giant.component <- function(graph, ...) {
  cl <- clusters(graph, ...)
  return(induced_subgraph(graph, which(cl$membership == which.max(cl$csize))))
}

oneminusberi_walktrap_noweight = c()
j=1
for(i in samplegraph){
  print(i)
  if(!is.null(truecircle_list[[i]] )){
    predicted = walktrap.community(giant.component(facebook_graphs_noego[[i]] ))
    oneminusberii_walktrap_noweight = associatedCommunitiesWithBER(truecircle_list[[i]],predicted,V(facebook_graphs_noego[[i]]),2)$maxber
    oneminusberi_walktrap_noweight[j] = mean(oneminusberii_walktrap_noweight)
    j=j+1
  }
}
accuracyBER.walktrap_noweight = mean(oneminusberi_walktrap_noweight)

oneminusberi_walktrap_weight = c()
j=1
for(i in samplegraph){
  print(i)
  if(!is.null(truecircle_list[[i]] )){
    predicted = walktrap.community(giant.component(facebook_graphs_noego_withweights[[i]] ))
    oneminusberii_walktrap_weight = associatedCommunitiesWithBER(truecircle_list[[i]],predicted,V(facebook_graphs_noego_withweights[[i]]),2)$maxber
    oneminusberi_walktrap_weight[j] = mean(oneminusberii_walktrap_noweight)
    j=j+1
  }
}
accuracyBER.walktrap_weight = mean(oneminusberi_walktrap_weight)
