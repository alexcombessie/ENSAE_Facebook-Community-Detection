
#Run quick classifier from igraph 3 clustering)
cc <- facebook_graphs_noego[[2]]
training <- truecircle_list[[2]]
nrow(training)
head(graph_list[[10]])
SimpleGraphPlot_nolabel(cc)

##Leading eigen vector
matrixcount = c()
for (i in 0:50){
  fstcom <- leading.eigenvector.community(cc,i)
  matrixcount[i+1] = length(table(fstcom$membership))
}
data<- data.frame(grp=matrixcount, occ=seq(0:50))
ggplot(data=data,aes(x=occ,y=grp)) + 
  geom_bar( stat="identity", colour = "lightblue", size = 0.2)+
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.3))+
  xlab("Number of iterations") +
  ylab("Number of communities") + 
  ggtitle("Leading eigen vector community convergence")

fstcom <- leading.eigenvector.community(cc,5)
length(table(fstcom$membership))
matrixcount[1]
plot(fstcom, cc, vertex.label=NA)

ber(training,fstcom,V(cc))


##walktrap
fstcom <- walktrap.community(cc, steps = 1)
plot(fstcom, cc, vertex.label=NA)

matrixcount = c()
for (i in 1:50){
  fstcom <- walktrap.community(cc, steps = i)
  matrixcount[i+1] = length(table(fstcom$membership))
}
data<- data.frame(grp=matrixcount, occ=seq(0:50))
ggplot(data=data,aes(x=occ,y=grp)) + 
  geom_bar( stat="identity", colour = "lightblue", size = 0.2)+
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.3))+
  xlab("Number of iterations") +
  ylab("Number of communities") + 
  ggtitle("Walktrap community convergence")
matrixcount


fstcom <- walktrap.community(cc)
table(fstcom$membership)
length(table(fstcom$membership))
plot(fstcom, cc, vertex.label=NA)

#fast greedy
fstcom <- fastgreedy.community(cc)
plot(fstcom, cc, vertex.label=NA)
dendPlot(fstcom, mode="phylo")
table(fstcom$membership)


#sizes(fstcom)
#membership(fstcom)
#end

#Unsupervised MCL - Monte Carlo Clustering
set.seed(41)
#Unsupervised MCL - Monte Carlo Clustering
clustering_MCL = function (cGraph)
{
  A <- as.matrix(get.adjacency(cGraph))  #transform as matrix to run MCL
  A<- A[-1,-1] #Delete the first colum andfirst row otherwise the MC will fail
  pp = mcl(A, addLoops = TRUE, expansion = 100, inflation = .2, allow1 = TRUE,max.iter = 10000, ESM = FALSE) #Run the MCL
  circleMCL = pp$Cluster
  circleMCL = append(circleMCL,0,0)
  return(circleMCL)
}
#cc = graph_list[[2]]
V(cc)$circle = clustering_MCL(graph_list[[2]])
length(table(V(cc)$circle))
SimpleGraphPlot_nolabel(cc) #plot without label for visibility