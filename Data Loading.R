#Creating a function to create an edge vector for one given path
CreateEdgeVector <- function(networkpath)  {
  egoDF<-read.csv(file = networkpath, header = FALSE, stringsAsFactors=FALSE, sep=':')
  colnames(egoDF)<-c("UserId","FriendList")
  n_ego<-unlist(regmatches(networkpath,gregexpr("[0-9]+", networkpath))) #id number of the ego
  n<-length(egoDF$UserId) #number of friends for a given ego
  edgevector<-c()
  for (j in 1:n)  {
    FriendListSplit = strsplit(egoDF$FriendList[j]," ")
    m<-length(FriendListSplit[[1]])
    edgevector=c(edgevector,n_ego,egoDF$UserId[j])  #edge between ego and its friends
    if (!is.na(FriendListSplit[[1]][2])){   #if no friends
      for (i in 2:m) {
        edgevector = c(edgevector,egoDF$UserId[j],FriendListSplit[[1]][i]) # edge between friends of the ego
      }
    }
  }
  return(edgevector)
}

#Simple plotgraph function with user id as label
SimpleGraphPlot<-function(g)    {
  plot(g, layout = layout.fruchterman.reingold,
       vertex.size = 5,
       vertex.color=V(g)$circle,
       vertex.frame.color= "white",
       vertex.label.color = "black",
       vertex.label.family = "sans",
       edge.width=0.5,  
       edge.color="grey")
}

#Simple plotgraph function with user id as label with no label
SimpleGraphPlot_nolabel<-function(g)    {
  plot(g, layout = layout.fruchterman.reingold,
       vertex.size = 10,
       vertex.color=V(g)$circle,
       vertex.label=NA,
       vertex.frame.color= "white",
       vertex.label.color = "black",
       vertex.label.family = "sans",
       edge.width=0.5,  
       edge.color="grey")
}

#Function to generate an igraph object from a given path. Simple combination of commands.
CreateGraph<-function(networkpath) {
  return(simplify(graph(CreateEdgeVector(networkpath), directed = FALSE),
                  remove.multiple = TRUE, remove.loops = TRUE,
                  edge.attr.comb = igraph_opt("edge.attr.comb")))
}

#Function to create an adjacency matrix as CSV file in a given folder from a given path, with ego number in the name
CreateAdjacency<-function(originpath,folder){
  dir.create(folder)
  n_ego<-unlist(regmatches(originpath,gregexpr("[0-9]+", originpath))) #id number of the ego
  g<-get.adjacency(CreateGraph(originpath))
  p<-paste(c("./", toString(folder), "/", "Adjacency", toString(n_ego),".csv"), collapse="")
  write.matrix(g, file = p, sep="\t")
}

#Function to turn an Adjacency matrix CSV into an igraph object again
FromCSVtoiGraph<-function(originpath){
  return(graph.adjacency(
    as.matrix(read.csv(originpath, header=TRUE, check.names=FALSE, sep="\t")),
    mode="undirected",weighted=NULL,add.colnames=c('id')))
}

loadCircles <- function(path){
  circlesFiles <- list.files(path)
  trueCircles <- vector("list", length=length(circlesFiles))
  i<-1
  for(file in circlesFiles[order(circlesFiles)]){
    pathToCircleFile <- str_replace_all(paste(path,file)," ","")
    n_col<-max(count.fields(pathToCircleFile, sep = " "))
    trueCircles[[i]] <- read.table(file = pathToCircleFile, sep=" ", fill=TRUE,col.names=1:n_col)
    # Delete the first column which stands for the circle's name
    row.names(trueCircles[[i]])<-trueCircles[[i]][,1]
    trueCircles[[i]] <- trueCircles[[i]][,-1]
    i<-i+1;
  }
  return(trueCircles)
}