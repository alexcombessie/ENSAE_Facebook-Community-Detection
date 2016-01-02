splitLine = function(x){
  return (strsplit(as.character(x[x!='']),';'))
}

toKeyValue = function(x){
  x= x[[1]]
  return (c(paste(x[-length(x)],collapse=';'),x[length(x)]))
}

parseFeatures = function(x){
  x = as.array(x[[1]])
  return (apply(x,1,toKeyValue))
}


#weighting

similarity<-function(g,e){
  edge_vector<-ends(g,e)
  list_attributes1<-lapply(sapply(allAttributes,function(x) get.vertex.attribute(g,x,edge_vector[1])),function(x) strsplit(x,split=","))
  list_attributes2<-lapply(sapply(allAttributes,function(x) get.vertex.attribute(g,x,edge_vector[2])),function(x) strsplit(x,split=","))
  sim<-rep(NA,length(list_attributes1))
  for(i in 1:length(list_attributes1)){
    sim[i]<-any(unlist(list_attributes1[[i]]) %in% unlist(list_attributes2[[i]])) && !is.na(list_attributes1[[i]]) && !is.na(list_attributes2[[i]])
  }
  return(sum(sim))
}
