ber = function(trueCircles,computedCircles,allVertices){
  numberComputedCircles = length(computedCircles)
  numberTrueCircles =nrow(trueCircles)
  allIds=allVertices$id
  ber = matrix(data = rep(0,numberTrueCircles*numberComputedCircles),nrow = numberComputedCircles,ncol=numberTrueCircles)
  
  ber = apply(trueCircles,1,trueCircleBer,allIds,computedCircles)
  return(ber)
}

trueCircleBer = function(trueCircle,allIds,computedCircles){
  # clean the true circle from NA values
  
  trueCircle = trueCircle[!is.na(trueCircle) ]
  # Compute the complementary values from the true circle
  ber=c()
  complementaryTrue = setdiff(allIds,trueCircle)
  if(length(groups(computedCircles))==0){
    ber=lapply(computedCircles,computedCircleBer,trueCircle=trueCircle,allIds=allIds,complementaryTrue=complementaryTrue)
  }else{
    ber=lapply(groups(computedCircles),computedCircleBer,trueCircle=trueCircle,allIds=allIds,complementaryTrue=complementaryTrue)
  }
  return(unlist(ber))
}

computedCircleBer = function(computedCircle,trueCircle,allIds,complementaryTrue){
  # Cpredicted\Ctrue
  diff = setdiff(computedCircle,trueCircle)
  
  complementaryPredicted = setdiff(allIds,computedCircle)
  
  diffComplementary = setdiff(complementaryPredicted,complementaryTrue)
  
  ber= 0.5 *(length(diff)/length(computedCircle) + length(diffComplementary)/length(complementaryPredicted))
  return(ber)
}

# Function associating the true circles and ones calculated (predicted), highlighting the corresponding max BER 
associatedCommunitiesWithBER = function(trueCircles,predicted,allVertices,...){
  method= as.numeric(...)
  oneminusber = 1- ber(trueCircles,predicted,allVertices)
  maxBer=apply(oneminusber,1,max)
  summaxBer = sum(maxBer)/length(predicted)
  associatedCommunities = apply(oneminusber,1,which.max)
  
  if(nrow(trueCircles) < length(predicted) & method==1){
    predictedFusion = list()
    j=1
    for(i in unique(associatedCommunities)){
      toBeFusionned = which(associatedCommunities==i)
      predictedFusion[[j]]=Reduce(union,predicted[toBeFusionned])
      j=j+1
    }
    oneminusber=1- ber(trueCircles,predictedFusion,allVertices)
    maxBer = sum(apply(oneminusber,1,max))/length(predicted)
    associatedCommunities = apply(oneminusber,1,which.max)
  }
  if(nrow(trueCircles) < length(predicted) & method==2){
    allPredictedCombinations = combinationsOfCircles(predicted,trueCircles)
    summaxBer=c()
    oneminusber_list=list()
    for(i in seq(1,length(allPredictedCombinations))){
      oneminusber_list[[i]] = 1- ber(trueCircles,allPredictedCombinations[[i]],allVertices)
      maxBer=apply(oneminusber_list[[i]],1,max)
      summaxBer[i] = sum(maxBer)/length(predicted)
    }
    imax = which(summaxBer==max(summaxBer))
    if(length(imax)>1){
      imax = imax[1]
    }
    oneminusber=oneminusber_list[[imax]]
    associatedCommunities = apply(oneminusber,1,which.max)
  }
  
  
  return(list(communities=associatedCommunities,oneminusber=oneminusber,maxber=maxBer))
}

combinationsOfCircles = function(predicted,trueCircles){
  result=fusionCombinations(length(predicted),nrow(trueCircles))
  return(lapply(result,fusion,predicted=predicted))
}

fusion= function(x,predicted){
  resulta =list()
  for(i in seq(1,length(x))){
    resulta = append(resulta,list(Reduce(union,predicted[x[[i]]])))
  }
  return(resulta)
}

f <- function(i,set,m.start,m.end) {
  L <- mapply(seq,m.start[,i],m.end[,i])
  lapply(L,function(i) if (any(i>length(set))) NA else set[i])
}

fusionCombinations= function(predictedLength,trueLength){
  predictedSeq <- seq(1,predictedLength)
  m<- compositions(predictedLength)
  m.end <- apply(m,2,cumsum)
  m.start <- rbind(row(m.end)[1,],m.end+1)[1:nrow(m.end),]
  a=lapply(lapply(1:ncol(m),f,predictedSeq,m.start,m.end),function(x) if(sum(is.na(unlist(x)))==(predictedLength-trueLength)) return(x[1:trueLength]))
  return(a[lapply(a,length)>0])
}

facebookMCL = function(i){
  A <- as.matrix(get.adjacency(facebook_graphs_noego[[i]]))  #transform as matrix to run MCL
  A<- A[-1,-1] #Delete the first colum andfirst row otherwise the MC will fail
  
  pp = mcl(A, addLoops = TRUE, expansion = 2, inflation = 2, allow1 = FALSE,max.iter = 200, ESM = FALSE)
  ids=V(facebook_graphs_noego[[i]])$id         
  clusters= sapply(unique(pp$Cluster),function(x) ids[which(pp$Cluster==x)] )
  return(clusters)
}