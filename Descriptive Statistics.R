#Function to plot the degree distribution in log scale
log_degree_distribution = function(g) {
  d<-degree(g, mode = "all")
  degree_dist<-degree.distribution(g, mode = "all", cumulative = FALSE)
  degree<-1:max(d)
  proba<-degree_dist[-1]
  nonzero <- which(proba != 0)  # delete blank values
  proba <- proba[nonzero]
  degree <- degree[nonzero]
  plot(proba ~ degree, log = "xy", xlab = "Degree distribution (log)", 
       ylab = "Probability (log)", col = 1)
}

#Adaptation of the previous function to plot all the networks and fitting a power law
log_degree_powerlaw = function(g_list) {
  proba_all<-c()
  degree_all<-c()
  #Loop on all the graphs
  for (i in 1:length(g_list)){
    d<-degree(g_list[[i]], mode = "all")
    degree_dist<-degree.distribution(g_list[[i]], mode = "all", cumulative = FALSE)
    degree<-1:max(d)
    proba<-degree_dist[-1]
    nonzero <- which(proba != 0)  # delete blank values
    proba <- proba[nonzero]
    degree <- degree[nonzero]
    proba_all<-c(proba_all,proba)
    degree_all<-c(degree_all,degree)
  }
  #Fit a linear regression model
  regression <- lm(log(proba_all) ~ log(degree_all))
  coeff <- coef(regression)
  powerlaw <- function(x){exp(coeff[[1]] + coeff[[2]] * log(x))}
  print(paste("Alpha =", round(-coeff[[2]], 3)))
  print(paste("R square =", round(summary(regression)$r.squared, 3)))
  #plot the distribution 
  plot(proba_all ~ degree_all, log = "xy", main="",
       xlab = "Degree distribution (log)", ylab = "Frequency (log)", col = "darkgrey")
  curve(powerlaw, col = "red", add = TRUE, lwd=2)
}

#Function to return the regression object
powerlaw_reg = function(g_list) {
  proba_all<-c()
  degree_all<-c()
  #Loop on all the graphs
  for (i in 1:length(g_list)){
    d<-degree(g_list[[i]], mode = "all")
    degree_dist<-degree.distribution(g_list[[i]], mode = "all", cumulative = FALSE)
    degree<-1:max(d)
    proba<-degree_dist[-1]
    nonzero <- which(proba != 0)  # delete blank values
    proba <- proba[nonzero]
    degree <- degree[nonzero]
    proba_all<-c(proba_all,proba)
    degree_all<-c(degree_all,degree)
  }
  #Fit a linear regression model
  regression <- lm(log(proba_all) ~ log(degree_all))
  return(regression)
}

#Function to count the number of lines in each file, for a specific directory. 
#Can be used to count the number of circles or friends for egonets
circlecount<-function(path){
  circlesFiles <- list.files(path)
  count<-vector("list", length=length(circlesFiles))
  i<-1
  for(file in circlesFiles){
    pathToCircleFile = str_replace_all(paste(path,file)," ","")
    id<-unlist(regmatches(pathToCircleFile,gregexpr("[0-9]+", pathToCircleFile)))
    count[[id]]<-length(readLines(pathToCircleFile))
    i<-i+1
  }
  return(count)
}

#Function to count the number of friends of each circle for all networks
circle_order<-function(path){
  numfriends<-c()
  for(x in loadCircles(path)){
    numfriends<-c(numfriends,apply(x, 1, function(x) length(which(!is.na(x)))))
  }
  return(numfriends)
}

#Function to return the set of friends in social circles
set_friends<-function(path){
  circlesFiles <- list.files(path)
  set <- vector("list", length=length(circlesFiles))
  i<-1
  for(file in circlesFiles){
    pathToCircleFile = str_replace_all(paste(path,file)," ","")
    id<-unlist(regmatches(pathToCircleFile,gregexpr("[0-9]+", pathToCircleFile)))
    for(l in readLines(pathToCircleFile)){
         set[[id]]<-c(set[[id]],tail(unlist(str_split(l," ")),-1))
       }
    i<-i+1
  }
  return(set)
}

#Function to count the number of friends belonging to k circles, from 1 to 10
overlap<-function(path){
  t<-c()
  t2<-c()
  set<-set_friends(path)
  for(l in set){
    t<-c(t, table(l))
  }
  for(i in 1:10){
    t2<-c(t2,length(which(t==i)))
  }
  return(t2)
}

#Function to see how many friends have been labeled to a social circle
circlecompletion<-function(path_training,path_egonets){
  countunique<-lapply(set_friends(path_training),function(x) length(unique(x)))
  numfriends<-circlecount(path_egonets)
  keys<-c(unique(names(countunique)))
  combination<-setNames(mapply(c, countunique[keys], numfriends[keys]), keys)
  return(lapply(combination,function(x) x[1]/x[2]))
}



