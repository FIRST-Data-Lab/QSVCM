
# Quantile loss function:

rhotau=function(u,tau){
  return(u*(tau-(u<0)))
}

