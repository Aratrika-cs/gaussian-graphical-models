# Flow cytometry data imported from sparsebn package in R 

library(sparsebn)
data = cytometryContinuous$data
colnames(data)=c("Raf","Mek","Plcg","PIP2","PIP3","Erk","Akt","PKA","PKC","P38","Jnk")

# Applying Graphical Lasso algorithm for a grid of penalty parameter values 

library(glasso)
rho =c(0,0.1,0.2,0.5,0.7,1)
glasso_fits = lapply(rho,function(x)glasso(cov(data),x,nobs=nrow(data)))

# Function to calculate L1 norm of estimated covaraince matrix

L1_norm = function(matrix)
{
  return(sum(abs(matrix)))
}

# Computation of L1 norm of estimated covariance matrix 
lapply(glasso_fits,function(x)L1_norm(x$wi))

# Ploting the estimated graph structures for a grid of penalty parameter values

library(network)

for(i in 1:6)
{
  P <- round(glasso_fits[[i]]$wi,10)
  colnames(P)=colnames(data)
  A <- ifelse(P!=0 & row(P)!=col(P),1,0)
  net=network(A)
  plot(net,displaylabels=T,mode="circle",vertex.cex=3,usearrows=F)
}