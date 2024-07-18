#code for finding a multi-stage drop-the-losers design with given FWER and power.

#Example:
#finddesign(J=3,K=4,ns=c(4,2,1),delta1=0.545,delta0=0.178,requiredFWER=0.05,requiredpower=0.9,treatmentsigmas=c(1,1,1,100,1))



library(mvtnorm)

#integral_typeIerror is a function called by finddesign that finds the difference between the FWER of a given design and the required FWER


integral_typeIerror=function(c,requiredtypeIerror,mean,var,combinations,K)
{

lower=rep(0,length(mean))
lower[length(lower)]=c

upper=rep(Inf,length(mean))

int=pmvnorm(lower=lower,upper=upper,mean=mean,sigma=var)

return(as.double(int)*factorial(K)-requiredtypeIerror)
}


#integral_power is a function called by finddesign that finds the difference between the power under the LFC of a given design and the required power


integral_power=function(n,requiredpower,c,A,var,combinations,J,K,delta1,delta0,cumgroupsizes)
{
#get mu under LFC (assuming standardised data
mu=rep(0,length(A[1,]))
for(i in 1:(K-1))
{
mu[((i-1)*J+1):(i*J)]=sqrt(n*cumgroupsizes/2)*delta0
    
}


mu[((K-1)*J+1):(K*J)]=sqrt(n*cumgroupsizes/2)*delta1


mean=as.double(A%*%mu)

lower=rep(0,length(mean))
lower[length(lower)]=c

upper=rep(Inf,length(mean))

int=pmvnorm(lower=lower,upper=upper,mean=mean,sigma=var)

return(as.double(int)*factorial(K-1)-requiredpower)
}




#finddesign returns a design with specified family-wise error rate and power. Arguments:

#J - number of stages
#K - number of experimental arms
#ns - vector of length J stating the number of treatments at each stage - first entry must be K and last entry must be 1, and entries must be strictly decreasing
#delta1 - standardised clinically relevant difference for power
#delta0 - standardised uninteresting treatment threshold for power
#requiredfwer - required family-wise error rate at the global null
#requiredpower - required power at the least favourable configuration
#groupsizes - relative proportion of patients recruited per arm at each stage (e.g. (1,1,1) for J=3 would mean equal numbers recruited per arm per stage:
#treatmentsigmas - vector of length K+1 with standard deviation of effect for each arm (1 is control, j+1th entry is jth experimental arm)

#output is a list with entries n, c and totalSS. n is the required group-size, c is the required critical value and totalSS is the totalSS used by the design.

finddesign=function(J,K,ns,delta1,delta0,requiredfwer,requiredpower,groupsizes,treatmentsigmas)
{


if(J<2 || J>10)
{
stop("J must be an integer between 2 and 10")
}

if(K<2 || K>100)
{
stop("K must be an integer between 2 and 100")
}

#test ns:

if(length(ns)!=J)
{
stop("ns must be a vector of length J")
}

if(ns[1]!=K || ns[length(ns)]!=1 || min(ns[-length(ns)]-ns[-1])<1)
{
stop("ns must have first entry K, last entry 1 and each entry be strictly less than the previous one")
}
if(requiredfwer<1e-8 || requiredfwer>1-1e-8)
{
stop("requiredfwer must be strictly between 0 and 1")
}


if(requiredpower<1e-8 || requiredpower>1-1e-8)
{
stop("requiredpower must be strictly between 0 and 1")
}


cumgroupsizes=cumsum(groupsizes)

#normalise so that first stage is 1:

cumgroupsizes=cumgroupsizes/cumgroupsizes[1]


sigma=matrix(0,J*K,J*K)

#fill in sigma blocks

for(i in 1:J)
{
for(j in 1:J)
{
sigma[i,j]=sqrt(min(cumgroupsizes[i],cumgroupsizes[j])/max(cumgroupsizes[i],cumgroupsizes[j]))
}
}

for(i in 1:K)
{
sigma[((i-1)*J+1):(i*J),((i-1)*J+1):(i*J)]=sigma[1:J,1:J]

for(j in (1:K)[-i])
{
sigma[((i-1)*J+1):(i*J),((j-1)*J+1):(j*J)]=sigma[1:J,1:J]*sqrt((treatmentsigmas[1]^4)/((treatmentsigmas[i+1]^2+treatmentsigmas[1]^2)*(treatmentsigmas[j+1]^2+treatmentsigmas[1]^2)))
     }

}


rowsofA=ns-1
rowsofA[length(rowsofA)]=1

A=matrix(0,sum(rowsofA),length(sigma[,1]))

treatments=1:K
whichstagedropped=rep(1:J,times=c(ns[1:(J-1)]-ns[-1],1))

#go through each stage; each treatment that drops out should have a row in A together with each treatment dropping out in a subsequent stage

tempint=0

for(i in 1:(J-1))
{

treatments_thisstage=treatments[which(whichstagedropped==i)]
treatments_futurestage=treatments[which(whichstagedropped>i)]

#each treatment in treatments_futurestage must beat last entry in treatments_thistage

for(k1 in treatments_futurestage)
{
tempint=tempint+1
A[tempint,(treatments_thisstage[length(treatments_thisstage)]-1)*J+i]=-1
A[tempint,(k1-1)*J+i]=1
}

#each treatment in treatments_thisstage must beat the one below it:
#check if more than one treatment is dropped
if(length(treatments_thisstage)>1)
{


for(k2 in 1:(length(treatments_thisstage)-1))
{

tempint=tempint+1
A[tempint,(treatments_thisstage[k2]-1)*J+i]=-1
A[tempint,(treatments_thisstage[k2+1]-1)*J+i]=1



}

}



}

A[length(A[,1]),length(A[1,])]=1

#get mean vector under HG:

mu=rep(0,length(sigma[,1]))
mean=as.double(A%*%mu)

var=A%*%sigma%*%t(A)




c=uniroot(integral_typeIerror,lower=-2,upper=5,requiredtypeIerror=requiredfwer,mean=mean,var=var,combinations=combinations,K=K)$root

#find sample size for given power power



n=uniroot(integral_power,lower=0,upper=2000,requiredpower=requiredpower,c=c,A=A,var=var,combinations=combinations,J=J,K=K,delta1=delta1,delta0=delta0,cumgroupsizes=cumgroupsizes)$root


return(list(n=n,c=c,totalSS=n*sum((ns+1)*groupsizes/groupsizes[1])))

}



finddesign(J=3,K=4,ns=c(4,2,1),delta1=0.545,delta0=0.178,requiredfwer=0.05,requiredpower=0.9,treatmentsigmas=c(1,1,1,100,1))






