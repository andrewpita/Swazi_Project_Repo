

library(dplyr)
library(BiasedUrn)
library(truncnorm)
library(coda)
library(stringr)
library(hash)
library(gdata)

loc_hash = hash(c("ME", "MM"), c("Mbabane","Manzini"))

results = read.csv("summary.csv", row.names = 1)
nmes = rownames(results)

results = results[,"Mean"]
names(results) = nmes

make_city = function(string) {
  
  #grepl returns a vector of logical values, true
  #when a parameter name matches the string argument.
  #The $ means that we want a match at the end 
  indeces = grepl(paste(string,"$", sep = ""), names(results))
  
  city = results[indeces]
  names(city)=NULL
  return(city)
}


ME = make_city("ME")
MM = make_city("MM")

#true Mbabane parameters
true.N.ME = round(ME[6])

true.p.srv.ME = ME[1]
true.p.uid.ME = ME[2]
true.p.cpn.ME = ME[3]
#true.N.cpn.ME = ME[4]


#true Manzini parameters
true.N.MM = round(MM[6])

true.p.srv.MM = MM[1]

true.p.uid.MM = MM[2]

true.p.cpn.MM = MM[3]

#names of the unknown parameters that we iteratively simulate in 
#the Gibbs sampler

parnames=c("p.srv.ME","p.uid.ME","p.cpn.ME","N.cpn.ME","r.ME","N.ME",
           "p.srv.MM","p.uid.MM","p.cpn.MM","N.cpn.MM","r.MM","N.MM")


probpars = parnames[which(startsWith(parnames,"p"))]

sim_summary = function(x) {
  
  
  summ=c(min(x),quantile(x,c(0.01,0.025,0.05,0.1,0.2,0.25)),
         median(x),mean(x), sd(x),quantile(x,1-rev(c(0.01,0.025,0.05,0.1,0.2,0.25))),max(x))
  names(summ)=c("min","1%","2.5%","5%","10%","20%","25%","Median",
                "Mean", "sd","75%","80%","90%","95%","97.5%","99%","max")
  summ
  
}

#function that performs the gibb sampling algorithm. Begins by simulating
#data 
sim_data_gibbs_samp = function(seed, N = 10000, a.prior = 1, b.prior = 1) {
  
  ### Simulate data for Mbabane/Ezulwini ###
  
  ME.sample = rmultinom(1, size = true.N.ME, prob = 
                          c(true.p.srv.ME*(1 - true.p.uid.ME)*(1 - true.p.cpn.ME), #only in survey [1]
                            true.p.srv.ME*true.p.uid.ME*(1- true.p.cpn.ME), #in only survey and uid [2]
                            true.p.srv.ME*true.p.uid.ME * true.p.cpn.ME, #in all three [3]
                            true.p.uid.ME*(1 - true.p.srv.ME)*(1 - true.p.cpn.ME), #in only uid [4]
                            true.p.uid.ME*(1 - true.p.srv.ME)*true.p.cpn.ME, #in uid and cpn [5]
                            true.p.cpn.ME*(1 - true.p.srv.ME)*(1 - true.p.uid.ME), #in only cpn [6]
                            true.p.cpn.ME*true.p.srv.ME*(1 - true.p.uid.ME),#in cpn and srv [7]
                            (1- true.p.srv.ME)*(1 - true.p.uid.ME)* (1 - true.p.cpn.ME)#in none [8]
                          ))
  
  ####unknown
  #N.srv.uid.cpn.ME = ME.sample[3]
  true.N.cpn.ME = as.numeric(ME.sample[6] + ME.sample[7] + ME.sample[5] + ME.sample[3])
  ####
  
  #storing simulated data
  N.srv.uid.ME = ME.sample[2] + ME.sample[3]
  N.srv.cpn.ME = ME.sample[7] + ME.sample[3]
  N.srv.ME = ME.sample[1] + ME.sample[2] + ME.sample[3] + ME.sample[7]
  N.uid.ME = ME.sample[4] + ME.sample[2] + ME.sample[3] + ME.sample[5]
  true.r.ME = as.numeric(true.N.ME - ME.sample[8])
  
  ### Simulate data for Manzini/Matsapha ###
  
  MM.sample = rmultinom(1, size = true.N.MM, prob = 
                          c(true.p.srv.MM*(1 - true.p.uid.MM)*(1 - true.p.cpn.MM), #only in survey [1]
                            true.p.srv.MM*true.p.uid.MM*(1- true.p.cpn.MM), #in only survey and uid [2]
                            true.p.srv.MM* true.p.uid.MM * true.p.cpn.MM, #in all three [3]
                            true.p.uid.MM*(1 - true.p.srv.MM)*(1 - true.p.cpn.MM), #in only uid [4]
                            true.p.uid.MM*(1 - true.p.srv.MM)*true.p.cpn.MM, #in uid and cpn [5]
                            true.p.cpn.MM*(1 - true.p.srv.MM)*(1 - true.p.uid.MM), #in only cpn [6]
                            true.p.cpn.MM*true.p.srv.MM*(1 - true.p.uid.MM),#in cpn and srv [7]
                            (1- true.p.srv.MM)*(1 - true.p.uid.MM)* (1 - true.p.cpn.MM)#in none [8]
                          ))
  
  
  ####unknown
  #N.srv.uid.cpn.MM = MM.sample[3]
  true.N.cpn.MM = as.numeric(MM.sample[6] + MM.sample[7] + MM.sample[5] + MM.sample[3])
  ####
  
  
  #storing simulated data
  N.srv.uid.MM = MM.sample[2] + MM.sample[3]
  N.srv.cpn.MM = MM.sample[7] + MM.sample[3]
  N.srv.MM = MM.sample[1] + MM.sample[2] + MM.sample[3] + MM.sample[7]
  N.uid.MM = MM.sample[4] + MM.sample[2] + MM.sample[3] + MM.sample[5]
  true.r.MM = as.numeric(true.N.MM - MM.sample[8])
  
  ## Corridor: sum of MM and ME
  N.cpn.Co = true.N.cpn.MM + true.N.cpn.ME
  

  #empty matrix to store values for parameters at each iteration
  M = matrix(nrow = N, ncol = length(parnames))
  colnames(M) = parnames
  
  #constants for beta posteriors
  a.srv = a.uid = a.rnb = a.cpn = a.prior
  b.srv = b.uid = b.rnb = b.cpn = b.prior
  
  for (par in probpars) {
    
    #par = 0.1
    assign(par,0.1)
  }
  
  N.ME = round(N.srv.ME/p.srv.ME)
  N.MM = round(N.srv.MM/p.srv.MM)
  N.cpn.ME = max(N.srv.cpn.ME, round(N.cpn.Co*N.ME/(N.ME + N.MM)))
  N.cpn.MM = N.cpn.Co - N.cpn.ME
  
  #begin Gibbs sampler algorithm
  set.seed(seed)
  
  for (i in 1:N) {
    
    ### Updating Mbabane/Ezulwini
    
    ov.ME = rhyper(1, N.uid.ME - N.srv.uid.ME, N.ME - N.srv.ME - 
                     N.uid.ME + N.srv.uid.ME, N.cpn.ME - N.srv.cpn.ME)
    
    r.ME = N.srv.ME + N.uid.ME + N.cpn.ME - N.srv.uid.ME - 
      N.srv.cpn.ME - ov.ME
    
    x = (1 - p.srv.ME)*(1 - p.uid.ME) * (1 - p.cpn.ME)
    
    N.ME = r.ME + rnbinom(1, r.ME, 1 - x)
    
    p.srv.ME = rbeta(1, a.srv + N.srv.ME, b.srv + N.ME - N.srv.ME)
    
    p.uid.ME = rbeta(1, a.uid + N.uid.ME, b.uid + N.ME - N.uid.ME)
    
    p.cpn.ME = rbeta(1, a.cpn + N.cpn.ME, b.cpn + N.ME - N.cpn.ME)
    
    ### Updating Manzini/Matsapha
    
    ov.MM = rhyper(1, N.uid.MM - N.srv.uid.MM, 
                   N.MM - N.srv.MM - N.uid.MM + N.srv.uid.MM, 
                   N.cpn.MM - N.srv.cpn.MM)
    
    r.MM = N.srv.MM + N.uid.MM + N.cpn.MM - N.srv.uid.MM - 
      N.srv.cpn.MM - ov.MM
    
    x = (1 - p.srv.MM)*(1 - p.uid.MM) * (1 - p.cpn.MM)
    
    N.MM = r.MM + rnbinom(1, r.MM, 1- x)
    
    p.srv.MM = rbeta(1, a.srv + N.srv.MM, b.srv + N.MM - N.srv.MM)
    
    p.uid.MM = rbeta(1, a.uid + N.uid.MM, b.uid + N.MM - N.uid.MM)
    
    p.cpn.MM = rbeta(1, a.cpn + N.cpn.MM, b.cpn + N.MM - N.cpn.MM)
    
    ### Updating coupon distribution in corridor
    
    #we define a variable z = N.cpn.ME - N.srv.cpn.ME
    # - ov.ME. The posterior distribution for z 
    # is fisher's noncentral hypegeometric
    
    N.cpn.ME = N.srv.cpn.ME + ov.ME + rFNCHypergeo(1, 
                                                   N.ME - N.srv.ME - N.uid.ME + N.srv.uid.ME,
                                                   N.MM - N.srv.MM - N.uid.MM + N.srv.uid.MM,
                                                   N.cpn.Co - N.srv.cpn.ME - ov.ME - 
                                                     N.srv.cpn.MM - ov.MM,
                                                   ( p.cpn.ME/(1-p.cpn.ME) ) / (( p.cpn.MM) / 
                                                                                  (1 - p.cpn.MM) ) )
    
    N.cpn.MM = N.cpn.Co - N.cpn.ME
    
    
    for (par in parnames) {
      
      M[i, par] = get(par)
      
    }
    
    if ( (i %% 1) == 200) {
      
      print(i)
    }
    
  }
  
  #make a matrix containing the estimates for Hhohho and 
  #Shiselweni, and then join it to our original matrix
  #M1 = cbind(M[,"N.ME"] + M[,"N.PP"], M[,"N.Nh"] + M[,"N.Lv"])
  #colnames(M1) = c("N.Hh", "N.Sw")
  #M2=cbind(M, M1)
  
  #apply the sim_summary function to each column of 
  #the gibbs samples
  summary=apply(M,2,sim_summary)
  
  #get the true values. Some of the true values
  #were defined in the global environment at the beginning of 
  #the script. Others can only be set once the data has 
  #been simulated. envir = sys.frame(sys.parent(0)) is an 
  #argument to get that specifies we want the true values
  #that have been defined within this function
  
  truth=sapply(paste0("true.",parnames),get,envir=sys.frame(sys.parent(0)), simplify=FALSE) 
  names(truth)=colnames(summary)
  
  #a 1 or 0 value for each parameter that tells whether or 
  #not it was contained in the 95% confidence interval of the 
  #gibbs samples
  cov=matrix((summary["2.5%",]< truth)*(summary["97.5%",]>truth),nrow=1)
  row.names(cov)="coverage"
  
  #burn the first 5,000 values so that a posterior distribution
  #can be constructed from multiple iterations
  #of the gibbs sampler
  M.burn = M[5000:10000,]
  
  #return the values that weren't burned, the summary and 
  #coverage values, and the true values for this gibbs sampler
  return.list = list(M.burn, rbind(summary,cov), truth)
  
  return(return.list)
  
  
}

#run the gibbs sampler nseed times
nseed=400
seeds=.Random.seed[1:nseed]
Mlist=lapply(seeds,function(seed) {print(seed);sim_data_gibbs_samp(seed)})

#construct the posterior distribution and 
#construct a matrix of the true values from each
#run of the gibbs sampler. Later we average
#these values 

posterior = Mlist[[1]][[1]]
truth.vals = unlist(Mlist[[1]][[3]])

for (i in 2:length(Mlist)) {
  
  posterior = rbind(posterior, Mlist[[i]][[1]])
  truth.vals = rbind(truth.vals, unlist(Mlist[[i]][[3]]))
}


#average the true values so that we have an 
#average truth to compare to the posterior mean
truth = apply(truth.vals, 2, mean)

#get the summary and coverage matrices from Mlist
sumlist=lapply(Mlist, "[[",2)

#average them. This gives the coverage percentages
avgsum = Reduce("+",sumlist)/nseed
avgsum = rbind(avgsum, truth)

#function that plots the density of the posterior distribution
mydensityplot=function(col,posterior){
  
  colname=str_replace_all(col,"[.]","_")
  imagename=paste0("results/densityplots/",colname,".png")
  png(imagename)
  plot(density(posterior[,col]),type="l",
       main=loc_hash[[str_sub(col,-2)]],xlab="x",ylab="")
  abline(v = truth[col], col = "blue")
  abline(v = mean(posterior[,col]), col = "black")
  legend(x = "topright",legend = c("truth", "post mean"), 
         col = c("blue", "black"), lty=1, cex=0.8)
  dev.off()
  
}

#create directories for the results and plots

dir.create("results")
dir.create("results/densityplots")


sapply(colnames(posterior),mydensityplot,posterior)

write.csv(avgsum,paste("output_","summary.csv",sep = ""))








