

library(dplyr)
library(BiasedUrn)
library(truncnorm)
library(coda)
library(stringr)
library(hash)
library(gdata)

loc_hash = hash(c("Nh"), c("Nhlangano"))

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


Nh = make_city("Nh")

true.N.Nh = round(Nh[5])

true.p.srv.Nh = Nh[1]
true.p.uid.Nh = Nh[2]
true.p.rnb.Nh = Nh[3]
#this time we actually save the true
#r.Nh
true.r.Nh = round(Nh[4])



parnames = c("p.srv.Nh", "p.uid.Nh","p.rnb.Nh", "r.Nh", "N.Nh")


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
  
  ### Simulate Nhlangano data ###
  
  Nh.sample = rmultinom(1, size = true.r.Nh, prob = 
                          c(true.p.srv.Nh*(1 - true.p.uid.Nh)*(1 - true.p.rnb.Nh), #only in survey [1]
                            true.p.srv.Nh*true.p.uid.Nh*(1- true.p.rnb.Nh), #in only survey and uid [2]
                            true.p.srv.Nh*true.p.uid.Nh * true.p.rnb.Nh, #in all three [3]
                            true.p.uid.Nh*(1 - true.p.srv.Nh)*(1 - true.p.rnb.Nh), #in only uid [4]
                            true.p.uid.Nh*(1 - true.p.srv.Nh)*true.p.rnb.Nh, #in uid and rainbow [5]
                            true.p.rnb.Nh*(1 - true.p.srv.Nh)*(1 - true.p.uid.Nh), #in only rnb [6]
                            true.p.rnb.Nh*true.p.srv.Nh*(1 - true.p.uid.Nh)#in rnb and srv [7]
                            #we leave out the bin of people in none, because we are drawing from
                            #the people who are in at least one survey
                            #(1- true.p.srv.Nh)*(1 - true.p.uid.Nh)* (1 - true.p.rnb.Nh)#in none [8]
                          ))
  
  ####unknown in the context of the data
  #N.srv.uid.rnb.Nh = Nh.sample[3]
  ####
  
  #storing simulated data
  N.srv.uid.Nh = Nh.sample[2] + Nh.sample[3]
  N.srv.rnb.Nh = Nh.sample[7] + Nh.sample[3]
  N.srv.Nh = Nh.sample[1] + Nh.sample[2] + Nh.sample[3] + Nh.sample[7]
  N.uid.Nh = Nh.sample[4] + Nh.sample[2] + Nh.sample[3] + Nh.sample[5]
  N.rnb.Nh = Nh.sample[6] + Nh.sample[7] + Nh.sample[5] + Nh.sample[3]
  Nh.in.none = true.N.Nh - true.r.Nh
  #true.r.Nh = as.numeric(true.N.Nh - Nh.sample[8] )
  
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
  
  N.Nh = round(N.srv.Nh/p.srv.Nh)
  
  #begin Gibbs sampler algorithm
  set.seed(seed)
  
  for (i in 1:N) {
    
    ### updating Nhlangano
    
    #r.Nh is the number of individuals who participated in 
    #at least 1 survey.  This was known in Lavumisa and Piggs Peak
    #but is unknown here. We draw from the posterior for ov.Nh, 
    #the unknown N.uid.rnb.Nh, which is hypergeometric
    
    r.Nh = N.srv.Nh + N.uid.Nh + N.rnb.Nh - N.srv.uid.Nh - N.srv.rnb.Nh - 
      rhyper(1, N.uid.Nh - N.srv.uid.Nh, 
             N.Nh - N.srv.Nh - N.uid.Nh + N.srv.uid.Nh, 
             N.rnb.Nh - N.srv.rnb.Nh)
    
    x = (1 - p.srv.Nh) * ( 1 - p.uid.Nh) * (1 - p.rnb.Nh)
    
    N.Nh = r.Nh + rnbinom(1, r.Nh, 1 - x)
    
    p.srv.Nh = rbeta(1, a.srv + N.srv.Nh, b.srv + N.Nh - N.srv.Nh)
    p.uid.Nh = rbeta(1, a.uid + N.uid.Nh, b.uid + N.Nh - N.uid.Nh)
    p.rnb.Nh = rbeta(1, a.rnb + N.rnb.Nh, b.rnb + N.Nh - N.rnb.Nh)
    
    
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








