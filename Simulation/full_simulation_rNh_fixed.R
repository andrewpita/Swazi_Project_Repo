library(dplyr)
library(BiasedUrn)
library(truncnorm)
library(coda)
library(stringr)
library(hash)
library(gdata)



#tie the city abbreviations to their full names

loc_hash=hash(c("Lv","PP","Nh","ME","MM","Hh","Sw"),
              c("Lavumisa","Piggs Peak","Nhlangano",
                "Mbabane","Manzini","Hhohho","Shiselweni"))


#read in output from Gibbs sampler to use as true values
#in simulation
results = read.csv("summary.csv", row.names = 1)
nmes = rownames(results)

#we use the posterior means as the true values
results = results[,"Mean"]
names(results) = nmes

#a function that takes as argument a string that is 
#the two letter abbreviation for each city.
#It returns a subset that contains the "true"
#values for each cities parameters
make_city = function(string) {
  
  #grepl returns a vector of logical values, true
  #when a parameter name matches the string argument.
  #The $ means that we want a match at the end 
  indeces = grepl(paste(string,"$", sep = ""), names(results))
  
  city = results[indeces]
  names(city)=NULL
  return(city)
}

####first things first, grab the true parameters
#subset so we have "true" parameters for each city
Lv = make_city("Lv")
PP = make_city("PP")
Nh = make_city("Nh")
ME = make_city("ME")
MM = make_city("MM")

#grab true parameters for Lavumisa
true.N.Lv = round(Lv[3])

true.p.srv.Lv = Lv[1]

true.p.uid.Lv = Lv[2]

#grab true parameters for Piggs Peak

true.N.PP = round(PP[3])

true.p.srv.PP = PP[1]
true.p.uid.PP = PP[2]


#Nhlangano true parameters

true.N.Nh = round(Nh[5])

true.p.srv.Nh = Nh[1]
true.p.uid.Nh = Nh[2]
true.p.rnb.Nh = Nh[3]
#saving true value this time
true.r.Nh = round(Nh[4])



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

parnames=c("p.srv.Lv","p.uid.Lv","N.Lv",
           "p.srv.PP","p.uid.PP","N.PP",
           "p.srv.Nh","p.uid.Nh","p.rnb.Nh","r.Nh","N.Nh",
           "p.srv.ME","p.uid.ME","p.cpn.ME","N.cpn.ME","r.ME","N.ME",
           "p.srv.MM","p.uid.MM","p.cpn.MM","N.cpn.MM","r.MM","N.MM")


#index of parameters that are probabilities 
#(start with p)

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
  
  ### simulate data for Lavumisa ####
  
  Lv.sample = rmultinom(1, size = true.N.Lv, prob = 
                          c(true.p.srv.Lv*(1 - true.p.uid.Lv), #in survey not in uid
                            true.p.srv.Lv*true.p.uid.Lv, #in both survey and uid
                            true.p.uid.Lv*(1 - true.p.srv.Lv), #in uid and not in survey
                            (1- true.p.srv.Lv)*(1 - true.p.uid.Lv))) #not in either 
  
  
  #store simulated data for use in Gibbs sampler
  N.srv.uid.Lv = Lv.sample[2]
  N.srv.Lv = Lv.sample[1] + N.srv.uid.Lv
  N.uid.Lv = Lv.sample[3] + N.srv.uid.Lv
  true.r.Lv = r.Lv = N.srv.Lv + N.uid.Lv - N.srv.uid.Lv
  
  ### Simulate data for Pigg's Peak ###
  
  PP.sample = rmultinom(1, size = true.N.PP, prob = 
                          c(true.p.srv.PP*(1 - true.p.uid.PP), #in survey not in uid
                            true.p.srv.PP*true.p.uid.PP, # in both
                            true.p.uid.PP*(1 - true.p.srv.PP),  #in uid not in survey
                            (1- true.p.srv.PP)*(1 - true.p.uid.PP))) #not in either
  
  
  #store simulated data
  N.srv.uid.PP = PP.sample[2]
  N.srv.PP = PP.sample[1] + N.srv.uid.PP
  N.uid.PP = PP.sample[3] + N.srv.uid.PP
  true.r.PP = r.PP = N.srv.PP + N.uid.PP - N.srv.uid.PP
  
  ### Simulate Nhlangano data ###
  
  Nh.sample = rmultinom(1, size = true.r.Nh, prob = 
                          c(true.p.srv.Nh*(1 - true.p.uid.Nh)*(1 - true.p.rnb.Nh), #only in survey
                            true.p.srv.Nh*true.p.uid.Nh*(1- true.p.rnb.Nh), #in only survey and uid
                            true.p.srv.Nh*true.p.uid.Nh * true.p.rnb.Nh, #in all three
                            true.p.uid.Nh*(1 - true.p.srv.Nh)*(1 - true.p.rnb.Nh), #in only uid
                            true.p.uid.Nh*(1 - true.p.srv.Nh)*true.p.rnb.Nh, #in uid and rainbow
                            true.p.rnb.Nh*(1 - true.p.srv.Nh)*(1 - true.p.uid.Nh), #in only rnb
                            true.p.rnb.Nh*true.p.srv.Nh*(1 - true.p.uid.Nh)#in rnb and srv
                            #(1- true.p.srv.Nh)*(1 - true.p.uid.N)* (1 - true.p.rnb.Nh)#in none
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
  #true.r.Nh = as.numeric(true.N.Nh - Nh.sample[8] )
  N.in.None = true.N.Nh - true.r.Nh
  
  
  
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
  
  #initializing the chain with start values for
  #the probability full conditionals
  
  for (par in probpars) {
    
    #par = 0.1
    assign(par,0.1)
  }
  
  #create initial estimates of the total 
  #MSM populations by dividing the number 
  #we observed in survey by the intial value for 
  #the probabilities
  
  N.Nh = round(N.srv.Nh/p.srv.Nh)
  N.ME = round(N.srv.ME/p.srv.ME)
  N.MM = round(N.srv.MM/p.srv.MM)
  N.cpn.ME = max(N.srv.cpn.ME, round(N.cpn.Co*N.ME/(N.ME + N.MM)))
  N.cpn.MM = N.cpn.Co - N.cpn.ME
  
  #begin Gibbs sampler algorithm
  set.seed(seed)
  
  for (i in 1:N) {
    
    ### updating Lavumisa
    
    #x is the probability of not being in either survey
    #we use this as a parameter in the negative binomial 
    #posterior for N.Lv, the total population of MSM in
    #Lavumisa
    
    x = ((1 - p.srv.Lv) * (1 - p.uid.Lv))
    N.Lv = r.Lv + rnbinom(1, r.Lv, 1 - x)
    
    #draw from posterior distribution for the probabilities
    p.srv.Lv = rbeta(1, a.srv + N.srv.Lv, b.srv + N.Lv - N.srv.Lv)
    p.uid.Lv = rbeta(1, a.uid + N.uid.Lv, b.uid + N.Lv - N.uid.Lv)  
    
    ### updating Piggs Peak
    
    x = (1- p.srv.PP) * (1 - p.uid.PP)
    
    N.PP = r.PP + rnbinom(1, r.PP, 1- x)
    
    p.srv.PP = rbeta(1, a.srv + N.srv.PP, b.srv + N.PP - N.srv.PP)
    p.uid.PP = rbeta(1, a.uid + N.uid.PP, b.uid + N.PP - N.uid.PP)
    
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
    
    #assign parameter values from this iteration to M
    #matrix
    
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

write.csv(avgsum,"1summary.csv")

