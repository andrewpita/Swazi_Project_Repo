
msmoverlapdata=read.csv("msmoverlapdata.csv")
msmregionaldata=read.csv("msmregionaldata.csv")

library(dplyr)
library(BiasedUrn)
library(truncnorm)
library(coda)
library(stringr)
library(hash)

############
### data ###

loc_hash=hash(c("Lv","PP","Nh","ME","MM","Hh","Sw"),
              c("Lavumisa","Piggs Peak","Nhlangano",
                "Mbabane","Manzini","Hhohho","Shiselweni"))


## Lavumisa ##

#grab Lavumisa row from msmregionaldata, select the survey column, pull the data from
#that column (its just a scalar in this case)
#effectively its the number who took the survey

N.srv.Lv=msmregionaldata %>% filter(location=="Lavumisa") %>% select(survey) %>% pull

#same thing but getting the number who were given a user id

N.uid.Lv=msmregionaldata %>% filter(location=="Lavumisa") %>% select(uid) %>% pull

#again the number who were given a user id?
N.srv.uid.Lv=msmoverlapdata %>% filter(location=="Lavumisa",uid=="Yes") %>% summarise(n=sum(n)) %>% pull

#number who took the survey + number who were given the uid - the number who have both
#is 
r.Lv=N.srv.Lv+N.uid.Lv-N.srv.uid.Lv
P.Lv=msmregionaldata %>% filter(location=="Lavumisa") %>% select(mpop90pc_2014) %>% pull

## Pigg's Peak ##
N.srv.PP=msmregionaldata %>% filter(location=="Piggs Peak") %>% select(survey) %>% pull
N.uid.PP=msmregionaldata %>% filter(location=="Piggs Peak") %>% select(uid) %>% pull
N.srv.uid.PP=msmoverlapdata %>% filter(location=="Piggs Peak",uid=="Yes") %>% summarise(n=sum(n)) %>% pull
r.PP=N.srv.PP+N.uid.PP-N.srv.uid.PP
P.PP=msmregionaldata %>% filter(location=="Piggs Peak") %>% select(mpop90pc_2014) %>% pull

## Nhlangano ##
N.srv.Nh=msmregionaldata %>% filter(location=="Nhlangano") %>% select(survey) %>% pull
N.uid.Nh=msmregionaldata %>% filter(location=="Nhlangano") %>% select(uid) %>% pull
N.rnb.Nh=msmregionaldata %>% filter(location=="Nhlangano") %>% select(rainbow) %>% pull
N.srv.uid.Nh=msmoverlapdata %>% filter(location=="Nhlangano",uid=="Yes") %>% summarise(n=sum(n)) %>% pull
N.srv.rnb.Nh=msmoverlapdata %>% filter(location=="Nhlangano",rainbow=="Yes") %>% summarise(n=sum(n)) %>% pull
P.Nh=msmregionaldata %>% filter(location=="Nhlangano") %>% select(mpop90pc_2014) %>% pull

## Mbabane/Ezulwini ##
N.srv.ME=msmregionaldata %>% filter(location=="Mbabane/Ezulwini") %>% 
  select(survey) %>% pull
N.uid.ME=msmregionaldata %>% filter(location=="Mbabane/Ezulwini") %>% 
  select(uid) %>% pull
N.srv.uid.ME=msmoverlapdata %>% filter(location=="Mbabane/Ezulwini",uid=="Yes") %>% 
  summarise(n=sum(n)) %>% pull
N.srv.cpn.ME=msmoverlapdata %>% filter(location=="Mbabane/Ezulwini",coupon=="Yes") %>% 
  summarise(n=sum(n)) %>% pull
P.ME=msmregionaldata %>% filter(location=="Mbabane/Ezulwini") %>% select(mpop90pc_2014) %>% 
  pull

## Manzini/Matsapha ##
N.srv.MM=msmregionaldata %>% filter(location=="Manzini/Matsapha") %>% select(survey) %>% pull
N.uid.MM=msmregionaldata %>% filter(location=="Manzini/Matsapha") %>% select(uid) %>% pull
N.srv.uid.MM=msmoverlapdata %>% filter(location=="Manzini/Matsapha",uid=="Yes") %>% summarise(n=sum(n)) %>% pull
N.srv.cpn.MM=msmoverlapdata %>% filter(location=="Manzini/Matsapha",coupon=="Yes") %>% summarise(n=sum(n)) %>% pull
P.MM=msmregionaldata %>% filter(location=="Manzini/Matsapha") %>% select(mpop90pc_2014) %>% pull

## Corridor ##
N.cpn.Co=msmregionaldata %>% filter(location=="Corridor") %>% select(coupon) %>% pull

parnames=c("p.srv.Lv","p.uid.Lv","N.Lv",
           "p.srv.PP","p.uid.PP","N.PP",
           "p.srv.Nh","p.uid.Nh","p.rnb.Nh","r.Nh","N.Nh",
           "p.srv.ME","p.uid.ME","p.cpn.ME","N.cpn.ME","r.ME","N.ME",
           "p.srv.MM","p.uid.MM","p.cpn.MM","N.cpn.MM","r.MM","N.MM")

probpars=parnames[which(startsWith(parnames,"p"))]

### Gibbs sampler 
misaligned.caprecap=function(seed,N=10000,atune=0.25,btune=1){
    
    print(seed)
    M=matrix(0,N,length(parnames)) ## stores MCMC samples 
    colnames(M)=parnames
    
    ### constants ###
    a.srv=a.uid=a.rnb=a.cpn=1
    b.srv=b.uid=b.rnb=b.cpn=1
    
    ### initializing the chain ###
    for(par in probpars) assign(par,0.1)
    N.Nh=round(N.srv.Nh/p.srv.Nh)
    N.ME=round(N.srv.ME/p.srv.ME)
    N.MM=round(N.srv.MM/p.srv.MM)
    N.cpn.ME=max(N.srv.cpn.ME,round(N.cpn.Co*N.ME/(N.ME+N.MM)))
    N.cpn.MM=N.cpn.Co-N.cpn.ME
    #a.phi=b.phi=1.5
    
    set.seed(seed)     
    for(i in 1:N){    
        
        #### updating Lavumisa ####
        x=(1-p.srv.Lv)*(1-p.uid.Lv)
        N.Lv=r.Lv+rnbinom(1,r.Lv,1-x)
        p.srv.Lv=rbeta(1,a.srv+N.srv.Lv,b.srv+N.Lv-N.srv.Lv)
        p.uid.Lv=rbeta(1,a.uid+N.uid.Lv,b.uid+N.Lv-N.uid.Lv)
        #phi.Lv=rbeta(1,a.phi+N.Lv,b.phi+P.Lv-N.Lv)
        
        #### updating Piggs Peak ####
        x=(1-p.srv.PP)*(1-p.uid.PP)
        N.PP=r.PP+rnbinom(1,r.PP,1-x)
        p.srv.PP=rbeta(1,a.srv+N.srv.PP,b.srv+N.PP-N.srv.PP)
        p.uid.PP=rbeta(1,a.uid+N.uid.PP,b.uid+N.PP-N.uid.PP)
        #phi.PP=rbeta(1,a.phi+N.PP,b.phi+P.PP-N.PP)
        
        ### updating Nhlangano ####
        r.Nh=N.srv.Nh+N.uid.Nh+N.rnb.Nh-N.srv.uid.Nh-N.srv.rnb.Nh - 
            rhyper(1,N.uid.Nh-N.srv.uid.Nh,
                   N.Nh-N.srv.Nh-N.uid.Nh+N.srv.uid.Nh,
                   N.rnb.Nh-N.srv.rnb.Nh)
        x=(1-p.srv.Nh)*(1-p.uid.Nh)*(1-p.rnb.Nh)
        N.Nh=r.Nh+rnbinom(1,r.Nh,1-x)
        p.srv.Nh=rbeta(1,a.srv+N.srv.Nh,b.srv+N.Nh-N.srv.Nh)
        p.uid.Nh=rbeta(1,a.uid+N.uid.Nh,b.uid+N.Nh-N.uid.Nh)
        p.rnb.Nh=rbeta(1,a.rnb+N.rnb.Nh,b.rnb+N.Nh-N.rnb.Nh)
        #phi.Nh=rbeta(1,a.phi+N.Nh,b.phi+P.Nh-N.Nh)
        
        ### updating Mbabane/Ezulwini ####
        ov.ME=rhyper(1,N.uid.ME-N.srv.uid.ME,N.ME-N.srv.ME-N.uid.ME+N.srv.uid.ME,
                     N.cpn.ME-N.srv.cpn.ME)
        r.ME=N.srv.ME+N.uid.ME+N.cpn.ME-N.srv.uid.ME-N.srv.cpn.ME - ov.ME
        x=(1-p.srv.ME)*(1-p.uid.ME)*(1-p.cpn.ME)
        N.ME=r.ME+rnbinom(1,r.ME,1-x)
        p.srv.ME=rbeta(1,a.srv+N.srv.ME,b.srv+N.ME-N.srv.ME)
        p.uid.ME=rbeta(1,a.uid+N.uid.ME,b.uid+N.ME-N.uid.ME)
        p.cpn.ME=rbeta(1,a.cpn+N.cpn.ME,b.cpn+N.ME-N.cpn.ME)
        #phi.ME=rbeta(1,a.phi+N.ME,b.phi+P.ME-N.ME)
        
        #### updating Manzini/Matsapha ####
        ov.MM=rhyper(1,N.uid.MM-N.srv.uid.MM,N.MM-N.srv.MM-N.uid.MM+N.srv.uid.MM,N.cpn.MM-N.srv.cpn.MM)
        r.MM=N.srv.MM+N.uid.MM+N.cpn.MM-N.srv.uid.MM-N.srv.cpn.MM - ov.MM
        x=(1-p.srv.MM)*(1-p.uid.MM)*(1-p.cpn.MM)
        N.MM=r.MM+rnbinom(1,r.MM,1-x)
        p.srv.MM=rbeta(1,a.srv+N.srv.MM,b.srv+N.MM-N.srv.MM)
        p.uid.MM=rbeta(1,a.uid+N.uid.MM,b.uid+N.MM-N.uid.MM)
        p.cpn.MM=rbeta(1,a.cpn+N.cpn.MM,b.cpn+N.MM-N.cpn.MM)
        #phi.MM=rbeta(1,a.phi+N.MM,b.phi+P.MM-N.MM)
        
        #### updating coupon distribution in the corridor ####
        N.cpn.ME=N.srv.cpn.ME+ov.ME+rFNCHypergeo(1,N.ME-N.srv.ME-N.uid.ME+N.srv.uid.ME,
            N.MM-N.srv.MM-N.uid.MM+N.srv.uid.MM,
            N.cpn.Co-N.srv.cpn.ME-ov.ME-N.srv.cpn.MM-ov.MM,
            (p.cpn.ME/(1-p.cpn.ME))/(p.cpn.MM/(1-p.cpn.MM)))
        N.cpn.MM=N.cpn.Co-N.cpn.ME
        
        for(par in parnames) M[i,par]=get(par)
        
        if((i %% 1)==200) print(i)    
        #print(i)
    }    
    M1=cbind(M[,"N.ME"]+M[,"N.PP"],M[,"N.Nh"]+M[,"N.Lv"])
    colnames(M1)=c("N.Hh","N.Sw")
    cbind(M,M1)
}

### trying out different seeds ###
N=10000
seeds=.Random.seed[1:10]
Mlist=lapply(seeds,misaligned.caprecap,N)

#M1=misaligned.caprecap(1,N)
#M12345=misaligned.caprecap(12345,N)

### gelman rubin diagnostics for all independent variables (removing N.MM and the method hyper-params)
MCMClist=lapply(Mlist,function(x) mcmc(x[,setdiff(colnames(x),c("N.MM",tail(colnames(x),-9)))]))
grdiag=gelman.diag(MCMClist)
print(grdiag)

mytraceplot=function(col,Mlist){
    colname=str_replace_all(col,"[.]","_")
    imagename=paste0("results/traceplots/",colname,".png")
    png(imagename)
    plot(Mlist[[1]][,col],type="l",xlab="Iterations",ylab=col)
    lines(Mlist[[2]][,col],col="red")
    lines(Mlist[[3]][,col],col="green")
    dev.off()
}

mydensityplot=function(col,Mlist){
    Nburn=nrow(Mlist[[1]])/2
    colname=str_replace_all(col,"[.]","_")
    imagename=paste0("results/densityplots/",colname,".png")
    png(imagename)
    plot(density(Mlist[[1]][-(1:Nburn),col]),type="l",
         main=loc_hash[[str_sub(col,-2)]],xlab="x",ylab="")
    lines(density(Mlist[[2]][-(1:Nburn),col]),col="red")
    lines(density(Mlist[[3]][-(1:Nburn),col]),col="green")
    dev.off()
}

myhistograms=function(col,M){
    Nburn=nrow(Mlist[[1]])/2
    colname=str_replace_all(col,"[.]","_")
    imagename=paste0("results/histograms/",colname,".png")
    png(imagename)
    plot(hist(Mlist[[1]][-(1:Nburn),col]),
         main=loc_hash[[str_sub(col,-2)]],xlab="value",ylab="count")
    dev.off()
}

mysummary=function(x){ 
    summ=c(min(x),quantile(x,c(0.01,0.025,0.05,0.1,0.2,0.25)),median(x),mean(x),sd(x),quantile(x,1-rev(c(0.01,0.025,0.05,0.1,0.2,0.25))),max(x))
    names(summ)=c("min","1%","2.5%","5%","10%","20%","25%","Median","Mean","sd","75%","80%","90%","95%","97.5%","99%","max")
    summ
    }

sapply(colnames(Mlist[[1]]),mytraceplot,Mlist)
sapply(colnames(Mlist[[1]]),mydensityplot,Mlist)

popsizecols=paste0("N.",names(loc_hash))
sapply(popsizecols,myhistograms,Mlist[[1]])

r.vec = c(r.Lv, r.PP, N.srv.uid.Nh, N.srv.rnb.Nh, N.srv.uid.ME,
          N.srv.cpn.ME, N.srv.uid.MM, N.srv.cpn.MM,N.cpn.Co)
names(r.vec) = c("r.Lv", "r.PP", "N.srv.uid.Nh", "N.srv.rnb.Nh", 
                 "N.srv.uid.ME","N.srv.cpn.ME", "N.srv.uid.MM",
                 "N.srv.cpn.MM", "N.cpn.Co")
r.vec = as.matrix(r.vec)
write.csv(r.vec,"results/r.csv")

sumtab=t(apply(Mlist[[1]],2,mysummary))
write.csv(sumtab,"results/summary.csv")

save.image(file="results/results.Rdata")
