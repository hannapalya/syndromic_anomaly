rm(list=ls(all=TRUE))

# Disable plotting to prevent script from hanging
pdf(NULL)  # Redirect plots to null device

# FUNCTIONS THAT PRODUCE THE DATA
# DEFINING FUNCTION h

#==============
# 5-day systems
#==============

h1=function(N,k,k2,alpha,beta,gama1,gama2,gama3,gama4,shift,shift2){
  t=1:N
  if(k==0 & k2==0){h1=alpha+beta*t}
  else{
    if(k==0)
    {
      l=1:k2
      h1=rep(0,N)
      for(i in 1:N){
        h1[i]=alpha+beta*(t[i]+shift)+sum(gama3*cos((2*pi*l*(t[i]+shift))/5)+gama4*sin((2*pi*l*(t[i]+shift))/5))
      }
    }
    else{
      j=1:k
      l=1:k2
      h1=rep(0,N)
      for(i in 1:N){
        h1[i]=alpha+beta*(t[i]+shift)+sum(gama1*cos((2*pi*j*(t[i]+shift))/(52*5))+gama2*sin((2*pi*j*(t[i]+shift2))/(52*5)))+sum(gama3*cos((2*pi*l*(t[i]+shift))/5)+gama4*sin((2*pi*l*(t[i]+shift))/5))
      }
    }
  }
  h1
}

negbinNoise1=function(N,k,k2,alpha,beta,gama1,gama2,gama3,gama4,phi,shift,shift2){
  mu <- exp(h1(N,k,k2,alpha,beta,gama1,gama2,gama3,gama4,shift,shift2))
  if(phi==1){yi <- rpois(N,mu)}
  else{
    prob <- 1/phi
    size <- mu/(phi-1)
    yi <- rnbinom(N,size=size,prob=prob)
  }
  yi
}


outbreak5=function(currentday,weeklength,wtime,yi,interval,k,k2,alpha,beta,gama1,gama2,gama3,gama4,shift,shift2,phi,numoutbk,peakoutbk,meanlog,sdlog){
  # theta, beta, gama1 and gama2 are the parameters of the equation for mu in Section 3.1

  N=length(yi)
  t=1:N

  mu <- exp(h1(N,k,k2,alpha,beta,gama1,gama2,gama3,gama4,shift,shift2))
  s=sqrt(mu*phi)

  #wtime = (currentday-49*5+1):currentday         # current outbreaks

  # GENERATING OUTBREAKS

  # STARTING TIMES OF OUTBREAKS

  startoutbk <- sample(wtime, numoutbk, replace = FALSE)

  # OUTBREAK SIZE OF CASES

  sizeoutbk=rep(0,numoutbk)
  for(i in 1:numoutbk){
    set.seed(i)
    soutbk=1
    sou=1
    while(soutbk<2){
      set.seed(sou)
      soutbk=rpois(1,s[startoutbk[i]]*peakoutbk)
      sou=sou+1
    }
    sizeoutbk[i]=soutbk
  }

  # DISTRIBUTE THESE CASES OVER TIME USING LOGNORMAL

  outbreak=rep(0,2*N)
  for( j in 1:numoutbk){
    set.seed(j)
    outbk <-rlnorm(sizeoutbk[j], meanlog = meanlog, sdlog = sdlog)
    #outbk <-rnorm(sizeoutbk[j], mean = meanlog2, sd = sdlog)
    #h<- hist(outbk,breaks=seq(0,ceiling(max(outbk)),1),plot=FALSE)
    h<- hist(outbk,breaks=seq(0,ceiling(max(outbk)),interval),plot=FALSE)
    cases <- h$counts
    weight=rep(0,length(cases))
    duration<-startoutbk:(startoutbk+length(cases)-1)
    dayofweek<-duration%%5 # 0 is friday; 1 is monday; 2 is tuesday etc.
    for(i in 1:length(cases)){
      if(dayofweek[i]==0){weight[i]=1.1}
      if(dayofweek[i]==1){weight[i]=1.5}
      if(dayofweek[i]==2){weight[i]=1.1}
      if(dayofweek[i]==3){weight[i]=1}
      if(dayofweek[i]==4){weight[i]=1}
    }
    cases2 <- cases*weight
    for (l in 1:(length(cases2))){
      outbreak[startoutbk[j]+(l-1)]= cases2[l]+outbreak[startoutbk[j]+(l-1)]
    }# l loop
  }# j loop

  #for(v in 1:(currentday-49*5)){if(outbreak[v]>0){outbreak[v]=0}}
  for(v in currentday:(currentday+100)){if(outbreak[v]>0){outbreak[v]=0}}
  outbreak=outbreak[1:N]

  # ADD NOISE AND OUTBREAKS

  yitot=yi+outbreak
  result=list(yitot=yitot,outbreak=outbreak,startoutbk=startoutbk,sizeoutbk=sizeoutbk,sd=s,mean=mu)
  #return(result)
}

#==============
# 7-day systems
#==============

h2=function(N,k,k2,alpha,beta,gama1,gama2,gama3,gama4,shift){
  t=1:N
  if(k==0 & k2==0){h2=alpha+beta*t}
  else{
    if(k==0)
    {
      l=1:k2
      h2=rep(0,N)
      for(i in 1:N){
        h2[i]=alpha+beta*(t[i]+shift)+sum(gama3*cos((2*pi*l*(t[i]+shift))/7)+gama4*sin((2*pi*l*(t[i]+shift))/7))
      }
    }
    else{
      j=1:k
      l=1:k2
      h2=rep(0,N)
      for(i in 1:N){
        h2[i]=alpha+beta*(t[i]+shift)+sum(gama1*cos((2*pi*j*(t[i]+shift))/(52*7))+gama2*sin((2*pi*j*(t[i]+shift))/(52*7)))+sum(gama3*cos((2*pi*l*(t[i]+shift))/7)+gama4*sin((2*pi*l*(t[i]+shift))/7))
      }
    }
  }
  h2
}


negbinNoise2=function(N,k,k2,alpha,beta,gama1,gama2,gama3,gama4,phi,shift){
  mu <- exp(h2(N,k,k2,alpha,beta,gama1,gama2,gama3,gama4,shift))
  if(phi==1){yi <- rpois(N,mu)}
  else{
    prob <- 1/phi
    size <- mu/(phi-1)
    yi <- rnbinom(N,size=size,prob=prob)
  }
  yi
}


outbreak7=function(currentday,weeklength,wtime,yi,interval,k,k2,alpha,beta,gama1,gama2,gama3,gama4,shift,phi,numoutbk,peakoutbk,meanlog,sdlog){
  # theta, beta, gama1 and gama2 are the parameters of the equation for mu in Section 3.1

  N=length(yi)
  t=1:N

  mu <- exp(h2(N,k,k2,alpha,beta,gama1,gama2,gama3,gama4,shift))
  s=sqrt(mu*phi)


  #wtime = (currentday-49*7+1):currentday         # current outbreaks
  # wtime = 350*1:7         # current outbreaks

  # GENERATING OUTBREAKS

  # STARTING TIMES OF OUTBREAKS

  startoutbk <- sample(wtime, numoutbk, replace = FALSE)

  # OUTBREAK SIZE OF CASES

  sizeoutbk=rep(0,numoutbk)
  for(i in 1:numoutbk){
    set.seed(i)
    soutbk=1
    sou=1
    while(soutbk<2){
      set.seed(sou)
      soutbk=rpois(1,s[startoutbk[i]]*peakoutbk)
      sou=sou+1
    }
    sizeoutbk[i]=soutbk
  }

  # DISTRIBUTE THESE CASES OVER TIME USING LOGNORMAL

  outbreak=rep(0,2*N)
  for( j in 1:numoutbk){
    set.seed(j)
    outbk <-rlnorm(sizeoutbk[j], meanlog = meanlog, sdlog = sdlog)
    #outbk <-rnorm(sizeoutbk[j], mean = meanlog2, sd = sdlog)
    #h<- hist(outbk,breaks=seq(0,ceiling(max(outbk)),1),plot=FALSE)
    h<- hist(outbk,breaks=seq(0,ceiling(max(outbk)),interval),plot=FALSE)
    cases <- h$counts
    weight=rep(0,length(cases))
    duration<-startoutbk:(startoutbk+length(cases)-1)
    dayofweek<-duration%%7 # 0 is sunday; 1 is monday; 2 is tuesday etc.
    for(i in 1:length(cases)){
      if(dayofweek[i]==0){weight[i]=2}
      if(dayofweek[i]==1){weight[i]=1}
      if(dayofweek[i]==2){weight[i]=1}
      if(dayofweek[i]==3){weight[i]=1}
      if(dayofweek[i]==4){weight[i]=1}
      if(dayofweek[i]==5){weight[i]=1}
      if(dayofweek[i]==6){weight[i]=2}
    }
    cases2 <- cases*weight
    for (l in 1:(length(cases2))){
      outbreak[startoutbk[j]+(l-1)]= cases2[l]+outbreak[startoutbk[j]+(l-1)]
    }# l loop
  }# j loop

  #for(v in (currentday-49*7):currentday){if(outbreak[v]>0){outbreak[v]=0}}
  for(v in currentday:(currentday+100)){if(outbreak[v]>0){outbreak[v]=0}}
  outbreak=outbreak[1:N]

  # ADD NOISE AND OUTBREAKS

  yitot=yi+outbreak
  result=list(yitot=yitot,outbreak=outbreak,startoutbk=startoutbk,sizeoutbk=sizeoutbk,sd=s,mean=mu)
  #return(result)
}

#==========================
# Specify the bank holidays
#==========================

years=7
bankholidays=read.csv("Bankholidays.csv",sep=',')

#fix(bankholidays)
bankhols7=bankholidays$bankhol
bankhols7=as.numeric(bankhols7)
length(bankhols7)
#fix(bankhols7)

bankhols5=bankhols7[-seq(6,length(bankhols7),7)]
bankhols5=bankhols5[-seq(6,length(bankhols5),6)]
bankhols5=as.numeric(bankhols5)
length(bankhols5)
#fix(bankhols5)

#=======================
# Define the data frames
#=======================

nsim=100

simulateddata1=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata2=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata3=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata4=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata5=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata6=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata7=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata8=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata9=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata10=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata11=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata12=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata13=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata14=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata15=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata16=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))


simulatedtotals1=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals2=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals3=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals4=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals5=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals6=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals7=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals8=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals9=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals10=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals11=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals12=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals13=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals14=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals15=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals16=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))

simulatedoutbreak1=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak2=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak3=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak4=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak5=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak6=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak7=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak8=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak9=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak10=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak11=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak12=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak13=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak14=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak15=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak16=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))


simulatedzseasoutbreak5=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedzseasoutbreak6=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedzseasoutbreak15=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))

#################################
#SIMULATE SYNDROMES AND OUTBREAKS
#################################

########################################
# SIMULATIONS WITH VERY SMALL(w=2), SMALL(w=3), MEDIUM(w=5) AND LARGE(w=10) OUTBREAKS
########################################

# VERY SMALL w=2

w=5

#=====================
# 5-day week syndromes
#=====================


days5=5
N=52*days5*years

#sigid5
print("Starting sigid5 simulations...")
for(i in 1:nsim){
  print(paste("sigid5: Simulation", i, "of", nsim))
  set.seed(i)
  yt=negbinNoise1(N=N,k=1,k2=1,alpha=6,beta=0,gama1=0.3,gama2=2,gama3=0.3,gama4=0.5,phi=1.5,shift=-50,shift2=-50)/10
  #mu=exp(h1(N=N,k=1,k2=1,alpha=6,beta=0,gama1=0.3,gama2=2,gama3=0.3,gama4=0.5,shift=-50,shift2=-50))

  out1=rep(0,N)
  for(j in 1:years){
    print(paste("  sigid5: Year", j, "of", years))
    set.seed(j+years*i)
    out=outbreak5(currentday=days5*52*years,weeklength=52*days5*years,wtime=((1+(j-1)*days5*52):(20+(j-1)*days5*52)),yi=yt,interval=0.02,k=1,
                  k2=1,alpha=6,beta=0,gama1=0.3,gama2=2,gama3=0.3,gama4=0.5,phi=1.5,shift=-50,shift2=-50,numoutbk=1,peakoutbk=3*days5*80,meanlog=0,sdlog=0.5)
    out1=out1+out$outbreak
  }

  set.seed(i)
  out2=outbreak5(currentday=days5*52*years,weeklength=52*days5*years,wtime=(length(yt)-49*days5+1):length(yt),yi=yt,interval=0.25,k=1,k2=1,alpha=6,beta=0,gama1=0.3,
                 gama2=2,gama3=0.3,gama4=0.5,phi=1.5,shift=-50,shift2=-50,numoutbk=1,peakoutbk=w*days5,meanlog=0,sdlog=0.5)

  zoutbreak=out2$outbreak
  zseasoutbreak=out2$outbreak +out1

  zt=yt +out1
  zitot=yt + out2$outbreak +out1

  #zitot[(bankhols5==1)]=0
  #zitot[(bankhols5==1)+1]=1.5*zitot[i+1]

  for(b in 1:length(zitot)){
    if(bankhols5[b]==1){
      zitot[b]=0
      zitot[b+1]=1.5*zitot[b+1]
    }
  }

  zeros=rep(0,2)
  weekend=seq(days5,days5*years*52,days5)
  #weekend=seq(0,days5*years*52-1,days5)
  for(s in 1:length(weekend)){
    yt=append(yt,zeros,after=2*(s-1)+weekend[s])
    zt=append(zt,zeros,after=2*(s-1)+weekend[s])
    zitot=append(zitot,zeros,after=2*(s-1)+weekend[s])
    zoutbreak=append(zoutbreak,zeros,after=2*(s-1)+weekend[s])
    zseasoutbreak=append(zseasoutbreak,zeros,after=2*(s-1)+weekend[s])
  }

  simulateddata5[,i]=round(zt)
  simulatedtotals5[,i]=round(zitot)
  simulatedoutbreak5[,i]=round(zoutbreak)
  simulatedzseasoutbreak5[,i]=round(zseasoutbreak)
}

#----------------------------------------------------
# Plot the datasets and outbreaks using the following
#----------------------------------------------------

#plot(1:N,yt,typ='l')
#plot(1:(52*years*7),zt,typ='l',xlim=c(1,364))
#plot(1:(52*years*7),zitot,typ='l',xlim=c(1,364))
#plot(1:(52*years*7),zitot,typ='l',xlim=c(365,728))
#plot(1:(52*years*7),zitot,typ='l',xlim=c(729,1092))
#plot(1:(52*years*7),zitot,typ='l',xlim=c(1093,1456))
#plot(1:(52*years*7),zitot,typ='l',xlim=c(1457,1820))
#plot(1:(52*years*7),zitot,typ='l',xlim=c(1821,2184))
#plot(1:(52*years*7),zitot,typ='l',xlim=c(2185,2548))
#lines(1:(52*years*7),zoutbreak,col='green')

plot(1:(52*years*7),simulatedtotals5[,1],typ='l',xlim=c(1,7*364))
lines(1:(52*years*7),simulatedzseasoutbreak5[,1],col='green')
lines(1:(52*years*7),simulatedoutbreak5[,1],col='red')

plot(1:(52*years*7),zitot,typ='l',xlim=c(1,7*364))
lines(1:(52*years*7),yt,col='blue')
lines(1:(52*years*7),zseasoutbreak,col='green')
lines(1:(52*years*7),zoutbreak,col='red')

plot(1:(52*years*7),simulatedzseasoutbreak5[,1],col='green',typ='l')

#sigid6
print("Starting sigid6 simulations...")
for(i in 1:nsim){
  print(paste("sigid6: Simulation", i, "of", nsim))
  set.seed(i)
  yt=negbinNoise1(N=N,k=1,k2=1,alpha=1,beta=0,gama1=0.1,gama2=2,gama3=0.05,gama4=0.05,phi=1,shift=-50,shift2=-50)
  #mu=exp(h1(N=N,k=1,k2=1,alpha=1.5,beta=0,gama1=0.1,gama2=2,gama3=0.1,gama4=0.1,shift=-50,shift2=-50))

  out1=rep(0,N)
  for(j in 1:years){
    set.seed(j+years*i)
    out=outbreak5(currentday=days5*52*years,weeklength=52*days5*years,wtime=((1+(j-1)*days5*52):(20+(j-1)*days5*52)),yi=yt,interval=0.02,k=1,k2=1,alpha=1,beta=0,gama1=0.1,
                  gama2=2,gama3=0.05,gama4=0.05,phi=1,shift=-50,shift2=-50,numoutbk=1,peakoutbk=3*days5*50,meanlog=0,sdlog=0.5)
    out1=out1+out$outbreak
  }

  set.seed(i)
  out2=outbreak5(currentday=days5*52*years,weeklength=52*days5*years,wtime=(length(yt)-49*days5+1):length(yt),yi=yt,interval=0.25,k=1,k2=1,alpha=1,beta=0,gama1=0.1,
                 gama2=2,gama3=0.05,gama4=0.05,phi=1,shift=-50,shift2=-50,numoutbk=1,peakoutbk=w*days5,meanlog=0,sdlog=0.5)

  zoutbreak=out2$outbreak
  zseasoutbreak=out2$outbreak+out1

  zt=yt +out1
  zitot=yt + out2$outbreak  +out1

  for(b in 1:length(zitot)){
    if(bankhols5[b]==1){
      zitot[b]=0
      zitot[b+1]=1.5*zitot[b+1]
    }
  }

  zeros=rep(0,2)
  weekend=seq(days5,days5*years*52,days5)
  #weekend=seq(0,days5*years*52-1,days5)
  for(s in 1:length(weekend)){
    yt=append(yt,zeros,after=2*(s-1)+weekend[s])
    zt=append(zt,zeros,after=2*(s-1)+weekend[s])
    zitot=append(zitot,zeros,after=2*(s-1)+weekend[s])
    zoutbreak=append(zoutbreak,zeros,after=2*(s-1)+weekend[s])
    zseasoutbreak=append(zseasoutbreak,zeros,after=2*(s-1)+weekend[s])
  }

  simulateddata6[,i]=round(zt)
  simulatedtotals6[,i]=round(zitot)
  simulatedoutbreak6[,i]=round(zoutbreak)
  simulatedzseasoutbreak6[,i]=round(zseasoutbreak)
}

plot(1:(52*years*7),simulatedtotals6[,1],typ='l',xlim=c(1,7*364))
lines(1:(52*years*7),simulatedzseasoutbreak6[,1],col='green')
lines(1:(52*years*7),simulatedoutbreak6[,1],col='red')

plot(1:(52*years*7),zitot,typ='l',xlim=c(1,7*364))
lines(1:(52*years*7),yt,col='blue')
lines(1:(52*years*7),zseasoutbreak,col='green')
lines(1:(52*years*7),zoutbreak,col='red')

plot(1:(52*years*7),simulatedzseasoutbreak6[,1],col='green',typ='l')


#sigid7
for(i in 1:nsim){
  set.seed(i)
  yt=negbinNoise1(N=N,k=0,k2=1,alpha=6,beta=0.0001,gama1=0,gama2=0,gama3=0.6,gama4=0.9,phi=1.5,shift=0,shift2=0)/10
  #mu=exp(h1(N=N,k=0,k2=1,alpha=6,beta=0,gama1=0,gama2=0,gama3=0.6,gama4=0.9,shift=0,shift2=0))

  set.seed(i)
  out2=outbreak5(currentday=days5*52*years,weeklength=52*days5*years,wtime=(length(yt)-49*days5+1):length(yt),yi=yt,interval=0.25,k=0,k2=1,alpha=6,beta=0,gama1=0,
                 gama2=0,gama3=0.6,gama4=0.9,phi=1.5,shift=0,shift2=0,numoutbk=1,peakoutbk=w*days5,meanlog=0,sdlog=0.5)

  zoutbreak=out2$outbreak

  zt=yt#+out1
  zitot=yt + out2$outbreak #+out1

  #zitot[(bankhols5==1)]=0
  #zitot[(bankhols5==1)+1]=1.5*zitot[i+1]

  for(b in 1:length(zitot)){
    if(bankhols5[b]==1){
      zitot[b]=0
      zitot[b+1]=1.5*zitot[b+1]
    }
  }

  zeros=rep(0,2)
  weekend=seq(days5,days5*years*52,days5)
  #weekend=seq(0,days5*years*52-1,days5)
  for(s in 1:length(weekend)){
    zt=append(zt,zeros,after=2*(s-1)+weekend[s])
    zitot=append(zitot,zeros,after=2*(s-1)+weekend[s])
    zoutbreak=append(zoutbreak,zeros,after=2*(s-1)+weekend[s])
  }

  simulateddata7[,i]=round(zt)
  simulatedtotals7[,i]=round(zitot)
  simulatedoutbreak7[,i]=round(zoutbreak)
}


plot(1:(52*years*7),simulateddata7[,1],typ='l',xlim=c(2185,2548),col='green')
lines(1:(52*years*7),simulateddata7[,1]+simulatedoutbreak7[,1])



#sigid8
for(i in 1:nsim){
  set.seed(i)
  yt=negbinNoise1(N=N,k=1,k2=1,alpha=3,beta=0,gama1=1.5,gama2=0.1,gama3=0.2,gama4=0.3,phi=1,shift=-150,shift2=-150)
  mu=exp(h1(N=N,k=1,k2=1,alpha=3,beta=0,gama1=1.5,gama2=0.1,gama3=0.6,gama4=0.8,shift=-150,shift2=-150))

  set.seed(i)
  out2=outbreak5(currentday=days5*52*years,weeklength=52*days5*years,wtime=(length(yt)-49*days5+1):length(yt),interval=0.25,yi=yt,k=1,k2=1,alpha=3,beta=0,gama1=1.5,
                 gama2=0.1,gama3=0.2,gama4=0.3,phi=1,shift=-150,shift2=-150,numoutbk=1,peakoutbk=w*days5,meanlog=0,sdlog=0.5)

  zoutbreak=out2$outbreak

  zt=yt#+out1
  zitot=yt + out2$outbreak #+out1

  #zitot[(bankhols5==1)]=0
  #zitot[(bankhols5==1)+1]=1.5*zitot[i+1]

  for(b in 1:length(zitot)){
    if(bankhols5[b]==1){
      zitot[b]=0
      zitot[b+1]=1.5*zitot[b+1]
    }
  }

  zeros=rep(0,2)
  weekend=seq(days5,days5*years*52,days5)
  #weekend=seq(0,days5*years*52-1,days5)
  for(s in 1:length(weekend)){
    zt=append(zt,zeros,after=2*(s-1)+weekend[s])
    zitot=append(zitot,zeros,after=2*(s-1)+weekend[s])
    zoutbreak=append(zoutbreak,zeros,after=2*(s-1)+weekend[s])
  }

  simulateddata8[,i]=round(zt)
  simulatedtotals8[,i]=round(zitot)
  simulatedoutbreak8[,i]=round(zoutbreak)
}


plot(1:(52*years*7),simulateddata8[,1],typ='l',xlim=c(2185,2548),col='green')
lines(1:(52*years*7),simulateddata8[,1]+simulatedoutbreak8[,1])


#sigid9
for(i in 1:nsim){
  set.seed(i)
  yt=negbinNoise1(N=N,k=1,k2=1,alpha=3,beta=0,gama1=0.2,gama2=0.1,gama3=0.05,gama4=0.15,phi=1,shift=-200,shift2=-200)
  #mu=exp(h1(N=N,k=1,k2=1,alpha=3,beta=0,gama1=0.2,gama2=0.1,gama3=0.05,gama4=0.15,shift=-200,shift2=-200))

  set.seed(i)
  out2=outbreak5(currentday=days5*52*years,weeklength=52*days5*years,wtime=(length(yt)-49*days5+1):length(yt),yi=yt,interval=0.25,k=1,k2=1,alpha=3,beta=0,gama1=0.2,
                 gama2=0.1,gama3=0.05,gama4=0.15,phi=1,shift=-200,shift2=-200,numoutbk=1,peakoutbk=w*days5,meanlog=0,sdlog=0.5)

  zoutbreak=out2$outbreak

  zt=yt#+out1
  zitot=yt + out2$outbreak #+out1

  #zitot[(bankhols5==1)]=0
  #zitot[(bankhols5==1)+1]=1.5*zitot[i+1]

  for(b in 1:length(zitot)){
    if(bankhols5[b]==1){
      zitot[b]=0
      zitot[b+1]=1.5*zitot[b+1]
    }
  }

  zeros=rep(0,2)
  weekend=seq(days5,days5*years*52,days5)
  #weekend=seq(0,days5*years*52-1,days5)
  for(s in 1:length(weekend)){
    zt=append(zt,zeros,after=2*(s-1)+weekend[s])
    zitot=append(zitot,zeros,after=2*(s-1)+weekend[s])
    zoutbreak=append(zoutbreak,zeros,after=2*(s-1)+weekend[s])
  }

  simulateddata9[,i]=round(zt)
  simulatedtotals9[,i]=round(zitot)
  simulatedoutbreak9[,i]=round(zoutbreak)
}



#sigid10
for(i in 1:nsim){
  set.seed(i)
  yt=negbinNoise1(N=N,k=1,k2=1,alpha=5,beta=0,gama1=0.2,gama2=0.1,gama3=0.05,gama4=0.1,phi=1,shift=0,shift2=0)
  mu=exp(h1(N=N,k=1,k2=1,alpha=5,beta=0,gama1=0.2,gama2=0.1,gama3=0.05,gama4=0.1,shift=0,shift2=0))

  set.seed(i)
  out2=outbreak5(currentday=days5*52*years,weeklength=52*days5*years,wtime=(length(yt)-49*days5+1):length(yt),interval=0.25,yi=yt,k=1,k2=1,alpha=5,beta=0,gama1=0.2,
                 gama2=0.1,gama3=0.05,gama4=0.1,phi=1,shift=0,shift2=0,numoutbk=1,peakoutbk=w*days5,meanlog=0,sdlog=0.5)

  zoutbreak=out2$outbreak

  zt=yt#+out1
  zitot=yt + out2$outbreak #+out1

  #zitot[(bankhols5==1)]=0
  #zitot[(bankhols5==1)+1]=1.5*zitot[i+1]

  for(b in 1:length(zitot)){
    if(bankhols5[b]==1){
      zitot[b]=0
      zitot[b+1]=1.5*zitot[b+1]
    }
  }

  zeros=rep(0,2)
  weekend=seq(days5,days5*years*52,days5)
  #weekend=seq(0,days5*years*52-1,days5)
  for(s in 1:length(weekend)){
    zt=append(zt,zeros,after=2*(s-1)+weekend[s])
    zitot=append(zitot,zeros,after=2*(s-1)+weekend[s])
    zoutbreak=append(zoutbreak,zeros,after=2*(s-1)+weekend[s])
  }

  simulateddata10[,i]=round(zt)
  simulatedtotals10[,i]=round(zitot)
  simulatedoutbreak10[,i]=round(zoutbreak)
}


#sigid11
for(i in 1:nsim){
  set.seed(i)
  yt=negbinNoise1(N=N,k=2,k2=1,alpha=0.5,beta=0,gama1=0.4,gama2=0,gama3=0.05,gama4=0.15,phi=1,shift=0,shift2=0)
  #mu=exp(h1(N=N,k=2,k2=1,alpha=0.5,beta=0,gama1=0.4,gama2=0,gama3=0.05,gama4=0.15,shift=0,shift2=0))

  set.seed(i)
  out2=outbreak5(currentday=days5*52*years,weeklength=52*days5*years,wtime=(length(yt)-49*days5+1):length(yt),yi=yt,interval=0.25,k=2,k2=1,alpha=0.5,beta=0,gama1=0.4,
                 gama2=0,gama3=0.05,gama4=0.15,phi=1,shift=0,shift2=0,numoutbk=1,peakoutbk=w*days5,meanlog=0,sdlog=0.5)

  zoutbreak=out2$outbreak

  zt=yt#+out1
  zitot=yt + out2$outbreak #+out1

  #zitot[(bankhols5==1)]=0
  #zitot[(bankhols5==1)+1]=1.5*zitot[i+1]

  for(b in 1:length(zitot)){
    if(bankhols5[b]==1){
      zitot[b]=0
      zitot[b+1]=1.5*zitot[b+1]
    }
  }

  zeros=rep(0,2)
  weekend=seq(days5,days5*years*52,days5)
  #weekend=seq(0,days5*years*52-1,days5)
  for(s in 1:length(weekend)){
    zt=append(zt,zeros,after=2*(s-1)+weekend[s])
    zitot=append(zitot,zeros,after=2*(s-1)+weekend[s])
    zoutbreak=append(zoutbreak,zeros,after=2*(s-1)+weekend[s])
  }

  simulateddata11[,i]=round(zt)
  simulatedtotals11[,i]=round(zitot)
  simulatedoutbreak11[,i]=round(zoutbreak)
}


#sigid12
for(i in 1:nsim){
  set.seed(i)
  yt=negbinNoise1(N=N,k=1,k2=1,alpha=9,beta=0,gama1=0.5,gama2=0.2,gama3=0.2,gama4=0.5,phi=1,shift=0,shift2=0)/100
  #mu=exp(h1(N=N,k=1,k2=1,alpha=9,beta=0,gama1=0.5,gama2=0.2,gama3=0.2,gama4=0.5,shift=0,shift2=0))

  set.seed(i)
  out2=outbreak5(currentday=days5*52*years,weeklength=52*days5*years,wtime=(length(yt)-49*days5+1):length(yt),yi=yt,interval=0.25,k=1,k2=1,alpha=9,beta=0,gama1=0.5,
                 gama2=0.2,gama3=0.2,gama4=0.5,phi=1,shift=0,shift2=0,numoutbk=1,peakoutbk=w*days5,meanlog=0,sdlog=0.5)

  zoutbreak=out2$outbreak

  zt=yt#+out1
  zitot=yt + out2$outbreak #+out1

  #zitot[(bankhols5==1)]=0
  #zitot[(bankhols5==1)+1]=1.5*zitot[i+1]

  for(b in 1:length(zitot)){
    if(bankhols5[b]==1){
      zitot[b]=0
      zitot[b+1]=1.5*zitot[b+1]
    }
  }

  zeros=rep(0,2)
  weekend=seq(days5,days5*years*52,days5)
  #weekend=seq(0,days5*years*52-1,days5)
  for(s in 1:length(weekend)){
    zt=append(zt,zeros,after=2*(s-1)+weekend[s])
    zitot=append(zitot,zeros,after=2*(s-1)+weekend[s])
    zoutbreak=append(zoutbreak,zeros,after=2*(s-1)+weekend[s])
  }

  simulateddata12[,i]=round(zt)
  simulatedtotals12[,i]=round(zitot)
  simulatedoutbreak12[,i]=round(zoutbreak)
}

plot(1:length(simulatedtotals12[,1]),simulatedtotals12[,1],typ='l')

#plot(1:N,simulatedtotals12[,1],typ='l',xlim=c(2206,2548),col='green')
#lines(1:N,simulateddata12[,1],typ='l')


#=====================
# 7-day week syndromes
#=====================

years=7
days7=7
N=52*days7*years


#sigid1
print("Starting sigid1 simulations...")
for(i in 1:nsim){
  print(paste("sigid1: Simulation", i, "of", nsim))
  set.seed(i)
  yt=negbinNoise2(N=N,k=1,k2=2,alpha=6,beta=0,gama1=0.2,gama2=0.2,gama3=0.5,gama4=0.4,phi=2,shift=29)
  #mu=exp(h2(N=N,k=1,k2=2,alpha=6,beta=0,gama1=0.2,gama2=0.2,gama3=0.5,gama4=0.4,shift=29))

  set.seed(i)
  out2=outbreak7(currentday=N,weeklength=52*days7*years,wtime=(length(yt)-49*days7+1):length(yt),yi=yt,interval=0.25,k=1,k2=2,alpha=6,beta=0,gama1=0.2,gama2=0.2,
                 gama3=0.5,gama4=0.4,phi=2,shift=29,numoutbk=1,peakoutbk=w*days7,meanlog=0,sdlog=0.5)

  zoutbreak=out2$outbreak

  zt=yt#+out1
  zitot=yt + out2$outbreak #+out1

  for(b in 1:length(zitot)){
    if(bankhols7[b]==1){
      zitot[b]=2*zitot[b]
    }
  }

  simulateddata1[,i]=round(zt)
  simulatedtotals1[,i]=round(zitot)
  simulatedoutbreak1[,i]=round(zoutbreak)
}

#sigid2
for(i in 1:nsim){
  set.seed(i)
  yt=negbinNoise2(N=N,k=1,k2=2,alpha=0.5,beta=0,gama1=1.5,gama2=1.4,gama3=0.5,gama4=0.4,phi=1,shift=-167)
  #mu=exp(h2(N=N,k=1,k2=2,alpha=0.5,beta=0,gama1=1.5,gama2=1.4,gama3=0.5,gama4=0.4,shift=-167))

  set.seed(i)
  out2=outbreak7(currentday=N,weeklength=52*days7*years,wtime=(length(yt)-49*7+1):length(yt),yi=yt,interval=0.25,k=1,k2=2,alpha=0.5,beta=0,gama1=1.5,
                 gama2=1.4,gama3=0.5,gama4=0.4,phi=1,shift=-167,numoutbk=1,peakoutbk=w*days7,meanlog=0,sdlog=0.5)

  zoutbreak=out2$outbreak

  zt=yt#+out1
  zitot=yt + out2$outbreak #+out1

  for(b in 1:length(zitot)){
    if(bankhols7[b]==1){
      zitot[b]=2*zitot[b]
    }
  }

  simulateddata2[,i]=round(zt)
  simulatedtotals2[,i]=round(zitot)
  simulatedoutbreak2[,i]=round(zoutbreak)
}


#sigid3
for(i in 1:nsim){
  set.seed(i)
  yt=negbinNoise2(N=N,k=0,k2=2,alpha=5.5,beta=0,gama1=0,gama2=0,gama3=0.3,gama4=0.25,phi=1,shift=1)
  #mu=exp(h2(N=N,k=0,k2=2,alpha=5.5,beta=0,gama1=0,gama2=0,gama3=0.3,gama4=0.25,shift=1))

  set.seed(i)
  out2=outbreak7(currentday=N,weeklength=52*days7*12,wtime=(length(yt)-49*7+1):length(yt),yi=yt,interval=0.25,k=0,k2=2,alpha=5.5,beta=0,gama1=0,
                 gama2=0,gama3=0.3,gama4=0.25,phi=1,shift=1,numoutbk=1,peakoutbk=w*days7,meanlog=0,sdlog=0.5)

  zoutbreak=out2$outbreak

  zt=yt#+out1
  zitot=yt + out2$outbreak #+out1

  for(b in 1:length(zitot)){
    if(bankhols7[b]==1){
      zitot[b]=2*zitot[b]
    }
  }

  simulateddata3[,i]=round(zt)
  simulatedtotals3[,i]=round(zitot)
  simulatedoutbreak3[,i]=round(zoutbreak)
}


#sigid4
for(i in 1:nsim){
  set.seed(i)
  yt=negbinNoise2(N=N,k=0,k2=2,alpha=2,beta=0,gama1=0,gama2=0,gama3=0.3,gama4=0.25,phi=1,shift=1)
  #mu=exp(h2(N=N,k=0,k2=2,alpha=2,beta=0,gama1=0,gama2=0,gama3=0.3,gama4=0.25,shift=1))

  set.seed(i)
  out2=outbreak7(currentday=N,weeklength=52*days7*years,wtime=(length(yt)-49*days7+1):length(yt),yi=yt,interval=0.25,k=0,k2=2,alpha=2,beta=0,gama1=0,
                 gama2=0,gama3=0.3,gama4=0.25,phi=1,shift=1,numoutbk=1,peakoutbk=w*days7,meanlog=0,sdlog=0.5)

  zoutbreak=out2$outbreak

  zt=yt#+out1
  zitot=yt + out2$outbreak #+out1

  for(b in 1:length(zitot)){
    if(bankhols7[b]==1){
      zitot[b]=2*zitot[b]
    }
  }

  simulateddata4[,i]=round(zt)
  simulatedtotals4[,i]=round(zitot)
  simulatedoutbreak4[,i]=round(zoutbreak)
}


#sigid13
for(i in 1:nsim){
  set.seed(i)
  yt=negbinNoise2(N=N,k=1,k2=2,alpha=2,beta=0.0005,gama1=0.8,gama2=0.8,gama3=0.8,gama4=0.4,phi=4,shift=57)
  #mu=exp(h2(N=N,k=1,k2=2,alpha=2,beta=0,gama1=0.8,gama2=0.8,gama3=0.8,gama4=0.4,shift=57))

  set.seed(i)
  out2=outbreak7(currentday=N,weeklength=52*days7*years,wtime=(length(yt)-49*days7+1):length(yt),yi=yt,interval=0.25,k=1,k2=2,alpha=6,beta=0,gama1=0.2,gama2=0.2,
                 gama3=0.5,gama4=0.4,phi=2,shift=29,numoutbk=1,peakoutbk=w*days7,meanlog=0,sdlog=0.5)

  zoutbreak=out2$outbreak

  zt=yt#+out1
  zitot=yt + out2$outbreak #+out1

  for(b in 1:length(zitot)){
    if(bankhols7[b]==1){
      zitot[b]=2*zitot[b]
    }
  }

  simulateddata13[,i]=round(zt)
  simulatedtotals13[,i]=round(zitot)
  simulatedoutbreak13[,i]=round(zoutbreak)
}


#sigid14
for(i in 1:nsim){
  set.seed(i)
  #yt=0.1*(negbinNoise2(N=N,k=4,k2=1,alpha=1.5,beta=0,gama1=0.1,gama2=0.1,gama3=1.8,gama4=0.1,phi=1,shift=-85)+2)

  yt=1*(negbinNoise2(N=N,k=4,k2=1,alpha=0.05,beta=0,gama1=0.01,gama2=0.01,gama3=1.8,gama4=0.1,phi=1,shift=-85)+0)

  set.seed(i)
  out2=outbreak7(currentday=N,weeklength=52*days7*years,wtime=(length(yt)-49*days7+1):length(yt),yi=yt,interval=0.25,k=1,k2=2,alpha=2,beta=0,gama1=0.8,
                 gama2=0.8,gama3=0.8,gama4=0.4,phi=4,shift=57,numoutbk=1,peakoutbk=w*days7,meanlog=0,sdlog=0.5)

  zoutbreak=out2$outbreak

  zt=yt#+out1
  zitot=yt + out2$outbreak #+out1

  for(b in 1:length(zitot)){
    if(bankhols7[b]==1){
      zitot[b]=2*zitot[b]
    }
  }

  simulateddata14[,i]=round(zt)
  simulatedtotals14[,i]=round(zitot)
  simulatedoutbreak14[,i]=round(zoutbreak)
}


#plot(1:N,yt,typ='l')
#plot(1:(52*years*7),zt,typ='l',xlim=c(1,364))
#plot(1:(52*years*7),zitot,typ='l',xlim=c(1,364))
#plot(1:(52*years*7),zitot,typ='l',xlim=c(365,728))
#plot(1:(52*years*7),zitot,typ='l',xlim=c(729,1092))
#plot(1:(52*years*7),zitot,typ='l',xlim=c(1093,1456))
#plot(1:(52*years*7),zitot,typ='l',xlim=c(1457,1820))
#plot(1:(52*years*7),zitot,typ='l',xlim=c(1821,2184))
#plot(1:(52*years*7),zitot,typ='l',xlim=c(2185,2548))
#lines(1:(52*years*7),zoutbreak,col='green')

plot(1:(52*years*7),simulatedtotals5[,1],typ='l',xlim=c(1,7*364))
lines(1:(52*years*7),simulatedzseasoutbreak5[,1],col='green')
lines(1:(52*years*7),simulatedoutbreak5[,1],col='red')

plot(1:(52*years*7),zitot,typ='l',xlim=c(1,7*364))
lines(1:(52*years*7),yt,col='blue')
lines(1:(52*years*7),zseasoutbreak,col='green')
lines(1:(52*years*7),zoutbreak,col='red')

plot(1:(52*years*7),simulatedzseasoutbreak5[,1],col='green',typ='l')


#sigid15
print("Starting sigid15 simulations...")
for(i in 1:nsim){
  print(paste("sigid15: Simulation", i, "of", nsim))
  set.seed(i)
  yt=negbinNoise2(N=N,k=1,k2=2,alpha=3,beta=0,gama1=0.8,gama2=0.6,gama3=0.8,gama4=0.4,phi=4,shift=29)
  #mu=exp(h2(N=N,k=1,k2=2,alpha=3,beta=0,gama1=0.8,gama2=0.6,gama3=0.8,gama4=0.4,shift=29))

  out1=rep(0,N)
  for(j in 1:years){
    set.seed(j+years*i)
    out=outbreak5(currentday=days7*52*years,weeklength=52*days7*years,wtime=((210+(j-1)*days7*52):(230+(j-1)*days7*52)),yi=yt,interval=0.02,k=1,k2=1,alpha=1,beta=0,gama1=0.1,
                  gama2=2,gama3=0.05,gama4=0.05,phi=1,shift=-50,shift2=-50,numoutbk=1,peakoutbk=3*days7*150,meanlog=0,sdlog=0.5)
    out1=out1+out$outbreak
  }

  set.seed(i)
  out2=outbreak7(currentday=N,weeklength=52*days7*years,wtime=(length(yt)-49*days7+1):length(yt),yi=yt,interval=0.25,k=1,k2=2,alpha=3,beta=0,gama1=0.8,
                 gama2=0.6,gama3=0.8,gama4=0.4,phi=4,shift=29,numoutbk=1,peakoutbk=w*days7,meanlog=0,sdlog=0.5)

  zoutbreak=out2$outbreak
  zseasoutbreak=out2$outbreak+out1

  zt=yt +out1
  zitot=yt + out2$outbreak +out1

  for(b in 1:length(zitot)){
    if(bankhols7[b]==1){
      zitot[b]=2*zitot[b]
    }
  }

  simulateddata15[,i]=round(zt)
  simulatedtotals15[,i]=round(zitot)
  simulatedoutbreak15[,i]=round(zoutbreak)
  simulatedzseasoutbreak15[,i]=round(zseasoutbreak)
}


plot(1:(52*years*7),zitot,typ='l',xlim=c(1,7*364))
lines(1:(52*years*7),yt,col='blue')
lines(1:(52*years*7),zseasoutbreak,col='green')
lines(1:(52*years*7),zoutbreak,col='red')

plot(1:(52*years*7),simulatedtotals15[,1],typ='l',xlim=c(1,7*364))
lines(1:(52*years*7),simulatedzseasoutbreak15[,1],col='green')
lines(1:(52*years*7),simulatedoutbreak15[,1],col='red')


plot(1:(52*years*7),simulatedtotals15[,1],typ='l',xlim=c(1,7*364))
lines(1:(52*years*7),simulatedzseasoutbreak15[,1],col='green')
lines(1:(52*years*7),simulatedoutbreak15[,1],col='red')


#sigid16
for(i in 1:nsim){
  set.seed(i)
  yt=negbinNoise2(N=N,k=0,k2=2,alpha=6,beta=0,gama1=0,gama2=0,gama3=0.8,gama4=0.4,phi=4,shift=1)
  #mu=exp(h2(N=N,k=0,k2=2,alpha=6,beta=0,gama1=0,gama2=0,gama3=0.8,gama4=0.4,shift=1))

  set.seed(i)
  out2=outbreak7(currentday=N,weeklength=52*7*12,wtime=(length(yt)-49*days7+1):length(yt),yi=yt,interval=0.25,k=0,k2=2,alpha=6,beta=0,gama1=0,
                 gama2=0,gama3=0.8,gama4=0.4,phi=4,shift=1,numoutbk=1,peakoutbk=w*days7,meanlog=0,sdlog=0.5)

  zoutbreak=out2$outbreak

  zt=yt#+out1
  zitot=yt + out2$outbreak #+out1

  for(b in 1:length(zitot)){
    if(bankhols7[b]==1){
      zitot[b]=2*zitot[b]
    }
  }

  simulateddata16[,i]=round(zt)
  simulatedtotals16[,i]=round(zitot)
  simulatedoutbreak16[,i]=round(zoutbreak)
}



###############
#PLOTS###########
################
years=7
i=2
par(mfrow=(c(4,4)))
plot(1:(52*years*7),simulateddata1[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal1')
plot(1:(52*years*7),simulateddata2[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal2')
plot(1:(52*years*7),simulateddata3[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal3')
plot(1:(52*years*7),simulateddata4[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal4')
plot(1:(52*years*7),simulateddata5[,i]-simulatedzseasoutbreak5[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal5')
plot(1:(52*years*7),simulateddata6[,i]-simulatedzseasoutbreak6[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal6')
plot(1:(52*years*7),simulateddata7[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal7')
plot(1:(52*years*7),simulateddata8[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal8')
plot(1:(52*years*7),simulateddata9[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal9')
plot(1:(52*years*7),simulateddata10[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal10')
plot(1:(52*years*7),simulateddata11[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal11')
plot(1:(52*years*7),simulateddata12[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal12')
plot(1:(52*years*7),simulateddata13[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal13')
plot(1:(52*years*7),simulateddata14[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal14')
plot(1:(52*years*7),simulateddata15[,i]-simulatedzseasoutbreak15[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal15')
plot(1:(52*years*7),simulateddata16[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal16')


par(mfrow=(c(2,1)))
plot(1:(52*years*7),simulateddata3[,1],typ='l',xlim=c(1,3*7),xlab='day',ylab='count',main='signal3',xaxt="n")
axis(1,at=1:21,1:21,las=2)
plot(1:(52*years*7),simulateddata7[,1],typ='l',xlim=c(1,3*7),xlab='day',ylab='count',main='signal7',xaxt="n")
axis(1,at=1:21,1:21,las=2)

i=98
par(mfrow=(c(4,4)))
plot(1:(52*years*7),simulateddata1[,i]+simulatedoutbreak1[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal1',col='blue')
lines(1:(52*years*7),simulateddata1[,i],col='black')
plot(1:(52*years*7),simulateddata2[,i]+simulatedoutbreak2[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal2',col='blue')
lines(1:(52*years*7),simulateddata2[,i],col='black')
plot(1:(52*years*7),simulateddata3[,i]+simulatedoutbreak3[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal3',col='blue')
lines(1:(52*years*7),simulateddata3[,i],col='black')
plot(1:(52*years*7),simulateddata4[,i]+simulatedoutbreak4[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal4',col='blue')
lines(1:(52*years*7),simulateddata4[,i],col='black')
plot(1:(52*years*7),simulateddata5[,i]-simulatedzseasoutbreak5[,i]+simulatedoutbreak5[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal5',col='blue')
lines(1:(52*years*7),simulateddata5[,i]-simulatedzseasoutbreak5[,i],col='black')
plot(1:(52*years*7),simulateddata6[,i]-simulatedzseasoutbreak6[,i]+simulatedoutbreak6[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal6',col='blue')
lines(1:(52*years*7),simulateddata6[,i]-simulatedzseasoutbreak6[,i],col='black')
plot(1:(52*years*7),simulateddata7[,i]+simulatedoutbreak7[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal7',col='blue')
lines(1:(52*years*7),simulateddata7[,i],col='black')
plot(1:(52*years*7),simulateddata8[,i]+simulatedoutbreak8[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal8',col='blue')
lines(1:(52*years*7),simulateddata8[,i],col='black')
plot(1:(52*years*7),simulateddata9[,i]+simulatedoutbreak9[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal9',col='blue')
lines(1:(52*years*7),simulateddata9[,i],col='black')
plot(1:(52*years*7),simulateddata10[,i]+simulatedoutbreak10[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal10',col='blue')
lines(1:(52*years*7),simulateddata10[,i],col='black')
plot(1:(52*years*7),simulateddata11[,i]+simulatedoutbreak11[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal11',col='blue')
lines(1:(52*years*7),simulateddata11[,i],col='black')
plot(1:(52*years*7),simulateddata12[,i]+simulatedoutbreak12[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal12',col='blue')
lines(1:(52*years*7),simulateddata12[,i],col='black')
plot(1:(52*years*7),simulateddata13[,i]+simulatedoutbreak13[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal13',col='blue')
lines(1:(52*years*7),simulateddata13[,i],col='black')
plot(1:(52*years*7),simulateddata14[,i]+simulatedoutbreak14[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal14',col='blue')
lines(1:(52*years*7),simulateddata14[,i],col='black')
plot(1:(52*years*7),simulateddata15[,i]-simulatedzseasoutbreak15[,i]+simulatedoutbreak15[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal15',col='blue')
lines(1:(52*years*7),simulateddata15[,i]-simulatedzseasoutbreak15[,i],col='black')
plot(1:(52*years*7),simulateddata16[,i]+simulatedoutbreak16[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal16',col='blue')
lines(1:(52*years*7),simulateddata16[,i],col='black')



par(mfrow=(c(3,1)))
plot(1:(52*years*7),simulateddata5[,1],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal5',col='blue')
lines(1:(52*years*7),simulateddata5[,1]-simulatedzseasoutbreak5[,1],col='black')

plot(1:(52*years*7),simulateddata6[,1],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal6',col='blue')
lines(1:(52*years*7),simulateddata6[,1]-simulatedzseasoutbreak6[,1],col='black')

plot(1:(52*years*7),simulateddata15[,1],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal15',col='blue')
lines(1:(52*years*7),simulateddata15[,1]-simulatedzseasoutbreak15[,1],col='black')


par(mfrow=(c(3,1)))
plot(1:(52*years*7),simulateddata5[,1],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal5',col='blue')
lines(1:(52*years*7),simulateddata5[,1]-simulatedzseasoutbreak5[,1],col='black')

plot(1:(52*years*7),simulateddata6[,1],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal6',col='blue')
lines(1:(52*years*7),simulateddata6[,1]-simulatedzseasoutbreak6[,1],col='black')

plot(1:(52*years*7),simulateddata15[,1],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal15',col='blue')
lines(1:(52*years*7),simulateddata15[,1]-simulatedzseasoutbreak15[,1],col='black')






plot(1:(52*years*7),zitot,typ='l',xlim=c(1,7*364))
lines(1:(52*years*7),yt,col='blue')
lines(1:(52*years*7),zseasoutbreak,col='green')
lines(1:(52*years*7),zoutbreak,col='red')

plot(1:(52*years*7),simulatedtotals15[,1],typ='l',xlim=c(1,7*364))
lines(1:(52*years*7),simulatedzseasoutbreak15[,1],col='green')
lines(1:(52*years*7),simulatedoutbreak15[,1],col='red')


plot(1:(52*years*7),simulatedtotals15[,1],typ='l',xlim=c(1,7*364))
lines(1:(52*years*7),simulatedzseasoutbreak15[,1],col='green')
lines(1:(52*years*7),simulatedoutbreak15[,1],col='red')


years=7
i=2
par(mfrow=(c(4,4)))
plot(1:(52*years*7),simulatedtotals1[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal1')
plot(1:(52*years*7),simulatedtotals2[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal2')
plot(1:(52*years*7),simulatedtotals3[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal3')
plot(1:(52*years*7),simulatedtotals4[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal4')
plot(1:(52*years*7),simulatedtotals5[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal5')
plot(1:(52*years*7),simulatedtotals6[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal6')
plot(1:(52*years*7),simulatedtotals7[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal7')
plot(1:(52*years*7),simulatedtotals8[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal8')
plot(1:(52*years*7),simulatedtotals9[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal9')
plot(1:(52*years*7),simulatedtotals10[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal10')
plot(1:(52*years*7),simulatedtotals11[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal11')
plot(1:(52*years*7),simulatedtotals12[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal12')
plot(1:(52*years*7),simulatedtotals13[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal13')
plot(1:(52*years*7),simulatedtotals14[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal14')
plot(1:(52*years*7),simulatedtotals15[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal15')
plot(1:(52*years*7),simulatedtotals16[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal16')

# Create a directory to store the datasets
dir.create("signal_datasets", showWarnings = FALSE)

# Loop through all 16 signals
for (sig_num in 1:16) {
  # Define file names for the current signal
  totals_filename <- paste0("signal_datasets/simulated_totals_sig", sig_num, ".csv")
  outbreaks_filename <- paste0("signal_datasets/simulated_outbreaks_sig", sig_num, ".csv")

  # Construct the variable names dynamically
  simulatedtotals_var <- paste0("simulatedtotals", sig_num)
  simulatedoutbreak_var <- paste0("simulatedoutbreak", sig_num)

  # Write the totals and test-set outbreaks to CSV
  write.csv(get(simulatedtotals_var), totals_filename, row.names = FALSE)
  write.csv(get(simulatedoutbreak_var), outbreaks_filename, row.names = FALSE)

  cat("Saved:", totals_filename, "and", outbreaks_filename, "\n")
}

# Handle the special seasonal outbreak files
write.csv(simulatedzseasoutbreak5, "signal_datasets/simulated_seasonal_outbreaks_sig5.csv", row.names = FALSE)
write.csv(simulatedzseasoutbreak6, "signal_datasets/simulated_seasonal_outbreaks_sig6.csv", row.names = FALSE)
write.csv(simulatedzseasoutbreak15, "signal_datasets/simulated_seasonal_outbreaks_sig15.csv", row.names = FALSE)

cat("Saved: Seasonal outbreak files for signals 5, 6, and 15.\n")
print("All datasets have been saved to the 'signal_datasets' folder.")
