#THINKS thermokinetic software, version #01.07.20 Website: thinks.chemphys.ru. To cite: https://doi.org/10.3390/molecules24122298
#This routine performs the predition of thermal behavior using the model-fitting or the isoconversional kinetic data.
#
#        Δ  ****  __|__  ****  Δ
#       / \      /.....\      / \
#       |Π|_____|_______|_____|Π|
#      /   \  ■ | T T T | ■  /   \
#     | ICP | ■ | | | | | ■ | RAS |
#     |_____|___|_|Π|Π|_|___|_____|
#
#Created by Dr. N.V. Muravyev (2016-2020).
#--------------------------------------------------------------------------------------------------------------------------

library(shiny)
library(deSolve)

options(shiny.maxRequestSize=2*1024^2)  #change the file-size limit here----in the server there is 10 MB limit

shinyServer(function(input, output,session) {
  
  session$onSessionEnded(function() {
    observe({ 
      stopApp() })
  })
  
    observeEvent(input$go, {
      updateTabsetPanel(session, "MainTabset",
                        selected = "panel1")
    }, once=TRUE)
    
    observe( {
      if ((input$KinT==1)&&(input$KinMod<3)) hideTab(inputId = "tabs", target = "Second stage")
      if ((input$KinT==1)&&(input$KinMod<9)) hideTab(inputId = "tabs", target = "Third stage")
            
      if ((input$KinT==1)&&(input$KinMod>=3)) showTab(inputId = "tabs", target = "Second stage")
      if ((input$KinT==1)&&(input$KinMod>=9)) showTab(inputId = "tabs", target = "Third stage")
      
      if ((input$KinT==1)&&(input$KinMod>=4)) {
        if ((input$KinMod==6)||(input$KinMod==9)) {
          showConc=c("A"=0,"B"=1,"C"=2,"D"=3)
          selConc=c(0,1,2,3) }
        else {
          showConc=c("A"=0,"B"=1,"C"=2)
          selConc=c(0,1,2) }
        updateSelectInput(session, "ConcLeft", choices = showConc, selected = selConc)
        updateSelectInput(session, "ConcRight", choices = showConc) }
      
    })
    
  Isoconv <- reactive({ #----------------------------------reads file with the isoconversional results
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    data=read.table(inFile$datapath,header = TRUE, check.names = FALSE)
    data=data[data[[1]]>input$aMin & data[[1]]<input$aMax,]
    al=data[[1]]
    Ea=data[[2]]
    dEa=data[[3]]
    lnf=data[[4]]
    dlnf=data[[5]]
    R2=data[[6]]
    list(al,Ea,dEa,lnf,dlnf,R2)
  })

  Prediction=reactive({ #Thermal behavior prediction
    Npts=input$Npts
    Npts=ifelse(((Npts>=200)&&(Npts<=2000)),Npts,300) #----------Restriction only for online!
    #------------Building a sequence of (Temperature, Time) data
    if (input$Tprg==1) {
      if (input$TimeSc==1) Temp=seq(input$THRini,input$THRfin, length.out=Npts) # in °C #linear sequence
      if (input$TimeSc==2) Temp=10^seq(log10(input$THRini),log10(input$THRfin), length.out=Npts) # in °C #logarithmic sequence
      ti=(Temp-min(Temp))/input$THR  }
    if (input$Tprg==2) {
      log10timeMin=-12
      if (input$TimeSc==2) ti=10^seq(log10timeMin,log10(input$Timefin), length.out = Npts)   #logarithmic sequence
      if (input$TimeSc==1) ti=seq(0,input$Timefin, length.out = Npts)   #linear sequence
      Temp=ti*0+input$Tiso }
    if (input$Tprg==3) {
      inFile2 <- input$file2
      if (is.null(inFile2)) return(NULL)
      TprgData=read.table(inFile2$datapath, skip=1) # 2nd row - t[s], 1st - T[°C]
      ti=TprgData[,2]/60
      Temp=TprgData[,1]
      userTt=suppressWarnings(splinefun(ti,Temp,method = "fmm"))      #-----------check using of other functions?
      lgTmin=ifelse(min(ti)>0,log10(min(ti)),log10timeMin)
      if (input$TimeSc==2) ti=10^seq(lgTmin,log10(max(ti)), length.out = Npts)   #logarithmic sequence
      if (input$TimeSc==1) ti=seq(min(ti),max(ti), length.out = Npts)    #linear sequence
      Temp=userTt(ti)  }
    Temp=Temp+273.15 # Now temperature in K
    ti=ti*60 # Now time in s
    al=numeric(Npts)
    dadt=numeric(Npts)
    al[1]=input$a0
    dadt[1]=0
    if (input$KinT==2) { #---------------------------------------solver for isoconversional kinetic data
      inFile <- input$file1
      if (is.null(inFile)) return(NULL)
      res=Isoconv()
      if(input$KinFT==1) Ea=splinefun(res[[1]],res[[2]],method = "fmm")
      if(input$KinFT==2) Ea=splinefun(res[[1]],res[[2]],method = "natural")
      if(input$KinFT==3) Ea=approxfun(res[[1]],res[[2]],method = "linear",rule=2)
      if(input$KinFT==1) lnAf=splinefun(res[[1]],res[[4]],method = "fmm")
      if(input$KinFT==2) lnAf=splinefun(res[[1]],res[[4]],method = "natural")
      if(input$KinFT==3) lnAf=approxfun(res[[1]],res[[4]],method = "linear",rule=2) 
      for (i in 2:Npts) {
        al[i]=al[i-1]+dadt[i-1]*(ti[i]-ti[i-1])
        dadt[i]=exp(lnAf(al[i])-Ea(al[i])*1000/8.314/Temp[i]) 
        if(is.na(al[i])) {
          al[i]=1
          dadt[i]=0 }
        if(al[i]>1) {
          al[i]=1
          dadt[i]=0 }
        if(al[i]<0) {
          al[i]=0 
          dadt[i]=0 }
        }
      return(list(ti,Temp,al,dadt))
      }
    else {  #---------------------------------------------------------------External solver for model-fitting data
      TempF=suppressWarnings(splinefun(ti,Temp,method = "fmm")) 
      
#      if(input$KinMod==1)  #---------------------------------------------------(V1) Flexible single-step (ePT)
#        resT=funK1(input$lA1,input$Ea,input$n1,input$m1,
#                   Temp,ti,0,0,0,0,input$methodDiff,input$q,input$a0,1,input$suppression)
#      else 
        resT=NLR_Prediction(Npts,ti,TempF,input$methodDiff,input$suppression)
      return(resT)
      }
  })
  
  NLR_Prediction=function(Npts,ti,TempF,mD,supY) {
    MyNull=1E-11
    MyUnity=1-MyNull
    MyNull_for_Log=3E-18
    a0=input$a0
    KinModel=input$KinMod
    if(KinModel==1) {  #---------------------------------------------------(V1) Flexible single-step (ePT/KJMAE)
      if (input$m1!=10) {
        v1Func <- function(t, x, parms) {
          with(as.list(parms), {
            dx <- exp(lA1-Ea*1000/8.31446/TempF(t))*((1-x)^n1)*(1-q*(1-x))^m1
            list(dx)
          }) }
        parms=c(lA1=input$lA1, Ea=input$Ea, n1=input$n1, m1=input$m1, q=input$q) }
      else {
        v1Func <- function(t, x, parms) { #----reaction model in extended KJMAE form
          with(as.list(parms), {
            x=ifelse(x<MyUnity,x,MyUnity)
            dx <- exp(lA1-Ea*1000/8.31446/TempF(t))*n1*(1-x)*(-log(1-x))^((n1-1)/n1)
            list(dx)
          }) } 
        parms=c(lA1=input$lA1, Ea=input$Ea, n1=input$n1)
        }
      x=c(a = a0)
      vFunc=v1Func
      }
    if(KinModel==2) {  #--------------------------------------------------(V2) Simple reaction types
      rMod=input$reMod
      v2Func <- function(t, x, parms) {
        with(as.list(parms), {
          dx <- exp(lA1-Ea*1000/8.31446/TempF(t))*reactionModel(rMod,x)
          list(dx)
        }) }
      x=c(a = a0)
      parms=c(lA1=input$lA1, Ea=input$Ea) 
      vFunc=v2Func }
    if(KinModel==3) {  #---------------------------------------------------(V3) Autocatalytic (DMM) model
      v3Func <- function(t, x, parms) {
        with(as.list(parms), {
          dx = exp(lA1-Ea*1000/8.31446/TempF(t))*(1-x)+exp(lA2-Ea2*1000/8.31446/TempF(t))*(1-mu)*x*(1-x)/(1-mu*x)
          list(dx)
        }) }
      x=c(a = a0)
      parms=c(lA1=input$lA1, Ea=input$Ea, mu=input$mu, lA2=input$lA2, Ea2=input$Ea2)
      vFunc=v3Func }
    if(KinModel==4) {  #-------------------------------------------------(v4) Parallel Flexible steps (ePT/KJMAE || ePT) 
      q2=input$q2 #ifelse(input$q2f==T,input$q2,input$q)
#     v4Func <- function(t, x, parms) {
#        with(as.list(parms), {
#          dx <- exp(lA1-Ea*1000/8.31446/TempF(t))*((1-x)^n1)*(1-q*(1-x))^m1 +
#                  exp(lA2-Ea2*1000/8.31446/TempF(t))*((1-x)^n2)*(1-q2*(1-x))^m2
#          list(dx)
#        }) }
      if (input$m1!=10) {
          v4Func <- function(t, x, parms) {
            with(as.list(c(parms, x)),{
              A=ifelse(A>MyNull,A,MyNull)
              B=ifelse(B>MyNull,B,MyNull)
              dA = -exp(lA1-Ea*1000/8.31446/TempF(t))*(A^n1)*(1-q*A)^m1 -exp(lA2-Ea2*1000/8.31446/TempF(t))*(A^n2)*(1-q2*A)^m2
              dB=exp(lA1-Ea*1000/8.31446/TempF(t))*(A^n1)*(1-q*A)^m1
              return(list(c(dA,dB)))
            }) }
          parms=c(lA1=input$lA1, Ea=input$Ea, n1=input$n1, m1=input$m1, q=input$q,
                  lA2=input$lA2, Ea2=input$Ea2, n2=input$n2, m2=input$m2, q2=q2) }
      else {
        v4Func <- function(t, x, parms) { #----reaction model in extended KJMAE form
          with(as.list(c(parms, x)),{
            A=ifelse(A>MyNull,A,MyNull)
            B=ifelse(B>MyNull,B,MyNull)
            dA = -exp(lA1-Ea*1000/8.31446/TempF(t))*n1*A*(-log(A))^((n1-1)/n1) -exp(lA2-Ea2*1000/8.31446/TempF(t))*(A^n2)*(1-q2*A)^m2
            dB=exp(lA1-Ea*1000/8.31446/TempF(t))*n1*A*(-log(A))^((n1-1)/n1)
            return(list(c(dA,dB)))
          }) } 
        parms=c(lA1=input$lA1, Ea=input$Ea, n1=input$n1, 
                lA2=input$lA2, Ea2=input$Ea2, n2=input$n2, m2=input$m2, q2=q2) }
      x=c(A=1-a0, B=a0)
      vFunc=v4Func }
    if(KinModel==5) {  #--------------------------------------------------(v5) Consecutive Flexible steps (ePT/KJMAE->ePT)
      q2=input$q2 #ifelse(input$q2f==T,input$q2,input$q)
      if (input$m1!=10) {
        v5Func <- function(t, x, parms) {
          with(as.list(c(parms, x)),{
            A=ifelse(A>MyNull,A,MyNull)
            B=ifelse(B>MyNull,B,MyNull)
            C=ifelse(C>MyNull,C,MyNull)
            dA = -exp(lA1-Ea*1000/8.31446/TempF(t))*(A^n1)*(1-q*(1-B))^m1 #---note that here A = 1 - a !
            dB = exp(lA1-Ea*1000/8.31446/TempF(t))*A^n1*(1-q*(1-B))^m1 - exp(lA2-Ea2*1000/8.31446/TempF(t))*B^n2*(1-q2*(1-C))^m2
            dC = exp(lA2-Ea2*1000/8.31446/TempF(t))*B^n2*(1-q2*(1-C))^m2 #1-A-B
            return(list(c(dA,dB,dC)))
          }) }
        parms=c(lA1=input$lA1, Ea=input$Ea, n1=input$n1, m1=input$m1, q=input$q,
                lA2=input$lA2, Ea2=input$Ea2, n2=input$n2, m2=input$m2, q2=q2) }
      else {
        v5Func <- function(t, x, parms) {
          with(as.list(c(parms, x)),{
            A=ifelse(A>MyNull,A,MyNull)
            B=ifelse(B>MyNull,B,MyNull)
            C=ifelse(C>MyNull,C,MyNull)
            dA = -exp(lA1-Ea*1000/8.31446/TempF(t))*n1*A*(-log(A))^((n1-1)/n1) #---note that here A = 1 - a !
            dB = exp(lA1-Ea*1000/8.31446/TempF(t))*n1*A*(-log(A))^((n1-1)/n1) - exp(lA2-Ea2*1000/8.31446/TempF(t))*B^n2*(1-q2*(1-C))^m2
            dC = exp(lA2-Ea2*1000/8.31446/TempF(t))*B^n2*(1-q2*(1-C))^m2 #1-A-B
            return(list(c(dA,dB,dC)))
          }) }
        parms=c(lA1=input$lA1, Ea=input$Ea, n1=input$n1, 
                lA2=input$lA2, Ea2=input$Ea2, n2=input$n2, m2=input$m2, q2=q2) }
      x=c(A=1-a0, B=a0, C=a0)
      vFunc=v5Func }
    if(KinModel==6) {  #--------------------------------------------------(v6) Independent Flexible steps (ePT/KJMAE+ePT)
      q2=input$q2 #ifelse(input$q2f==T,input$q2,input$q)
      if (input$m1!=10) {
        v6Func <- function(t, x, parms) {
          with(as.list(c(parms, x)),{
            A=ifelse(A>MyNull,A,MyNull)
            C=ifelse(C>MyNull,C,MyNull)
            dA = -exp(lA1-Ea*1000/8.31446/TempF(t))*(A^n1)*(1-q*A)^m1 #---note that here A = 1 - a !
            dC = -exp(lA2-Ea2*1000/8.31446/TempF(t))*C^n2*(1-q2*C)^m2
            return(list(c(dA,dC)))
          }) }
        parms=c(lA1=input$lA1, Ea=input$Ea, n1=input$n1, m1=input$m1, q=input$q,
                lA2=input$lA2, Ea2=input$Ea2, n2=input$n2, m2=input$m2, q2=q2) }
      else {
        v6Func <- function(t, x, parms) {
          with(as.list(c(parms, x)),{
            A=ifelse(A>MyNull,A,MyNull)
            C=ifelse(C>MyNull,C,MyNull)
            dA = -exp(lA1-Ea*1000/8.31446/TempF(t))*n1*A*(-log(A))^((n1-1)/n1) #---note that here A = 1 - a !
            dC = -exp(lA2-Ea2*1000/8.31446/TempF(t))*C^n2*(1-q2*C)^m2
            return(list(c(dA,dC)))
          }) }
        parms=c(lA1=input$lA1, Ea=input$Ea, n1=input$n1, 
                lA2=input$lA2, Ea2=input$Ea2, n2=input$n2, m2=input$m2, q2=q2) }
      x=c(A=1-a0, C=1-a0)
      vFunc=v6Func }
    #7
    
    #8
    
    if(KinModel==9) {  #--------------------------------------------------(v9) Flexible steps (ePT/KJMAE->ePT||ePT)  
      if (input$m1!=10) {
        v9Func <- function(t, x, parms) {
          with(as.list(c(parms, x)),{
            A=ifelse(A>MyNull,A,MyNull)
            B=ifelse(B>MyNull,B,MyNull)
            C=ifelse(C>MyNull,C,MyNull)
            D=ifelse(D>MyNull,D,MyNull)
            dA = -exp(lnA1_ode-Ea1_ode*1000/8.31446/TempF(t))*(A^n1_ode)*(1-q1_ode*(1-B))^m1_ode #---note that here A = 1 - a !
            dB = exp(lnA1_ode-Ea1_ode*1000/8.31446/TempF(t))*A^n1_ode*(1-q1_ode*(1-B))^m1_ode - 
              exp(lnA2_ode-Ea2_ode*1000/8.31446/TempF(t))*B^n2_ode*(1-q2_ode*(1-C))^m2_ode -
              exp(lnA3_ode-Ea3_ode*1000/8.31446/TempF(t))*B^n3_ode*(1-q3_ode*(1-D))^m3_ode 
            dC = exp(lnA2_ode-Ea2_ode*1000/8.31446/TempF(t))*B^n2_ode*(1-q2_ode*(1-C))^m2_ode 
            dD = exp(lnA3_ode-Ea3_ode*1000/8.31446/TempF(t))*B^n3_ode*(1-q3_ode*(1-D))^m3_ode 
            #D=1-A-B-C
            return(list(c(dA,dB,dC,dD)))
          }) }
        parms=c(lnA1_ode=input$lA1, Ea1_ode=input$Ea, n1_ode=input$n1, m1_ode=input$m1, q1_ode=input$q,
                lnA2_ode=input$lA2, Ea2_ode=input$Ea2, n2_ode=input$n2, m2_ode=input$m2, q2_ode=input$q2,
                lnA3_ode=input$lA3, Ea3_ode=input$Ea3, n3_ode=input$n3, m3_ode=input$m3, q3_ode=input$q3)
        vFunc=v9Func
      }
      else {
        v9Func <- function(t, x, parms) {
          with(as.list(c(parms, x)),{
            A_corr=ifelse(A>=MyNull_for_Log,A,MyNull_for_Log)
            B_corr=ifelse(B>=MyNull,B,MyNull)
            #C_corr=C#ifelse(C>=a0,ifelse(C<1,C,1),a0)
            dA = -exp(lnA1_ode-Ea1_ode*1000/8.31446/TempF(t))*n1_ode*A_corr*(-log(A_corr))^((n1_ode-1)/n1_ode) #---note that here A = 1 - a !
            dB = exp(lnA1_ode-Ea1_ode*1000/8.31446/TempF(t))*n1_ode*A_corr*(-log(A_corr))^((n1_ode-1)/n1_ode) - 
              exp(lnA2_ode-Ea2_ode*1000/8.31446/TempF(t))*B_corr^n2_ode*(1-q2_ode*(1-C))^m2_ode -
              exp(lnA3_ode-Ea3_ode*1000/8.31446/TempF(t))*B_corr^n3_ode*(1-q3_ode*(1-D))^m3_ode 
            dC = exp(lnA2_ode-Ea2_ode*1000/8.31446/TempF(t))*B_corr^n2_ode*(1-q2_ode*(1-C))^m2_ode 
            dD = exp(lnA3_ode-Ea3_ode*1000/8.31446/TempF(t))*B_corr^n3_ode*(1-q3_ode*(1-D))^m3_ode 
            #D=1-A_corr-B_corr-C
            return(list(c(dA,dB,dC,dD),A_corr=A_corr,B_corr=B_corr))
          }) }
        parms=c(lnA1_ode=input$lA1, Ea1_ode=input$Ea, n1_ode=input$n1, 
                lnA2_ode=input$lA2, Ea2_ode=input$Ea2, n2_ode=input$n2, m2_ode=input$m2, q2_ode=input$q2,
                lnA3_ode=input$lA3, Ea3_ode=input$Ea3, n3_ode=input$n3, m3_ode=input$m3, q3_ode=input$q3)
        vFunc=v9Func
      }      
      x=c(A=1-a0, B=a0, C=a0, D=a0)
      #x=c(A=0.99995, B=0.00003, C=0.000015, D=0.000005) #for compatibility with results of Netzsch software
    }
    
      #----------------------------------------------------Selection of the solver (see details in deSolve package)
      if (mD==1) mD_N = "euler"
      if (mD==2) mD_N = "rk2"
      if (mD==4) mD_N = "rk4"
      if (mD==5) mD_N = "rk45dp7"
      if (mD<6) out=as.data.frame(rk(x, ti, vFunc, parms,method = mD_N)) #, atol = 1e-12, rtol = 1e-10
      else #{
        if (mD==6)out=as.data.frame(lsoda(x, ti, vFunc, parms, atol = 1e-12, rtol = 1e-10))
      #  if(mD==7) out=as.data.frame(ode(x, ti, vFunc, parms, method = "lsode"))
      #  if(mD==8) out=as.data.frame(ode(x, ti, vFunc, parms, method = "daspk")) #atol = 1e-10, rtol = 1e-10)
      #  }
    if(KinModel<4) { #-----------------------------------Defines the output (important for two step models)
      al=out$a
      if (supY==T) { #-------------------------------------Whether additionally suppress NaN's and outliers
        al[is.na(al)]=2
        unity=0
        for (j in 2:length(al)) {
          if (unity==0) {
            if (al[j]<0) al[j]=0 
            else {
              if (al[j]>1) 
                if (al[j-1]>=0.8) { 
                  al[j]=1
                  unity=1 } } }
        else al[j]=1 } }
      C_A=1-al
      C_B=al
      C_C=0*al
      C_D=0*al }
    if(KinModel==4) {
      al=1-out$A
      al2=out$B
      if (supY==T) { #-------------------------------------Whether additionally suppress NaN's and outliers
        al[is.na(al)]=2
        unity=0
        for (j in 2:length(al)) {
          if (unity==0) {
            if (al[j]<0) al[j]=0 
            else {
              if (al[j]>1) 
                if (al[j-1]>=0.8) { 
                  al[j]=1
                  unity=1 } } }
          else al[j]=1 } 
        al2[is.na(al2)]=2
        unity=0
        for (j in 2:length(al2)) {
          if (unity==0) {
            if (al2[j]<0) al2[j]=0 
            else {
              if (al2[j]>1) 
                if (al2[j-1]>=0.8) { 
                  al2[j]=1
                  unity=1 } } }
          else al2[j]=1 } }
      C_A=1-al
      C_B=al2
      C_C=1-C_A-C_B
      C_D=0*C_A     }
    if(KinModel==5) {
      al1=1-out$A
      al2=out$B
      C_A=out$A
      C_B=out$B
      C_C=out$C
      if (supY==T) { #-------------------------------------Whether additionally suppress NaN's and outliers
        al1[is.na(al1)]=2
        unity=0
        for (j in 2:length(al1)) {
          if (unity==0) {
            if (al1[j]<0) al1[j]=0 
            else {
              if (al1[j]>1) 
                if (al1[j-1]>=0.8) { 
                  al1[j]=1
                  unity=1 } } }
          else al1[j]=1 } 
        al2[is.na(al2)]=2
        unity=0
        for (j in 2:length(al2)) {
          if (unity==0) {
            if (al2[j]<0) al2[j]=0 
            else {
              if (al2[j]>1) 
               if (al2[j-1]>=0.8) { 
                  al2[j]=1
                  unity=1 } } }
          else al2[j]=1         } 
        C_A=1-al1
        C_B=al2
        al=input$cd1*(1-C_A)+(1-input$cd1)*(1-C_A-C_B)
        C_C=1-C_A-C_B
      }
      al=input$cd1*(1-C_A)+(1-input$cd1)*(1-C_A-C_B)
      C_D=0*C_A     }
    if(KinModel==6) {
      al1=1-out$A
      al2=1-out$C
      if (supY==T) { #-------------------------------------Whether additionally suppress NaN's and outliers
        al1[is.na(al1)]=2
        unity=0
        for (j in 2:length(al1)) {
          if (unity==0) {
            if (al1[j]<0) al1[j]=0 
            else {
              if (al1[j]>1) 
                if (al1[j-1]>=0.8) { 
                  al1[j]=1
                  unity=1 } } }
          else al1[j]=1 } 
        al2[is.na(al2)]=2
        unity=0
        for (j in 2:length(al2)) {
          if (unity==0) {
            if (al2[j]<0) al2[j]=0 
            else {
              if (al2[j]>1) 
                if (al2[j-1]>=0.8) { 
                  al2[j]=1
                  unity=1 } } }
          else al2[j]=1 } }
      C_A=1-al1
      C_C=1-al2
      al=input$cd1*(1-C_A)+(1-input$cd1)*(1-C_C)
      C_B=1-C_A
      C_D=1-C_C   }
    
    #7
    
    #8
    
    if(KinModel==9) {
      #write.table(cbind(out$D,(1-out$A-out$B-out$C),(1-out$A_corr-out$B_corr-out$C)), file="test_v9.txt")
      C_A=ifelse(input$m1!=10, out$A, out$A_corr)
      C_B=ifelse(input$m1!=10, out$B, out$B_corr)
      C_C=out$C
      C_D=out$D
      al=input$cd1*(1-C_A)+(1-input$cd1)*(input$cr21*C_C+input$cr22*C_D)
      if (supY==T) { #-------------------------------------Whether additionally suppress NaN's and outliers
        al[is.na(al)]=2
        unity=0
        for (j in 2:length(al)) {
          if (unity==0) {
            if (al[j]<0) al[j]=0 
            else {
              if (al[j]>1) 
                if (al[j-1]>=0.9) { 
                  al[j]=1 
                  C_A[j]=0
                  C_B[j]=0
                  C_C[j]=C_C[j-1]
                  C_D[j]=C_D[j-1]
                  unity=1 } } }
          else {
            al[j]=1 
            C_A[j]=0
            C_B[j]=0
            C_C[j]=C_C[j-1]
            C_D[j]=C_D[j-1]
            } } }
  }
    
#    if (supY==T) { #-------------------------------------Whether additionally suppress NaN's and outliers
#      al[is.na(al)]=2
#      for (i in 2:Npts) {
#        if (al[i]>0) {
#          if (al[i]>1) { if (al[i-1]<=1) al[i]=1  }
#        }
#        else 
#          if (al[i-1]>0.9) al[i]=1
#          else al[i]=0 } }
    alph=suppressWarnings(splinefun(ti,al,method = "fmm"))
    dadt=alph(ti,deriv=1)
    
    list(ti,TempF(ti),al,dadt,C_A,C_B,C_C,C_D)
  }
  
  
  output$downloadData2 <- downloadHandler( #-----------generates file with thermal prediction data 
    filename = function() {
      paste("KinPrediction-", Sys.Date(), ".txt", sep="") },
    content = function(file) {
      resP=Prediction() 
      if (input$KinT == 2) {
        IsoKinPre=cbind(resP[[1]],resP[[2]],resP[[3]],resP[[4]])
        write.table(IsoKinPre, file, sep=" ", eol = "\r\n", col.names=c("Time[s]","Temperature[K]",
                                                                        "Conversion","Conversion rate[1/s]"),
                    row.names=FALSE, quote = FALSE) }
      else { 
        #(ti,TempF(ti),al,dadt,C_A,C_B,C_C,C_D)
        IsoKinPre=cbind(resP[[1]],resP[[2]],resP[[3]],resP[[4]],resP[[5]],resP[[6]], resP[[7]], resP[[8]])
        write.table(IsoKinPre, file, sep=" ", eol = "\r\n", col.names=c("Time[s]","Temperature[K]",
                                                                        "Conversion","Conversion rate[1/s]",
                                                                        "Concentration [A]","Concentration [B]",
                                                                        "Concentration [C]","Concentration [D]"),
                    row.names=FALSE, quote = FALSE) }
    })
  
  output$Plot1 <- renderPlot( height= function(){ if (input$Compact==T) 300
    else input$BaseHeight   },
  { #------------------------builds the plot with isoconversional activation energy and its spline fit
      inFile <- input$file1
      if (is.null(inFile)) return(NULL)
      res=Isoconv()
      par (fig=c(0,0.8,0,1.0),mar=c(5.1,4.9,4.1,0))  #old parameters (fig=c(0,0.8,0,1.0),mar=c(6.1,4.9,4.1,0),mgp=c(2.6,0.7,0),tcl=-0.5)
      EaStep=max((max(res[[2]])-min(res[[2]]))/10,50)
      EaStep2=min((max(res[[2]])-min(res[[2]]))/10,20)
      plot(res[[1]],res[[2]], 
           main = "Loaded isoconversional data",
           type="n", xlim=c(0, 1), ylim=c(min(res[[2]])-EaStep, max(res[[2]])+EaStep),
           xlab=NA, ylab=NA,pch=21, cex=2.5, bg="white", col="red",las=1, cex.main=input$TitleSize, font.main=input$TitleFont)
      title(xlab=expression('Conversion '~alpha),ylab=expression(italic('E')[a]*', kJ/mol'), # / kJ mol'^-1),
            cex.lab=input$AxisSize)
      arrows(res[[1]],res[[2]]-res[[3]],res[[1]],res[[2]]+res[[3]],length=0,angle=90,code=3,col="red",lwd=0.5)
      if(input$KinFT==1) EaFunc=splinefun(res[[1]],res[[2]],method = "fmm")
      if(input$KinFT==2) EaFunc=splinefun(res[[1]],res[[2]],method = "natural")
      if(input$KinFT==3) EaFunc=approxfun(res[[1]],res[[2]],method = "linear",rule=2)
      alpSeq=seq(0,1,by=0.02)
      lines(alpSeq,EaFunc(alpSeq), col = "red",lwd=1.5) 
      #lines(lowess(res[[1]],res[[2]],f=.1), col = "red",lwd=1.5) 
      points(res[[1]],res[[2]], pch=21, cex=2.5, bg="white", col="red")
      par (fig=c(0.85,1.0,0,1.0), mar=c(5.1,0,4.1,1.5),new=TRUE) #(fig=c(0.85,1.0,0,1.0), mar=c(6.1,0,4.1,1.5),new=TRUE)
      hs=hist(res[[2]], plot=FALSE, breaks=seq(min(res[[2]])-EaStep, max(res[[2]])+EaStep, EaStep2)) #breaks="FD"
      barplot(hs$counts, space=0, horiz=T, col = "red",border="white", axes=F)
    })
  
  output$Plot3 <- renderPlot( height= function(){ if (input$Compact==T) ifelse(((input$CompSolv==T)&&(input$KinMod>3)),250,350) 
                                                  else input$BaseHeight   },
  { #-----------------generates the main plot with prediction data
    
    if (input$KinT==2) {
      inFile1 <- input$file1
      if (is.null(inFile1)) return(NULL) }
    if (input$Tprg==3) {
      inFile2 <- input$file2
      if (is.null(inFile2)) return(NULL) 
      else TprgData=read.table(inFile2$datapath, skip=1) } # 1st - T[°C], 2nd row - t[s], 3rd (if exist) - alpha[-], 4th - da/dt [1/s]
    res=Prediction()
    par (fig=c(0,0.8,0,1.0),mar=c(5.1,5.5,3.1,0),las=1)#, mgp=c(3.7,0.7,0))
    if (input$showTP==3) { 
        if (input$PlotT==1) { 
            xPr=res[[2]]-273.15
            axisT='Temperature, °C'}
        else {
            xPr=res[[2]] 
            axisT='Temperature, K'} }
    else xPr=res[[1]]/60
    if (input$showY==1)  { 
        yPr=res[[3]]
        yAx=3
        yLab=expression('Conversion '~alpha) }
    if (input$showY==2) {
        yPr=res[[4]]
        yAx=4
        yLab='Conversion rate, 1/s' }
    plot(xPr,yPr,
         type="l",col = "black",lwd=2.5, ylab=NA,xlab=NA,cex=2.5)
    title(main = "Prediction for the specified (T,t) program",
      xlab=ifelse(input$showTP==3,axisT,'Time, min'), xpd=NA, 
      font.main=input$TitleFont, cex.lab=input$AxisSize,cex.main=input$TitleSize) #ylab=yLab,
    mtext(yLab,side=2,line=yAx,cex=input$AxisSize, las=3) 
    #text(par("usr")[1]-15,mean(par("usr")[3:4]), yLab, xpd=NA, srt = 90, cex=2)
    if ((input$Tprg==3) && (input$showAlTprg==T)) {
      if (input$showY==1) yTpr=TprgData[,3]
      if (input$showY==2) yTpr=TprgData[,4]  
      if (input$showTP==3) { 
        if (input$PlotT==1) xTpr=TprgData[,1]
        else xTpr=TprgData[,1]+273.15 }
      else xTpr=TprgData[,2]/60
      points(xTpr,yTpr)
      legend("topleft",legend=c("Prediction","From loaded file"),lty=c(1,NA), lwd=c(2.5,NA),
             pch=c(NA,1), inset = .05, cex=input$LegendSize) }
    if ((input$showTt==T)||(input$showTP==1)) {
      par(new = T)
      if (input$PlotT==1) { 
        xT=res[[2]]-273.15
        axisT='Temperature, °C'}
      else {
        xT=res[[2]] 
        axisT='Temperature, K'}
      plot(res[[1]]/60,xT, type="l", col="red", lwd=1.5, axes=F, xlab=NA, ylab=NA, cex=2.5)
      if (input$Tprg==3) points(TprgData[,2]/60,TprgData[,1], pch=1, cex=1.5, col="red") 
      mtext(axisT,side=4,col="red",line=4,cex=input$AxisSize, las=3) 
      #text(par("usr")[2]*1.11,mean(par("usr")[3:4]), axisT, srt = 90, xpd = NA, col="red",cex=2)
      axis(side = 4, col="red",col.axis="red",las=1) }
  })

  output$Plot5 <- renderPlot( height= function(){ if (input$Compact==T) ifelse(input$KinMod>3,200,310)
                                                  else input$BaseHeight  }, { #----compares the various solvers
    if ((input$KinT==2)||(input$CompSolv==F)||(input$ODEs==input$methodDiff)) return(NULL) 
    res=Prediction()
    par (fig=c(0,0.8,0,1.0),mar=c(5.1,5.5,3.1,0),las=1)
    if (input$showTP==3) { 
      if (input$PlotT==1) { 
        xPr=res[[2]]-273.15
        axisT='Temperature, °C'}
      else {
        xPr=res[[2]] 
        axisT='Temperature, K'} }
    else xPr=res[[1]]/60
    if (input$showY==1)  { 
      yPr=res[[3]]
      yAx=3
      yLab=expression('Conversion '~alpha) }
    if (input$showY==2) {
      yPr=res[[4]]
      yAx=4
      yLab='Conversion rate, 1/s' }
    plot(xPr,yPr, 
         type="l",col = rgb(0, 0, 0, 0.3),lwd=6,ylab=NA,xlab=NA)
    title(xlab=ifelse(input$showTP==3,axisT,'Time, min'), 
          main = "Comparsion of two different solvers for ODEs", font.main=input$TitleFont,
          cex.main=input$TitleSize, cex.lab=input$AxisSize)
    mtext(yLab,side=2,line=yAx,cex=input$AxisSize, las=3) 
    #legend!
    legend("topleft",legend=c("Main solver","Another selected solver"),lty=c(1,1), lwd=c(6,1),
           col=c(rgb(0, 0, 0, 0.3), "red"), bty = "n", inset = .02, cex=input$LegendSize)

    TempF=suppressWarnings(splinefun(res[[1]],res[[2]],method = "fmm")) 
    resF=NLR_Prediction(length(res[[1]]),res[[1]],TempF,input$ODEs,input$suppression)
    
#    if(input$KinMod==1) resF=funK1(input$lA1,input$Ea,input$n1,input$m1,
#               res[[2]],res[[1]],0,0,0,0,input$ODEs,input$q,input$a0,1,input$suppression)
    ##----repair others-----------------------------------------------------------------------------------------
#        if(input$KinMod==2) resF=funK2(input$lA1,input$Ea,input$reMod,
#                                       1/res[[2]]/8.31446,res[[1]],0,0,0,0,input$ODEs[i],input$a0,1)
#        #3
#        if(input$KinMod==4) {
#          q2=input$q2 #ifelse(input$q2f==T,input$q2,input$q)
#          resF=funK4(input$lA1,input$Ea,input$n1,input$m1,input$lA2,input$Ea2,input$n2,input$m2,
#                                       1/res[[2]]/8.31446,res[[1]],0,0,0,0,input$ODEs[i],input$q,q2,input$a0,1) }
#        if(input$KinMod==5) 
#          resF=funK5(input$lA1,input$Ea,input$n1,input$m1,input$lA2,input$Ea2,input$n2,input$m2,
#                                   input$cd1,0,1/res[[2]]/8.31446,res[[1]],0,0,0,0,0,input$methodDiff,input$q,input$q,input$a0,1)
        #6
        #7
    ##----repair others-----------------------------------------------------------------------------------------      
      al=resF[[3]]
      dadt= resF[[4]] 
      if (input$showY==1) yODE=al
      if (input$showY==2) yODE=dadt
      lines(xPr,yODE, col ="red",lwd=1)
          
      par(new = T)  #-----------------------------Draws the difference between two outputs
      Diff=yPr-yODE
      colorD="blue"
      plot(xPr,Diff, type="l", col=colorD, lwd=1.5, axes=F, xlab=NA, ylab=NA, cex=2.5)
      mtext("Difference",side=4,col=colorD,line=5,xpd=NA,cex=input$AxisSize, las=3)
      axis(side = 4, col=colorD,col.axis=colorD,las=1) 
  })
  
  output$Plot4 <- renderPlot( height= function(){ if (input$Compact==T) ifelse(input$CompSolv==T,250,300) 
                                                  else input$BaseHeight }, { #--generates plot with conversion for stages
    if (input$KinT==2) return(NULL) 
    if (input$KinMod<4) return(NULL)
    res=Prediction()
    #-------------------plot for v4 and v5 (?) add updateVal?
      par (fig=c(0,0.8,0,1.0),mar=c(5.1,5.5,3.1,0),las=1)
      concNames=c("A","B","C","D")
      if (input$showTP==3) { 
        if (input$PlotT==1) { 
          xPr=res[[2]]-273.15
          axisT='Temperature, °C'}
        else {
          xPr=res[[2]] 
          axisT='Temperature, K'} }
      else xPr=res[[1]]/60
      LeftNumb=length(input$ConcLeft)
      if (LeftNumb==0) return(NULL)
        ltys=numeric(LeftNumb)
        lwds=numeric(LeftNumb)
        cols=numeric(LeftNumb)
        minYs=numeric(LeftNumb)
        maxYs=numeric(LeftNumb)
        ltys[1]=1
        cols[1]=2
        lwds[1]=2.5
        for (i in 1:LeftNumb) {
          minYs[i]=min(res[[5+as.numeric(input$ConcLeft[i])]])
          maxYs[i]=max(res[[5+as.numeric(input$ConcLeft[i])]]) }
        yLeftPr=res[[5+as.numeric(input$ConcLeft[1])]]
        axisYLeft=paste0('Concentration of ', paste(unlist(concNames[as.numeric(input$ConcLeft)+1]), collapse=', ')) 
#--------------The order in Prediction() output: time, temperature, conversion, da/dt, A, B, C, (D)
        plot(xPr,yLeftPr, ylim = c(min(minYs), max(maxYs)), 
           type="l",col = cols[1],lwd=2.5, ylab=NA,xlab=NA,cex=2.5)
        title(xlab=ifelse(input$showTP==3,axisT,'Time, min'), ylab=axisYLeft, xpd=NA, 
              main = "Predicted concentration of components", 
              font.main=input$TitleFont,cex.main=input$TitleSize, cex.lab=input$AxisSize)
      if (LeftNumb>1) {
        for (i in 2:LeftNumb) {
          ltys[i]=1
          cols[i]=i+1
          lwds[i]=2.5
          lines(xPr,res[[5+as.numeric(input$ConcLeft[i])]], col=cols[i],lwd=2.5) }
      }
      RightNumb=length(input$ConcRight)
      if (RightNumb>0) {
        par(new = T)  #-----------------------------Draws the concentration data over the secondary axis
        ltysR=numeric(RightNumb)
        lwdsR=numeric(RightNumb)
        colsR=numeric(RightNumb)
        minYsR=numeric(RightNumb)
        maxYsR=numeric(RightNumb)
        ltysR[1]=1
        colsR[1]=LeftNumb+2
        lwdsR[1]=2.5
        for (i in 1:RightNumb) {
          minYsR[i]=min(res[[5+as.numeric(input$ConcRight[i])]])
          maxYsR[i]=max(res[[5+as.numeric(input$ConcRight[i])]]) }
        yRightPr=res[[5+as.numeric(input$ConcRight[1])]]
        axisYRight=paste0('Concentration of ', paste(unlist(concNames[as.numeric(input$ConcRight)+1]), collapse=', ')) 
        plot(xPr,yRightPr, axes=F, type="l",col =colsR[1],lwd=2.5, ylab=NA,xlab=NA,cex=2.5)
        mtext(axisYRight,side=4,col=1,line=4,xpd=NA,cex=input$AxisSize, las=3)
        axis(side = 4, col=1,col.axis=1,las=1, ylim = c(min(minYsR), max(maxYsR)),cex.main=1.3) 
        if (RightNumb>1) {
          for (i in 2:RightNumb) {
            ltysR[i]=1
            colsR[i]=LeftNumb+i+1
            lwdsR[i]=2.5
            lines(xPr,res[[5+as.numeric(input$ConcRight[i])]], col=colsR[i],lwd=2.5) }
        }
        colU=c(cols,colsR)
        legU=c(concNames[as.numeric(input$ConcLeft)+1],concNames[as.numeric(input$ConcRight)+1])
       # legend("topleft",legend=legU, lty=c(ltys,ltysR), lwd=c(lwds,lwdsR), col=colU, bty = "n", inset = .02)
        }
        else {
          legU=concNames[as.numeric(input$ConcLeft)+1]
          colU=cols }
          #legend("topleft",legend=legU,lty=ltys, lwd=lwds, col=colU, bty = "n", inset = .02) }
        if(input$KinMod==4) {
          labelModelIni="Model v4: A -> B, A -> C"
          labelModel=expression("Model v4: "*phantom("A")*" -> "*phantom("B")*", "*phantom("A")*" -> "*phantom("C"))
          labelModelTr=expression(phantom("Model v4: ")*"A"*phantom(" -> ")*"B"*phantom(", ")*"A"*phantom(" -> ")*"C")
          inistr="phantom(\"Model v4: A -> B, A -> C\")" }
        if(input$KinMod==5) {
          labelModelIni="Model v5: A -> B -> C"
          labelModel=expression("Model v5: "*phantom("A")*" -> "*phantom("B")*" -> "*phantom("C"))
          labelModelTr=expression(phantom("Model v5: ")*"A"*phantom(" -> ")*"B"*phantom(" -> ")*"C")
          inistr="phantom(\"Model v5: A -> B -> C\")" }
        if(input$KinMod==6) {
          labelModelIni=expression("Model v6: A -> B, C -> D")
          labelModel=expression("Model v6: "*phantom("A")*" -> "*phantom("B")*", "*phantom("C")*" -> "*phantom("D"))
          labelModelTr=expression(phantom("Model v6: ")*"A"*phantom(" -> ")*"B"*phantom(", ")*"C"*phantom(" -> ")*"D")
          inistr="phantom(\"Model v6: A -> B, C -> D\")" }
      
      #7
      
      #8
      
      if(input$KinMod==9) {
        labelModelIni=expression("Model v9: A -> B -> C || D")
        labelModel=expression("Model v9: "*phantom("A")*" -> "*phantom("B")*" -> "*phantom("C")*" || "*phantom("D"))
        labelModelTr=expression(phantom("Model v6: ")*"A"*phantom(" -> ")*"B"*phantom(" -> ")*"C"*phantom(" || ")*"D")
        inistr="phantom(\"Model v6: A -> B -> C || D\")" }
        usr=par("usr")
        xCoord=(usr[2]-usr[1])*1/40+usr[1]
        text(xCoord,(usr[3]-usr[4])*0.15+usr[4], labelModel, pos=4, cex=input$LegendSize, col=1) 
        text(xCoord,(usr[3]-usr[4])*0.15+usr[4], labelModelTr, pos=4, cex=input$LegendSize, col=1) 
        for (i in 1:(LeftNumb+RightNumb)) {
          labelI=str2expression(gsub(legU[i], paste0("\")*'",legU[i],"'*phantom(\""), inistr))  #str2expression(gsub(legU[i], paste0("\")*bold('",legU[i],"')*phantom(\""), inistr))
          text(xCoord,(usr[3]-usr[4])*0.15+usr[4], labelI, pos=4, cex=input$LegendSize, col=colU[i])
        }
  })

  reactionModel=function(Rtype,alpha) {  #idealized reaction model types. Five group of models from:
    #[S. Vyazovkin, A. K. Burnham, J. M. Criado, et al, ICTAC Kinetics Committee recommendations for performing kinetic computations on thermal analysis data, Thermochimica Acta, Vol. 520 (1–2), 2011, pp. 1–19; 
    #P.E. Sanchez-Jimenez, L.A. Perez-Maqueda, A.Perejon, J.M. Criado, Constant rate thermal analysis for thermal stability studies of polymers, Polymer Degradation and Stability, 96, 2011, 974-981]:
    # (1) nucleation models - power law (P2,P3,P4,P23), Avrami-Erofeev (A2,A3,A4) and Prout-Tompkins (B1)
    # (2) geometrical contraction models - contracting area (R2) and sphere (R3)
    # (3) diffusion models (D1,D2,D3,D4)
    # (4) reaction-order models (F0,F1,F2,F3)
    # (5) random scission model (L2)
    if (Rtype=="F0") res=1-alpha+alpha
    if (Rtype=="F1") res=(1-alpha)
    if (Rtype=="F2") res=(1-alpha)^2
    if (Rtype=="F3") res=(1-alpha)^3
    if (Rtype=="R2") res=2*(1-alpha)^0.5
    if (Rtype=="R3") res=3*(1-alpha)^(2/3)
    if (Rtype=="A2") res=suppressWarnings(2*(1-alpha)*(-log(1-alpha))^0.5)
    if (Rtype=="A3") res=suppressWarnings(3*(1-alpha)*(-log(1-alpha))^(2/3))
    if (Rtype=="A4") res=suppressWarnings(4*(1-alpha)*(-log(1-alpha))^(3/4))
    if (Rtype=="P2") res=2*alpha^0.5
    if (Rtype=="P3") res=3*alpha^(2/3)
    if (Rtype=="P4") res=4*alpha^(3/4)
    if (Rtype=="P23") res=2/3*alpha^(-0.5)
    if (Rtype=="D1") res=0.5/alpha
    if (Rtype=="D2") res=suppressWarnings(-1/log(1-alpha))
    if (Rtype=="D3") res=1.5*(1-alpha)^(2/3)/(1-(1-alpha)^(1/3))
    if (Rtype=="D4") res=1.5/((1-alpha)^(-1/3)-1)
    if (Rtype=="L2") res=2*(alpha^0.5-alpha)
    if (Rtype=="B1") res=(1-alpha)*alpha
    res }
  

  #-------------------------FUNKY TIME (Old Model fitting functions start here)-------------------------------------------------------------------
  
  DMM58 = function (lA1,Ea1,mu,lA2,Ea2,oT,al) { #Dubovitskiy-Manelis-Merzhanov-1958 function
    if (mu==1) mpart=0
    else {
      if (mu==0) mpart=al*(1-al)
      else mpart=(1-mu)*al*(1-al)/(1-mu*al) }
    exp(lA1-Ea1*1000*oT)*(1-al)+exp(lA2-Ea2*1000*oT)*mpart
  }
  
  ePT = function (lA,Ea,n,m,q,oT,al) { #ePT (extended Prout-Tompkins or truncated Sestack-Berggren) function
    if (m==0) m_part=1
    else { 
      if (q==1) m_part=al^m
      else m_part=(1-q*(1-al))^m }
    if (n==0) n_part=1
    else n_part=(1-al)^n
    exp(lA-Ea*1000*oT)*n_part*m_part
  }
  
  funK1=function(lA1,Ea1,n1,m1,fT,ft,fa,nr,Dtype,WghDA,method,q,amin,amax,supY) {  #single reaction in flexible ePT form
    alpLim=amin
    Nal=length(ft)
    if (method==1) {     #------------------------crude Euler
      alpI=double(Nal)
      dadtI=double(Nal)
      alpI[1]=alpLim
      dadtI[1]=0
      fT1=1/fT/8.31446
        for (i in 2:Nal) {
          h=ft[i]-ft[i-1]
          alpI[i]=alpI[i-1]+h*dadtI[i-1]
          dadtI[i]=ePT(lA1,Ea1,n1,m1,q,fT1[i],alpI[i]) 
          if (supY==T) { #-------------------------------------Whether additionally suppress NaN's and outliers
            if ((is.na(alpI[i]))&&(alpI[i-1]>=amax-0.1)) {
              alpI[i]=1
              dadtI[i]=0  }
            else {
              if (is.na(alpI[i])) {
                alpI[i]=2
                dadtI[i]=0 }
              if (alpI[i]>=amax) {
                alpI[i]=1
                dadtI[i]=0 } } #------------------------------check it(!)
          } } }
    else {
    TempF=suppressWarnings(splinefun(ft,fT,method = "fmm")) 
    v1Func <- function(t, x, parms) {
      with(as.list(parms), {
        dx <- exp(lA1_ode-Ea_ode*1000/8.31446/TempF(t))*((1-x)^n1_ode)*(1-q_ode*(1-x))^m1_ode
        list(dx)
      }) }
    x=c(a = alpLim)
    parms=c(lA1_ode=lA1, Ea_ode=Ea1, n1_ode=n1, m1_ode=m1, q_ode=q)
    
   # if (mD==1) mD_N = "euler" #---------------------------Selection of the solver (see details in deSolve package)
    if (method==2) mD_N = "rk2"
    if (method==4) mD_N = "rk4"
    if (method==5) mD_N = "rk45dp7"
    if (method<6) out=as.data.frame(rk(x, ft, v1Func, parms,method = mD_N))
    else if (method==6)out=as.data.frame(lsoda(x, ft, v1Func, parms))
    alpI=out$a
    if (supY==T) { #-------------------------------------Whether additionally suppress NaN's and outliers
      alpI[is.na(alpI)]=2
      for (i in 2:length(alpI)) {
        if (alpI[i]>0) {
          if (alpI[i]>amax) { if (alpI[i-1]<=amax) alpI[i]=amax  } }
        else 
          if (alpI[i-1]>amax-0.1) alpI[i]=1
          else alpI[i]=0 } }
    alph=suppressWarnings(splinefun(ft,alpI,method = "fmm"))  #----computes da/dt
    dadtI=alph(ft,deriv=1)
    }
    res=list(ft,fT,alpI,dadtI,0*alpI,0*alpI,0*alpI,0*alpI)
  }
  
  funK2=function(lA1,Ea1,ReType,fT1,ft,fa,nr,Dtype,WghDA,method,amin,amax) {  #single reaction in a theoretical form
    alpLim=amin
    Nal=length(ft)
    alpI=double(Nal)
    dadtI=double(Nal)
    alpI[1]=alpLim
    dadtI[1]=0
    for (i in 2:Nal) {
        h=ft[i]-ft[i-1]
        if (method==1) {      # crude (??) Euler method
          alpI[i]=alpI[i-1]+h*dadtI[i-1]
          dadtI[i]=exp(lA1-Ea1*1000*fT1[i])*reactionModel(ReType,alpI[i]) }
        if (method==2) {      # Improved Euler method, a second order!
          approx_al=alpI[i-1]+h*dadtI[i-1]
          alpI[i]=alpI[i-1]+h*(dadtI[i-1]+exp(lA1-Ea1*1000*fT1[i])*reactionModel(ReType,approx_al))/2
          dadtI[i]=exp(lA1-Ea1*1000*fT1[i])*reactionModel(ReType,alpI[i])  }
        if (method==4) {      # fourth-order Runge-Kutta method!
          fTh2=2/(1/fT1[i]+1/fT1[i-1]) #fT=1/Te/8.314
          u2=alpI[i-1]+h*dadtI[i-1]/2
          k2=exp(lA1-Ea1*1000*fTh2)*reactionModel(ReType,u2) 
          u3=alpI[i-1]+h*k2/2            #alpI[i-1]+h*tSB(lA1,Ea1,n1,m1,fTh2,u2)/2
          k3=exp(lA1-Ea1*1000*fTh2)*reactionModel(ReType,u3)
          u4=alpI[i-1]+h*k3              #alpI[i-1]+h*tSB(lA1,Ea1,n1,m1,fTh2,u3)
          k4=exp(lA1-Ea1*1000*fT1[i])*reactionModel(ReType,u4) 
          alpI[i]=alpI[i-1]+h*(dadtI[i-1]+2*k2+2*k3+k4)/6
          dadtI[i]=exp(lA1-Ea1*1000*fT1[i])*reactionModel(ReType,alpI[i]) }
        
        #if(is.na(alpI[i])) cat(lA1,Ea1,n1,m1,q,i,alpI[i-1],dadtI[i-1],alpI[i],fT1[i],dadtI[i],"\n",file="aaaa1.txt", append = T)
        if ((is.na(alpI[i]))&&(alpI[i-1]>=amax-0.1)) {
          alpI[i]=1
          dadtI[i]=0  }
        else {
          if (is.na(alpI[i])) {
            alpI[i]=2
            dadtI[i]=0 }
          if (alpI[i]>=amax) {
            alpI[i]=1
            dadtI[i]=0 } } 
    }
    res=cbind(alpI,dadtI)
  }
  
  
  funK4=function(lA1,Ea1,n1,m1,lA2,Ea2,n2,m2,fT1,ft,fa,nr,Dtype,WghDA,method,q1,q2,amin,amax) {  #two independent reactions
    alpLim=amin
    Nal=length(ft)
    alpI=numeric(Nal)
    dadtI=numeric(Nal)
    alpI[1]=alpLim
    dadtI[1]=0
    if (Ea2==0) Ea2=Ea1
    for (i in 2:Nal) {
        h=ft[i]-ft[i-1]
        if (method==1) {      # crude (??) Euler method
          alpI[i]=alpI[i-1]+h*dadtI[i-1]
          dadtI[i]=tSB(lA1,Ea1,n1,m1,q1,fT1[i],alpI[i])+tSB(lA2,Ea2,n2,m2,q2,fT1[i],alpI[i]) 
        }
        if (method==2) {      # Improved Euler method, a second order!
          approx_al=alpI[i-1]+h*dadtI[i-1]
          alpI[i]=alpI[i-1]+h*(dadtI[i-1]+tSB(lA1,Ea1,n1,m1,q1,fT1[i],approx_al)+tSB(lA2,Ea2,n2,m2,q2,fT1[i],approx_al))/2
          dadtI[i]=tSB(lA1,Ea1,n1,m1,q1,fT1[i],alpI[i])+tSB(lA2,Ea2,n2,m2,q2,fT1[i],alpI[i])  }
        if ((is.na(alpI[i]))&&(alpI[i-1]>=amax-0.1)) { #-check for NaN's
          alpI[i]=1
          dadtI[i]=0  }
        else {
          if (is.na(alpI[i])) {
            alpI[i]=2
            dadtI[i]=0 }
          if (alpI[i]>=amax) {
            alpI[i]=1
            dadtI[i]=0 } }
    }
    res=cbind(alpI,dadtI)
  }
  
  funK5=function(lA1,Ea1,n1,m1,lA2,Ea2,n2,m2,cd1,cd2,fT1,ft,fa,flag,nr,Dtype,WghDA,method,q1,q2,amin,amax) {  #two consecutive reactions
    alpLim=amin
    Nal=length(ft)
    alpI=double(Nal)
    dadtI=double(Nal)
    alpI1=double(Nal)
    dadtI1=double(Nal)
    alpI2=double(Nal)
    dadtI2=double(Nal)
    alpI[1]=alpLim
    dadtI[1]=0
    alpI1[1]=alpLim
    dadtI1[1]=0
    alpI2[1]=alpLim
    dadtI2[1]=0
    if (Ea2==0) Ea2=Ea1
    for (i in 2:Nal) {
        h=ft[i]-ft[i-1]
        if (method==1)       # crude Euler method, delete after all(?)
          alpI1[i]=alpI1[i-1]+h*dadtI1[i-1]
        if (method==2) {      # Improved Euler method, a second order!
          approx_al1=alpI1[i-1]+h*dadtI1[i-1]
          alpI1[i]=alpI1[i-1]+h*(dadtI1[i-1]+tSB(lA1,Ea1,n1,m1,q1,fT1[i],approx_al1))/2 }
        
          if (alpI1[i]>=amax) {
            alpI1[i]=1
            dadtI1[i]=0 } 
          else dadtI1[i]=tSB(lA1,Ea1,n1,m1,q1,fT1[i],alpI1[i]) 
          
        if (method==1)  alpI2[i]=alpI2[i-1]+h*dadtI2[i-1]
        if (method==2) {
          approx_al2=alpI2[i-1]+h*dadtI2[i-1]
          alpI2[i]=alpI2[i-1]+h*(dadtI2[i-1]+tSB(lA2,Ea2,n2,m2,q2,fT1[i],1-alpI1[i]+approx_al2))/2 }
        
        if (is.na(alpI1[i])) alpI1[i]=0
        if (is.na(alpI2[i])) alpI2[i]=0
        
          if (alpI2[i]>=amax) {
            alpI2[i]=1
            dadtI2[i]=0 } 
          else {
            if (alpI2[i]>alpI1[i]) {
              alpI2[i]=alpI1[i] #ifelse(alpI2[i-1]>alpI1[i],alpI2[i-1],alpI1[i]) 
              dadtI2[i]=0 }
            else dadtI2[i]=tSB(lA2,Ea2,n2,m2,q2,fT1[i],1-alpI1[i]+alpI2[i]) }
        
        alpI[i]=alpI[i-1]+h*dadtI[i-1] 
        if (alpI[i]>=amax) {
          alpI[i]=1
          dadtI[i]=0 } 
        else {
          if (alpI[i]>alpI1[i]) { #----to prevent error when the second stage kinetics is exceptionally high(!???)
            alpI[i]=alpI1[i]
            dadtI[i]=0 }
          else dadtI[i]=cd1*dadtI1[i]+(1-cd1)*dadtI2[i]
          }
    }
    res=cbind(alpI,dadtI,alpI1,dadtI1,alpI2,dadtI2)
  }

  
})