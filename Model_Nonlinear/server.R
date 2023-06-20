#THINKS thermokinetic software, version #01.07.20 Website: thinks.chemphys.ru. To cite: 10.3390/molecules24122298
#This routine is used to perform the model fitting (specifically, the nonlinear regression)
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
library(shinyjs)
library(minpack.lm)
library(nlmrt)
library(deSolve)
library(devEMF)
#library(deTestSet)
#library(nlsr)

options(shiny.maxRequestSize=2*1024^2)  #change the file-size limit here----in the server there is 10 MB limit

shinyServer(function(input, output,session) {

  session$onSessionEnded(function() {
    observe({ 
      stopApp() })
  })
  
  output$cond <- reactive({
    is.null(input$file1)
  })
  outputOptions(output, "cond", suspendWhenHidden = FALSE)
  
  observeEvent(input$file1, {
    updateTabsetPanel(session, "MainTabset",
                      selected = "panel1")
  }, once=TRUE)
  
  observe({ #---------------generates the initial guesses to some of the kinetic parameters
    vallnA2=17.4
    valEa2=96
#    if (input$taskType==1) {
#      vallnA2=13.6
#      valEa2=81 }
#    if ((input$taskType==2)||(input$taskType==3)) {
#      vallnA2=3.3
#      valEa2=86.6 }
#    if (input$taskType>=4) {
#      vallnA2=12
#      valEa2=110 }
#    updateSliderInput(session, "lnA2", value = vallnA2)
#    updateSliderInput(session, "Ea2", value = valEa2)
#    if (input$taskType>=7) updateSliderInput(session, "n2", value = 2, min = 0.5, max = 4) 
    if (input$m1fixed==10) updateSliderInput(session, "n1", value = 2, min = 0.1, max = 6) 
    if (input$m2fixed==10) updateSliderInput(session, "n2", value = 2, min = 0.1, max = 6) 
    
#    if (input$taskType!=1) {
#        if (input$taskType==2) MethodsDiff=c("Euler"=1,"Improved Euler"=2,"Runge-Kutta 4th"=4)
#        if (input$taskType==5) MethodsDiff=c("Euler"=1)
#        if ((input$taskType>2)&&(input$taskType!=5)) MethodsDiff=c("Euler"=1,"Improved Euler"=2)
#        updateRadioButtons(session, "methodDiff", choices = MethodsDiff, inline =T) }
  })

  observeEvent(input$file1, {
    if ((is.null(input$file1))||(!input$go)) return(NULL)
    updateTabsetPanel(session, "MainTabset",
                      selected = "panel1")
  }, once=TRUE)
  
  output$setButton <- renderUI({
    if ((is.null(input$file1))||(!input$go)) return(NULL)
    if (reggr()[[3]][1] == "ERROR") return(NULL)
    else actionButton("setIni", "Set the current results", class="btn btn-secondary")
  })
  
  observeEvent(input$setIni, { #set the initial guess values using the NLR results
    res=isolate(reggr())
    if (res[[3]][1] != "ERROR") {
      coefs=coef(res[[1]])
      names=names(coefs)
      coefs=unname(coefs)
      for (i in 1:length(coefs)) {
        updateSliderInput(session, names[i], value=coefs[i])
      }
    }
  })  
  
al_min_global=0.00005
al_max_global=0.99995
Kt=100 #increase of the reaction rate after melting
MyNull=1E-11
MyUnity=1-MyNull
MyNull_for_Log=3E-18

ePT = function (lA,Ea,n,m,q,oT,al) { #-denoted below as ePT, the flexible equation (called also the extended Prout-Tompkins or truncated Sestak–Berggren function)
  if (m==0) m_part=1
  else { 
    if (q==1) m_part=al^m
    else m_part=(1-q*(1-al))^m }
  if (n==0) n_part=1
  else n_part=(1-al)^n
  exp(lA-Ea*1000*oT)*n_part*m_part
  }

#single reaction in flexible ePT form, only simple methods are here, the next function is applied for higher-order solvers
funK1s=function(lnA1,Ea1,n1,m1,fT1,ft,fa,nr,Dtype,WghDA,method,q,amin,amax,A_trigg) {  
  alpLim=amin
  Nal=length(fa)
  alpI=numeric(Nal)
  dadtI=numeric(Nal)
  alpI[1]=alpLim
  dadtI[1]=0
  for (i in 2:Nal) {
    if ((nr[i]==nr[i-1]+1)||((nr[i-1]>2)&&(nr[i]==1))) {  #-to switch between runs in row (conversion=1 to 0)
      alpI[i]=alpLim
      dadtI[i]=0 }
    else {
      h=ft[i]-ft[i-1]
      if (method==1) {      # crude Euler method
        alpI[i]=alpI[i-1]+h*dadtI[i-1]
        dadtI[i]=ePT(lnA1+A_trigg[i],Ea1,n1,m1,q,fT1[i],alpI[i]) }
      if (method==2) {      # Improved Euler method, a second order!
        approx_al=alpI[i-1]+h*dadtI[i-1]
        alpI[i]=alpI[i-1]+h*(dadtI[i-1]+ePT(lnA1+A_trigg[i],Ea1,n1,m1,q,fT1[i],approx_al))/2
        dadtI[i]=ePT(lnA1+A_trigg[i],Ea1,n1,m1,q,fT1[i],alpI[i]) }
      if (method==4) {      # fourth-order Runge-Kutta method!
        fTh2=2/(1/fT1[i]+1/fT1[i-1]) #fT=1/Te/8.314
        u2=alpI[i-1]+h*dadtI[i-1]/2
        k2=ePT(lnA1+A_trigg[i],Ea1,n1,m1,q,fTh2,u2)
        u3=alpI[i-1]+h*k2/2            #alpI[i-1]+h*ePT(lnA1,Ea1,n1,m1,fTh2,u2)/2
        k3=ePT(lnA1+A_trigg[i],Ea1,n1,m1,q,fTh2,u3)
        u4=alpI[i-1]+h*k3              #alpI[i-1]+h*ePT(lnA1,Ea1,n1,m1,fTh2,u3)
        k4=ePT(lnA1+A_trigg[i],Ea1,n1,m1,q,fT1[i],u4)
        alpI[i]=alpI[i-1]+h*(dadtI[i-1]+2*k2+2*k3+k4)/6
        dadtI[i]=ePT(lnA1+A_trigg[i],Ea1,n1,m1,q,fT1[i],alpI[i]) }
      #if(is.na(alpI[i])) cat(lnA1,Ea1,n1,m1,q,i,alpI[i-1],dadtI[i-1],alpI[i],fT1[i],dadtI[i],"\n",file="aaaa1.txt", append = T)
      if ((is.na(alpI[i]))&&(alpI[i-1]>=amax-0.1)) {
        alpI[i]=1
        dadtI[i]=0  }
      else {
        if (is.na(alpI[i])) {
          alpI[i]=2
          dadtI[i]=0 }
        if (alpI[i]>=amax) {
          alpI[i]=1
          dadtI[i]=0 } } } 
  }
  if (Dtype==3) res=dadtI*WghDA
  else  if(Dtype==1) res=dadtI
  else res=alpI
  res
}
#single reaction in flexible ePT form, function below uses the sunctions from library deSolve, 
  funK1=function(lnA1,Ea1,n1,m1,fT1,ft,fa,nr,Dtype,WghDA,method,q,amin,amax,Nruns,TempF,supY,ti_trigg) {  
    alpLim=amin
    for (i in 1:Nruns) {
      if (m1!=10) {
        v1Func <- function(t, x, parms) {
          with(as.list(parms), {
            x=ifelse(x<MyUnity,x,MyUnity)
            lnAt=ifelse(t>ti_trigg[i], lnA1_ode, lnA1_ode-log(Kt))
            dx <- exp(lnAt-Ea_ode*1000/8.31446/TempF[[i]](t))*((1-x)^n1_ode)*(1-q_ode*(1-x))^m1_ode
            list(dx)
          }) } 
        parms=c(lnA1_ode=lnA1, Ea_ode=Ea1, n1_ode=n1, m1_ode=m1, q_ode=q)
        }
      else {
        v1Func <- function(t, x, parms) { #----reaction model in extended KJMAE form
          with(as.list(parms), {
            x=ifelse(x<MyUnity,x,MyUnity)
            lnAt=ifelse(t>ti_trigg[i], lnA1_ode, lnA1_ode-log(Kt))
            dx <- exp(lnAt-Ea_ode*1000/8.31446/TempF[[i]](t))*n1_ode*(1-x)*(-log(1-x))^((n1_ode-1)/n1_ode)
            list(dx)
          }) } 
        parms=c(lnA1_ode=lnA1, Ea_ode=Ea1, n1_ode=n1)
        }
        x=c(a = alpLim)
        tI=ft[nr==i]
        if (method==1) mD_N = "euler" #---------------------------Selection of the solver (see details in deSolve package)
        if (method==2) mD_N = "rk2"
        if (method==4) mD_N = "rk4"
        if ((method==2)||(method==4)) out=as.data.frame(rk(x, tI, v1Func, parms,method = mD_N))
        if (method==5) out=as.data.frame(ode(x, tI, v1Func, parms, method = rkMethod("rk45ck", densetype = NULL, nknots = 0)))
        if (method==6) out=as.data.frame(ode(x, tI, v1Func, parms, method = rkMethod("rk45dp7", densetype = NULL, nknots = 0)))
        if (method==9) out=as.data.frame(lsoda(x, tI, v1Func, parms, atol = 1e-10, rtol = 1e-10))
        if (method==7) out=as.data.frame(ode(x, tI, v1Func, parms, method = rkMethod("rk45ck")))
        if (method==8) out=as.data.frame(ode(x, tI, v1Func, parms, method = rkMethod("rk45dp7")))
        alp=out$a
        if (supY==T) { #-------------------------------------Whether additionally suppress NaN's and outliers
          alp[is.na(alp)]=2
          unity=0
          for (j in 2:length(alp)) {
            if (unity==0) {
              if (alp[j]<0) alp[j]=0 
              else {
                if (alp[j]>amax) 
                  if (alp[j-1]>=amax-0.2) { 
                    alp[j]=amax
                    unity=1 }
              } }
            else alp[j]=amax  } }
        if (Dtype!=2) {
            alph=suppressWarnings(splinefun(tI,alp,method = "fmm"))  #----computes da/dt
            dadt=alph(tI,deriv=1) 
            if(i==1) dadtI=dadt
            else dadtI=c(dadtI,dadt) }
        if(i==1) alpI=alp
        else alpI=c(alpI,alp)
    }
    if (Dtype==3) res=dadtI*WghDA
    else  if(Dtype==1) res=dadtI
    else res=alpI
    res
  }
    
  reactionModel=function(Rtype,alpha) {  #idealized reaction model types. Five group of models from:
    #S. Vyazovkin, A. K. Burnham, J. M. Criado, et al, ICTAC Kinetics Committee recommendations for performing kinetic computations on thermal analysis data, Thermochimica Acta, Vol. 520 (1–2), 2011, pp. 1–19
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

#single reaction in a theoretical form  
  funK2s=function(lnA1,Ea1,ReType,fT1,ft,fa,nr,Dtype,WghDA,method,amin,amax,A_trigg) {  
    alpLim=amin
    Nal=length(fa)
    alpI=numeric(Nal)
    dadtI=numeric(Nal)
    alpI[1]=alpLim
    dadtI[1]=0
    for (i in 2:Nal) {
      if ((nr[i]==nr[i-1]+1)||((nr[i-1]>2)&&(nr[i]==1))) {  #-to switch between runs in row (conversion=1 to 0)
        alpI[i]=alpLim
        dadtI[i]=0 }
      else {
        h=ft[i]-ft[i-1]
        if (method==1) {      # crude (??) Euler method
          alpI[i]=alpI[i-1]+h*dadtI[i-1]
          dadtI[i]=exp(lnA1+A_trigg[i]-Ea1*1000*fT1[i])*reactionModel(ReType,alpI[i]) }
        if (method==2) {      # Improved Euler method, a second order!
          approx_al=alpI[i-1]+h*dadtI[i-1]
          alpI[i]=alpI[i-1]+h*(dadtI[i-1]+exp(lnA1+A_trigg[i]-Ea1*1000*fT1[i])*reactionModel(ReType,approx_al))/2
          dadtI[i]=exp(lnA1+A_trigg[i]-Ea1*1000*fT1[i])*reactionModel(ReType,alpI[i])  }
        if (method==4) {      # fourth-order Runge-Kutta method!
          fTh2=2/(1/fT1[i]+1/fT1[i-1]) #fT=1/Te/8.314
          u2=alpI[i-1]+h*dadtI[i-1]/2
          k2=exp(lnA1+A_trigg[i]-Ea1*1000*fTh2)*reactionModel(ReType,u2) 
          u3=alpI[i-1]+h*k2/2            #alpI[i-1]+h*ePT(lnA1,Ea1,n1,m1,fTh2,u2)/2
          k3=exp(lnA1+A_trigg[i]-Ea1*1000*fTh2)*reactionModel(ReType,u3)
          u4=alpI[i-1]+h*k3              #alpI[i-1]+h*ePT(lnA1,Ea1,n1,m1,fTh2,u3)
          k4=exp(lnA1+A_trigg[i]-Ea1*1000*fT1[i])*reactionModel(ReType,u4) 
          alpI[i]=alpI[i-1]+h*(dadtI[i-1]+2*k2+2*k3+k4)/6
          dadtI[i]=exp(lnA1+A_trigg[i]-Ea1*1000*fT1[i])*reactionModel(ReType,alpI[i]) }

        #if(is.na(alpI[i])) cat(lnA1,Ea1,n1,m1,q,i,alpI[i-1],dadtI[i-1],alpI[i],fT1[i],dadtI[i],"\n",file="aaaa1.txt", append = T)
        if ((is.na(alpI[i]))&&(alpI[i-1]>=amax-0.1)) {
          alpI[i]=1
          dadtI[i]=0  }
        else {
          if (is.na(alpI[i])) {
            alpI[i]=2
            dadtI[i]=0 }
          if (alpI[i]>=amax) {
            alpI[i]=1
            dadtI[i]=0 } } } 
    }
    if (Dtype==3) res=dadtI*WghDA
    else  if(Dtype==1) res=dadtI
    else res=alpI
    res
  }

  funK2=function(lnA1,Ea1,ReType,fT1,ft,fa,nr,Dtype,WghDA,method,amin,amax,Nruns,TempF,supY,ti_trigg) {  
    alpLim=amin
    for (i in 1:Nruns) {
      rMod=ReType
      v2Func <- function(t, x, parms) {
        with(as.list(parms), {
          lnAt=ifelse(t>ti_trigg[i], lnA1_ode, lnA1_ode-log(Kt))
          dx <- exp(lnAt-Ea_ode*1000/8.31446/TempF[[i]](t))*reactionModel(rMod,x)
          list(dx)
        }) }
      x=c(a = alpLim)
      parms=c(lnA1_ode=lnA1, Ea_ode=Ea1)
      tI=ft[nr==i]
      if (method==1) mD_N = "euler" #---------------------------Selection of the solver (see details in deSolve package)
      if (method==2) mD_N = "rk2"
      if (method==4) mD_N = "rk4"
      if ((method==2)||(method==4)) out=as.data.frame(rk(x, tI, v2Func, parms,method = mD_N))
      if (method==5) out=as.data.frame(ode(x, tI, v2Func, parms, method = rkMethod("rk45ck", densetype = NULL, nknots = 0)))
      if (method==6) out=as.data.frame(ode(x, tI, v2Func, parms, method = rkMethod("rk45dp7", densetype = NULL, nknots = 0)))
      if (method==9) out=as.data.frame(lsoda(x, tI, v2Func, parms, atol = 1e-10, rtol = 1e-10))
      if (method==7) out=as.data.frame(ode(x, tI, v2Func, parms, method = rkMethod("rk45ck")))
      if (method==8) out=as.data.frame(ode(x, tI, v2Func, parms, method = rkMethod("rk45dp7")))
      alp=out$a
      if (supY==T) { #-------------------------------------Whether additionally suppress NaN's and outliers
        alp[is.na(alp)]=2
        unity=0
        for (j in 2:length(alp)) {
          if (unity==0) {
            if (alp[j]<0) alp[j]=0 
            else {
              if (alp[j]>amax) 
                if (alp[j-1]>=amax-0.2) { 
                  alp[j]=amax
                  unity=1 }
            } }
          else alp[j]=amax  } }
      if (Dtype!=2) {
        alph=suppressWarnings(splinefun(tI,alp,method = "fmm"))  #----computes da/dt
        dadt=alph(tI,deriv=1) 
        if(i==1) dadtI=dadt
        else dadtI=c(dadtI,dadt) }
      if(i==1) alpI=alp
      else alpI=c(alpI,alp)
    }
    if (Dtype==3) res=dadtI*WghDA
    else  if(Dtype==1) res=dadtI
    else res=alpI
    res
  }
    
  DMM58 = function (lnA1,Ea1,mu,lnA2,Ea2,oT,al) { #-Dubovitskiy-Manelis-Merzhanov(1958) function
    if (mu==1) mpart=0
    else {
      if (mu==0) mpart=al*(1-al)
      else mpart=(1-mu)*al*(1-al)/(1-mu*al) }
    exp(lnA1-Ea1*1000*oT)*(1-al)+exp(lnA2-Ea2*1000*oT)*mpart
  }
  
  #-Dubovitskiy-Manelis-Merzhanov-1958 model  
  funK3s=function(lnA1,Ea1,mu,lnA2,lnA22,Ea2,fT1,ft,fa,flag,nr,Dtype,WghDA,method,amin,amax,A_trigg) {  
    alpLim=amin
    Nal=length(fa)
    alpI=numeric(Nal)
    dadtI=numeric(Nal)
    alpI[1]=alpLim
    dadtI[1]=0
    if (Ea2==0) Ea2=Ea1
    for (i in 2:Nal) {
      if ((nr[i]==nr[i-1]+1)||((nr[i-1]>2)&&(nr[i]==1))) {  #switch between runs in row (conversion=1 to 0)
        alpI[i]=alpLim
        dadtI[i]=0 }
      else {
        h=ft[i]-ft[i-1]
        if ((flag[i]=="d2")&&(lnA22!=0)) lnA2_one=lnA22
        else lnA2_one=lnA2
        if (method==1) {      # crude Euler method, delete after all!
          alpI[i]=alpI[i-1]+h*dadtI[i-1]
          dadtI[i]=DMM58(lnA1+A_trigg[i],Ea1,mu,lnA2_one+A_trigg[i],Ea2,fT1[i],alpI[i]) }
        if (method==2) {      # Improved Euler method, a second order!
          approx_al=alpI[i-1]+h*dadtI[i-1]
          alpI[i]=alpI[i-1]+h*(dadtI[i-1]+DMM58(lnA1+A_trigg[i],Ea1,mu,lnA2_one+A_trigg[i],Ea2,fT1[i],approx_al))/2
          dadtI[i]=DMM58(lnA1+A_trigg[i],Ea1,mu,lnA2_one+A_trigg[i],Ea2,fT1[i],alpI[i]) }
        #if(is.na(alpI[i])) cat(lnA1+A_trigg[i],Ea1,mu,lnA2,Ea2,i,alpI[i-1],dadtI[i-1],alpI[i],fT1[i],dadtI[i],"\n",file="funk2_error.txt", append = T)
        if ((is.na(alpI[i]))&&(alpI[i-1]>=amax-0.1)) {
          alpI[i]=1
          dadtI[i]=0  }
        else {
          if (is.na(alpI[i])) {
            alpI[i]=2
            dadtI[i]=0 }
          if (alpI[i]>=amax) {
            alpI[i]=1
            dadtI[i]=0 } } } 
    }
    if (Dtype==3) res=dadtI*WghDA
    else  if(Dtype==1) res=dadtI
    else res=alpI
    res
  }

  funK3=function(lnA1,Ea1,mu,lnA2,lnA22,Ea2,fT1,ft,fa,flag,nr,Dtype,WghDA,method,amin,amax,Nruns,TempF,supY,ti_trigg) {  
    alpLim=amin
    for (i in 1:Nruns) {
      v3Func <- function(t, x, parms) {
        with(as.list(parms), {
          lnAt1=ifelse(t>ti_trigg[i], lnA1_ode, lnA1_ode-log(Kt))
          lnAt2=ifelse(t>ti_trigg[i], lnA2_ode, lnA2_ode-log(Kt))
          dx = exp(lnAt1-Ea1_ode*1000/8.31446/TempF[[i]](t))*(1-x)+
            exp(lnAt2-Ea2_ode*1000/8.31446/TempF[[i]](t))*(1-mu_ode)*x*(1-x)/(1-mu_ode*x)
          list(dx)
        }) }
      x=c(a = alpLim)
      if (Ea2==0) Ea2=Ea1
      if ((flag[nr==i][1] == "d2")&&(lnA22!=0)) lnA2_one=lnA22
      else lnA2_one=lnA2
      parms=c(lnA1_ode=lnA1, Ea1_ode=Ea1, mu_ode=mu, lnA2_ode=lnA2_one, Ea2_ode=Ea2)
      tI=ft[nr==i]
      if (method==1) mD_N = "euler" #---------------------------Selection of the solver (see details in deSolve package)
      if (method==2) mD_N = "rk2"
      if (method==4) mD_N = "rk4"
      if ((method==2)||(method==4)) out=as.data.frame(rk(x, tI, v3Func, parms,method = mD_N))
      if (method==5) out=as.data.frame(ode(x, tI, v3Func, parms, method = rkMethod("rk45ck", densetype = NULL, nknots = 0)))
      if (method==6) out=as.data.frame(ode(x, tI, v3Func, parms, method = rkMethod("rk45dp7", densetype = NULL, nknots = 0)))
      if (method==9) out=as.data.frame(lsoda(x, tI, v3Func, parms, atol = 1e-10, rtol = 1e-10))
      if (method==7) out=as.data.frame(ode(x, tI, v3Func, parms, method = rkMethod("rk45ck")))
      if (method==8) out=as.data.frame(ode(x, tI, v3Func, parms, method = rkMethod("rk45dp7")))
      alp=out$a
      if (supY==T) { #-------------------------------------Whether additionally suppress NaN's and outliers
        alp[is.na(alp)]=2
        unity=0
        for (j in 2:length(alp)) {
            if (unity==0) {
              if (alp[j]<0) alp[j]=0 
              else {
                if (alp[j]>amax) 
                  if (alp[j-1]>=amax-0.2) { 
                    alp[j]=amax
                    unity=1 }
              } }
             else alp[j]=amax }   }
      if (Dtype!=2) {
        alph=suppressWarnings(splinefun(tI,alp,method = "fmm"))  #----computes da/dt
        dadt=alph(tI,deriv=1) 
        if(i==1) dadtI=dadt
        else dadtI=c(dadtI,dadt) }
      if(i==1) alpI=alp
      else alpI=c(alpI,alp)
    }
    if (Dtype==3) res=dadtI*WghDA
    else  if(Dtype==1) res=dadtI
    else res=alpI
    res
  }
    
  funK4s=function(lnA1,Ea1,n1,m1,lnA2,Ea2,n2,m2,fT1,ft,fa,nr,Dtype,WghDA,method,q1,q2,amin,amax,A_trigg) {  #two independent reactions
    alpLim=amin
    Nal=length(fa)
    alpI=numeric(Nal)
    dadtI=numeric(Nal)
    alpI[1]=alpLim
    dadtI[1]=0
    if (Ea2==0) Ea2=Ea1
    for (i in 2:Nal) {
      if ((nr[i]==nr[i-1]+1)||((nr[i-1]>2)&&(nr[i]==1))) {  #switch between runs in row (conversion=1 to 0)
        alpI[i]=alpLim
        dadtI[i]=0 }
      else {
        h=ft[i]-ft[i-1]
        if (method==1) {      # crude (??) Euler method
          alpI[i]=alpI[i-1]+h*dadtI[i-1]
          dadtI[i]=ePT(lnA1+A_trigg[i],Ea1,n1,m1,q1,fT1[i],alpI[i])+ePT(lnA2+A_trigg[i],Ea2,n2,m2,q2,fT1[i],alpI[i]) 
          }
        if (method==2) {      # Improved Euler method, a second order!
          approx_al=alpI[i-1]+h*dadtI[i-1]
          alpI[i]=alpI[i-1]+h*(dadtI[i-1]+ePT(lnA1+A_trigg[i],Ea1,n1,m1,q1,fT1[i],approx_al)+ePT(lnA2+A_trigg[i],Ea2,n2,m2,q2,fT1[i],approx_al))/2
          dadtI[i]=ePT(lnA1+A_trigg[i],Ea1,n1,m1,q1,fT1[i],alpI[i])+ePT(lnA2+A_trigg[i],Ea2,n2,m2,q2,fT1[i],alpI[i])  }
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
    }
    #    write.table(cbind(fT1,dadtI1,dadtI2,dadtI), file="SB+SB_diffStepsTGA.txt", eol="\n")
    if (Dtype==3) res=dadtI*WghDA
    else  if(Dtype==1) res=dadtI
    else res=alpI
    res
  }
  
funK4=function(lnA1,Ea1,n1,m1,lnA2,Ea2,n2,m2,fT1,ft,fa,nr,Dtype,WghDA,method,q1,q2,amin,amax,Nruns,TempF,supY,ti_trigg) {  
  alpLim=amin
  for (i in 1:Nruns) {
    if (Ea2==0) Ea2=Ea1
    if (m1!=10) {
      v4Func <- function(t, x, parms) {
        with(as.list(parms), {
          x=ifelse(x<MyUnity,x,MyUnity)
          lnAt1=ifelse(t>ti_trigg[i], lnA1_ode, lnA1_ode-log(Kt))
          lnAt2=ifelse(t>ti_trigg[i], lnA2_ode, lnA2_ode-log(Kt))
          dx = exp(lnAt1-Ea1_ode*1000/8.31446/TempF[[i]](t))*((1-x)^n1_ode)*(1-q1_ode*(1-x))^m1_ode +
            exp(lnAt2-Ea2_ode*1000/8.31446/TempF[[i]](t))*((1-x)^n2_ode)*(1-q2_ode*(1-x))^m2_ode
          list(dx)
        }) }
      parms=c(lnA1_ode=lnA1, Ea1_ode=Ea1, n1_ode=n1, m1_ode=m1, q1_ode=q1,
              lnA2_ode=lnA2, Ea2_ode=Ea2, n2_ode=n2, m2_ode=m2, q2_ode=q2)
    }
    else {
      v4Func <- function(t, x, parms) { #----reaction model in extended KJMAE form
        with(as.list(parms), {
          x=ifelse(x<MyUnity,x,MyUnity)
          lnAt1=ifelse(t>ti_trigg[i], lnA1_ode, lnA1_ode-log(Kt))
          lnAt2=ifelse(t>ti_trigg[i], lnA2_ode, lnA2_ode-log(Kt))
          dx = exp(lnAt1-Ea1_ode*1000/8.31446/TempF[[i]](t))*n1_ode*(1-x)*(-log(1-x))^((n1_ode-1)/n1_ode) +
            exp(lnAt2-Ea2_ode*1000/8.31446/TempF[[i]](t))*((1-x)^n2_ode)*(1-q2_ode*(1-x))^m2_ode
          list(dx)
        }) } 
      parms=c(lnA1_ode=lnA1, Ea1_ode=Ea1, n1_ode=n1, 
              lnA2_ode=lnA2, Ea2_ode=Ea2, n2_ode=n2, m2_ode=m2, q2_ode=q2)
    }
    x=c(a = alpLim)
    tI=ft[nr==i]
    if (method==1) mD_N = "euler" #---------------------------Selection of the solver (see details in deSolve package)
    if (method==2) mD_N = "rk2"
    if (method==4) mD_N = "rk4"
    if ((method==2)||(method==4)) out=as.data.frame(rk(x, tI, v4Func, parms,method = mD_N))
    if (method==5) out=as.data.frame(ode(x, tI, v4Func, parms, method = rkMethod("rk45ck", densetype = NULL, nknots = 0)))
    if (method==6) out=as.data.frame(ode(x, tI, v4Func, parms, method = rkMethod("rk45dp7", densetype = NULL, nknots = 0)))
    if (method==9) out=as.data.frame(lsoda(x, tI, v4Func, parms, atol = 1e-10, rtol = 1e-10))
    #these solvers only for initial guesses!
    if (method==7) out=as.data.frame(ode(x, tI, v4Func, parms, method = rkMethod("rk45ck")))
    if (method==8) out=as.data.frame(ode(x, tI, v4Func, parms, method = rkMethod("rk45dp7")))
    alp=out$a
    if (supY==T) { #-------------------------------------Whether additionally suppress NaN's and outliers
      alp[is.na(alp)]=2
      unity=0
      for (j in 2:length(alp)) {
        if (unity==0) {
          if (alp[j]<0) alp[j]=0 
          else {
            if (alp[j]>amax) 
              if (alp[j-1]>=amax-0.2) { 
                alp[j]=amax
                unity=1 }
          } }
        else alp[j]=amax }   }
    if (Dtype!=2) {
      alph=suppressWarnings(splinefun(tI,alp,method = "fmm"))  #----computes da/dt
      dadt=alph(tI,deriv=1) 
      if(i==1) dadtI=dadt
      else dadtI=c(dadtI,dadt) }
    if(i==1) alpI=alp
    else alpI=c(alpI,alp)
  }
  if (Dtype==3) res=dadtI*WghDA
  else  if(Dtype==1) res=dadtI
  else res=alpI
  res
}

  funK5s=function(lnA1,Ea1,n1,m1,lnA2,Ea2,n2,m2,cd1,cd2,fT1,ft,fa,flag,nr,Dtype,WghDA,method,q1,q2,amin,amax,A_trigg) {  #two consecutive reactions
    alpLim=amin
    Nal=length(fa)
    alpI=numeric(Nal)
    dadtI=numeric(Nal)
    alpI1=numeric(Nal)
    dadtI1=numeric(Nal)
    alpI2=numeric(Nal)
    dadtI2=numeric(Nal)
    alpI[1]=alpLim
    dadtI[1]=0
    alpI1[1]=alpLim
    dadtI1[1]=0
    alpI2[1]=alpLim
    dadtI2[1]=0
    if (Ea2==0) Ea2=Ea1
    for (i in 2:Nal) {
      if ((nr[i]==nr[i-1]+1)||((nr[i-1]>2)&&(nr[i]==1))) {  #switch between runs in row (conversion=1 to 0)
        alpI[i]=alpLim
        dadtI[i]=0
        alpI1[i]=alpLim
        dadtI1[i]=0
        alpI2[i]=alpLim
        dadtI2[i]=0 }
      else {
        h=ft[i]-ft[i-1]
        if (method==1) {      # crude Euler method, delete after all(?)
          alpI1[i]=alpI1[i-1]+h*dadtI1[i-1]
          if (alpI1[i]>=amax) {
            alpI1[i]=1
            dadtI1[i]=0 } 
          else dadtI1[i]=ePT(lnA1+A_trigg[i],Ea1,n1,m1,q1,fT1[i],alpI1[i]) 
          alpI2[i]=alpI2[i-1]+h*dadtI2[i-1]
          if (alpI2[i]>=amax) {
            alpI2[i]=1
            dadtI2[i]=0 } 
          else {
            if (alpI2[i]>alpI1[i]) {
              alpI2[i]=alpI1[i] #ifelse(alpI2[i-1]>alpI1[i],alpI2[i-1],alpI1[i]) 
              dadtI2[i]=0 }
            else dadtI2[i]=ePT(lnA2,Ea2,n2,m2,q2,fT1[i],1-alpI1[i]+alpI2[i]) }
        }
        if (is.na(alpI1[i])) alpI1[i]=0
        if (is.na(alpI2[i])) alpI2[i]=0
        alpI[i]=alpI[i-1]+h*dadtI[i-1] 
        if (alpI[i]>=amax) {
          alpI[i]=1
          dadtI[i]=0 } 
        else {
          if (alpI[i]>alpI1[i]) { #----to prevent error when the second stage kinetics is exceptionally high(!???)
            alpI[i]=alpI1[i]
            dadtI[i]=0 }
          else {
            if (flag[i]=="d2") dadtI[i]=cd2*dadtI1[i]+(1-cd2)*dadtI2[i]
            else dadtI[i]=cd1*dadtI1[i]+(1-cd1)*dadtI2[i]
            }
          }
        }
    }
    #    write.table(cbind(fT1,dadtI1,dadtI2,dadtI), file="SB+SB_diffStepsTGA.txt", eol="\n")
    if (Dtype==3) res=dadtI*WghDA
    else  if(Dtype==1) res=dadtI
    else res=alpI
    res
  }

  funK5=function(lnA1,Ea1,n1,m1,lnA2,Ea2,n2,m2,cd1,cd2,fT1,ft,fa,flag,nr,Dtype,WghDA,method,q1,q2,amin,amax,Nruns,TempF,supY,ti_trigg) {  
    alpLim=amin
    for (i in 1:Nruns) {
      if (Ea2==0) Ea2=Ea1
      if (m1!=10) {
        v5Func <- function(t, x, parms) {
          with(as.list(c(parms, x)),{
            lnAt=ifelse(t>ti_trigg[i], lnA1_ode, lnA1_ode-log(Kt))
            A=ifelse(A>MyNull,A,MyNull)
            B=ifelse(B>MyNull,B,MyNull)
            C=ifelse(C>MyNull,C,MyNull)
            dA = -exp(lnAt-Ea1_ode*1000/8.31446/TempF[[i]](t))*(A^n1_ode)*(1-q1_ode*(1-B))^m1_ode #---note that here A = 1 - a !
            dB = exp(lnAt-Ea1_ode*1000/8.31446/TempF[[i]](t))*A^n1_ode*(1-q1_ode*(1-B))^m1_ode - exp(lnA2_ode-Ea2_ode*1000/8.31446/TempF[[i]](t))*B^n2_ode*(1-q2_ode*(1-C))^m2_ode
            dC = 1-A-B #exp(lnA2_ode-Ea2_ode*1000/8.31446/TempF[[i]](t))*B^n2_ode*(1-q2_ode*(1-C))^m2_ode #1-A-B
            return(list(c(dA,dB,dC)))
          }) }
        parms=c(lnA1_ode=lnA1, Ea1_ode=Ea1, n1_ode=n1, m1_ode=m1, q1_ode=q1,
                lnA2_ode=lnA2, Ea2_ode=Ea2, n2_ode=n2, m2_ode=m2, q2_ode=q2)
      }
      else {
        v5Func <- function(t, x, parms) {#----reaction model in extended KJMAE form
          with(as.list(c(parms, x)),{
            lnAt=ifelse(t>ti_trigg[i], lnA1_ode, lnA1_ode-log(Kt))
            A=ifelse(A>MyNull,A,MyNull)
            B=ifelse(B>MyNull,B,MyNull)
            C=ifelse(C>MyNull,C,MyNull)
            dA = -exp(lnAt-Ea1_ode*1000/8.31446/TempF[[i]](t))*n1_ode*A*(-log(A))^((n1_ode-1)/n1_ode) #---note that here A = 1 - a !
            dB = exp(lnAt-Ea1_ode*1000/8.31446/TempF[[i]](t))*n1_ode*A*(-log(A))^((n1_ode-1)/n1_ode) - exp(lnA2_ode-Ea2_ode*1000/8.31446/TempF[[i]](t))*B^n2_ode*(1-q2_ode*(1-C))^m2_ode
            dC = exp(lnA2_ode-Ea2_ode*1000/8.31446/TempF[[i]](t))*B^n2_ode*(1-q2_ode*(1-C))^m2_ode #1-A-B
            return(list(c(dA,dB,dC)))
          }) }
        parms=c(lnA1_ode=lnA1, Ea1_ode=Ea1, n1_ode=n1, 
                lnA2_ode=lnA2, Ea2_ode=Ea2, n2_ode=n2, m2_ode=m2, q2_ode=q2)
      }
      #x=c(A=0.99995, B=0.00003, C=0.00002) #---to compare the outpur with Netzsch
      x=c(A=1-alpLim, B=alpLim, C=alpLim)
      tI=ft[nr==i]
      if (method==1) mD_N = "euler" #---------------------------Selection of the solver (see details in deSolve package)
      if (method==2) mD_N = "rk2"
      if (method==4) mD_N = "rk4"
      if ((method==2)||(method==4)) out=as.data.frame(rk(x, tI, v5Func, parms,method = mD_N))
      if (method==5) out=as.data.frame(ode(x, tI, v5Func, parms, method = rkMethod("rk45ck", densetype = NULL, nknots = 0)))
      if (method==6) out=as.data.frame(ode(x, tI, v5Func, parms, method = rkMethod("rk45dp7", densetype = NULL, nknots = 0)))
      if (method==9) out=as.data.frame(lsoda(x, tI, v5Func, parms, atol = 1e-10, rtol = 1e-10))
      if (method==7) out=as.data.frame(ode(x, tI, v5Func, parms, method = rkMethod("rk45ck")))
      if (method==8) out=as.data.frame(ode(x, tI, v5Func, parms, method = rkMethod("rk45dp7")))
      
      fraction=ifelse(flag[nr==i][1] == "d2", cd2, cd1)
      alp1=1-out$A
      alp2=1-out$A-out$B
      if (supY==T) { #-------------------------------------Whether additionally suppress NaN's and outliers
        alp1[is.na(alp1)]=2
        unity=0
        for (j in 2:length(alp1)) {
          if (unity==0) {
            if (alp1[j]<0) alp1[j]=0 
            else {
              if (alp1[j]>amax) 
                if (alp1[j-1]>=amax-0.2) { 
                  alp1[j]=amax
                  unity=1 }
            } }
          else alp1[j]=amax }
        alp2[is.na(alp2)]=2
        unity=0
        for (j in 2:length(alp2)) {
          if (unity==0) {
            if (alp2[j]<0) alp2[j]=0 
            else {
              if (alp2[j]>amax) 
                if (alp2[j-1]>=amax-0.2) { 
                  alp2[j]=amax
                  unity=1 }
            } }
          else alp2[j]=amax }
      }
      alp=fraction*alp1+(1-fraction)*alp2
      if (Dtype!=2) {
        alph=suppressWarnings(splinefun(tI,alp,method = "fmm"))  #----computes da/dt
        dadt=alph(tI,deriv=1) 
        if(i==1) dadtI=dadt
        else dadtI=c(dadtI,dadt) }
      if(i==1) alpI=alp
      else alpI=c(alpI,alp)
    }
    if (Dtype==3) res=dadtI*WghDA
    else  if(Dtype==1) res=dadtI
    else res=alpI
    res
  }
  
  funK6s=function(lnA1,Ea1,n1,m1,lnA2,Ea2,n2,m2,cd1,cd2,fT1,ft,fa,flag,nr,Dtype,WghDA,method,q1,q2,amin,amax,A_trigg) {  #two independent reactions
    alpLim=amin
    Nal=length(fa)
    alpI=numeric(Nal)
    dadtI=numeric(Nal)
    alpI1=numeric(Nal)
    dadtI1=numeric(Nal)
    alpI2=numeric(Nal)
    dadtI2=numeric(Nal)
    alpI[1]=alpLim
    dadtI[1]=0
    alpI1[1]=alpLim
    dadtI1[1]=0
    alpI2[1]=alpLim
    dadtI2[1]=0
    if (Ea2==0) Ea2=Ea1
    for (i in 2:Nal) {
      if ((nr[i]==nr[i-1]+1)||((nr[i-1]>2)&&(nr[i]==1))) {  #switch between runs in row (conversion=1 to 0)
        alpI[i]=alpLim
        dadtI[i]=0
        alpI1[i]=alpLim
        dadtI1[i]=0
        alpI2[i]=alpLim
        dadtI2[i]=0 }
      else {
        h=ft[i]-ft[i-1]
        if (method==1) {      # crude Euler method, delete after all!
          alpI1[i]=alpI1[i-1]+h*dadtI1[i-1]
          dadtI1[i]=ePT(lnA1+A_trigg[i],Ea1,n1,m1,q1,fT1[i],alpI1[i]) 
          alpI2[i]=alpI2[i-1]+h*dadtI2[i-1]
          dadtI2[i]=ePT(lnA2+A_trigg[i],Ea2,n2,m2,q2,fT1[i],alpI2[i]) 
          }
        if (method==2) {      # Improved Euler method, a second order!
          approx_al1=alpI1[i-1]+h*dadtI1[i-1]
          alpI1[i]=alpI1[i-1]+h*(dadtI1[i-1]+ePT(lnA1+A_trigg[i],Ea1,n1,m1,q1,fT1[i],approx_al1))/2
          dadtI1[i]=ePT(lnA1+A_trigg[i],Ea1,n1,m1,q1,fT1[i],alpI1[i]) 
          approx_al2=alpI2[i-1]+h*dadtI2[i-1]
          alpI2[i]=alpI2[i-1]+h*(dadtI2[i-1]+ePT(lnA2+A_trigg[i],Ea2,n2,m2,q2,fT1[i],approx_al2))/2
          dadtI2[i]=ePT(lnA2+A_trigg[i],Ea2,n2,m2,q2,fT1[i],alpI2[i])  }
        ###
        if (is.na(alpI1[i])) alpI1[i]=0
        if (is.na(alpI2[i])) alpI2[i]=0
        if (alpI1[i]>=amax) {
          alpI1[i]=1
          dadtI1[i]=0 } 
        if (alpI2[i]>=amax) {
          alpI2[i]=1
          dadtI2[i]=0 } 
        alpI[i]=alpI[i-1]+h*dadtI[i-1] 
        #if(is.na(alpI[i])) cat(lnA1,Ea1,n1,m1,q,q2,i,alpI[i-1],dadtI[i-1],alpI[i],fT1[i],dadtI[i],"\n",file="aaaa.txt", append = T)
        if (flag[i]=="d2") dadtI[i]=cd2*dadtI1[i]+(1-cd2)*dadtI2[i]
        else dadtI[i]=cd1*dadtI1[i]+(1-cd1)*dadtI2[i]
        }
    }
#    write.table(cbind(fT1,dadtI1,dadtI2,dadtI), file="SB+SB_diffStepsTGA.txt", eol="\n")
    if (Dtype==3) res=dadtI*WghDA
    else  if(Dtype==1) res=dadtI
    else res=alpI
    res
  }

  funK6=function(lnA1,Ea1,n1,m1,lnA2,Ea2,n2,m2,cd1,cd2,fT1,ft,fa,flag,nr,Dtype,WghDA,method,q1,q2,amin,amax,Nruns,TempF,supY,ti_trigg) {  
    alpLim=amin
    for (i in 1:Nruns) {
      if (Ea2==0) Ea2=Ea1
      if (m1!=10) {
        v6Func <- function(t, x, parms) {
          with(as.list(c(parms, x)),{
            lnAt1=ifelse(t>ti_trigg[i], lnA1_ode, lnA1_ode-log(Kt))
            #lnAt2=ifelse(t>ti_trigg[i], lnA2_ode, lnA2_ode-log(Kt))
            A=ifelse(A>MyNull,A,MyNull)
            C=ifelse(C>MyNull,C,MyNull)
            dA = -exp(lnAt1-Ea1_ode*1000/8.31446/TempF[[i]](t))*(A^n1_ode)*(1-q1_ode*A)^m1_ode #---note that here A = 1 - a !
            dC = -exp(lnA2_ode-Ea2_ode*1000/8.31446/TempF[[i]](t))*C^n2_ode*(1-q2_ode*C)^m2_ode
            return(list(c(dA,dC)))
          }) }
        parms=c(lnA1_ode=lnA1, Ea1_ode=Ea1, n1_ode=n1, m1_ode=m1, q1_ode=q1,
                lnA2_ode=lnA2, Ea2_ode=Ea2, n2_ode=n2, m2_ode=m2, q2_ode=q2)
      }
      else {
        v6Func <- function(t, x, parms) {#----reaction model in extended KJMAE form
          with(as.list(c(parms, x)),{
            lnAt1=ifelse(t>ti_trigg[i], lnA1_ode, lnA1_ode-log(Kt))
            #lnAt2=ifelse(t>ti_trigg[i], lnA2_ode, lnA2_ode-log(Kt))
            A=ifelse(A>MyNull,A,MyNull)
            C=ifelse(C>MyNull,C,MyNull)
            dA = -exp(lnAt1-Ea1_ode*1000/8.31446/TempF[[i]](t))*n1_ode*A*(-log(A))^((n1_ode-1)/n1_ode) #---note that here A = 1 - a !
            dC = -exp(lnA2_ode-Ea2_ode*1000/8.31446/TempF[[i]](t))*C^n2_ode*(1-q2_ode*C)^m2_ode
            return(list(c(dA,dC)))
          }) }
        parms=c(lnA1_ode=lnA1, Ea1_ode=Ea1, n1_ode=n1, 
                lnA2_ode=lnA2, Ea2_ode=Ea2, n2_ode=n2, m2_ode=m2, q2_ode=q2)
      }
      x=c(A=1-alpLim, C=1-alpLim)
      tI=ft[nr==i]
      if (method==1) mD_N = "euler" #---------------------------Selection of the solver (see details in deSolve package)
      if (method==2) mD_N = "rk2"
      if (method==4) mD_N = "rk4"
      if ((method==2)||(method==4)) out=as.data.frame(rk(x, tI, v6Func, parms,method = mD_N))
      if (method==5) out=as.data.frame(ode(x, tI, v6Func, parms, method = rkMethod("rk45ck", densetype = NULL, nknots = 0)))
      if (method==6) out=as.data.frame(ode(x, tI, v6Func, parms, method = rkMethod("rk45dp7", densetype = NULL, nknots = 0)))
      if (method==9) out=as.data.frame(lsoda(x, tI, v6Func, parms, atol = 1e-10, rtol = 1e-10))
      if (method==7) out=as.data.frame(ode(x, tI, v6Func, parms, method = rkMethod("rk45ck")))
      if (method==8) out=as.data.frame(ode(x, tI, v6Func, parms, method = rkMethod("rk45dp7")))
      
      fraction=ifelse(flag[nr==i][1] == "d2", cd2, cd1)
      alp1=1-out$A
      alp2=1-out$C
      if (supY==T) { #-------------------------------------Whether additionally suppress NaN's and outliers
        alp1[is.na(alp1)]=2
        unity=0
        for (j in 2:length(alp1)) {
          if (unity==0) {
            if (alp1[j]<0) alp1[j]=0 
            else {
              if (alp1[j]>amax) 
                if (alp1[j-1]>=amax-0.2) { 
                  alp1[j]=amax
                  unity=1 }
            } }
          else alp1[j]=amax }
        alp2[is.na(alp2)]=2
        unity=0
        for (j in 2:length(alp2)) {
          if (unity==0) {
            if (alp2[j]<0) alp2[j]=0 
            else {
              if (alp2[j]>amax) 
                if (alp2[j-1]>=amax-0.2) { 
                  alp2[j]=amax
                  unity=1 }
            } }
          else alp2[j]=amax }
        }
      alp=fraction*alp1+(1-fraction)*alp2
      if (Dtype!=2) {
        alph=suppressWarnings(splinefun(tI,alp,method = "fmm"))  #----computes da/dt
        dadt=alph(tI,deriv=1) 
        if(i==1) dadtI=dadt
        else dadtI=c(dadtI,dadt) }
      if(i==1) alpI=alp
      else alpI=c(alpI,alp)
    }
    if (Dtype==3) res=dadtI*WghDA
    else  if(Dtype==1) res=dadtI
    else res=alpI
    res
  }
  
  funK7s=function(lnA1,Ea1,n1,m1,lnA2,Ea2,n2,m2,cw1,fT1,ft,fa,flag,nr,Dtype,WghDA,method,q1,q2,amin,amax) {  
    #first reaction switches to second
    alpLim=amin
    cw2=cw1
    Nal=length(fa)
    alpI=numeric(Nal)
    dadtI=numeric(Nal)
    alpI[1]=alpLim
    dadtI[1]=0
    if (Ea2==0) Ea2=Ea1
    for (i in 2:Nal) {
      if ((nr[i]==nr[i-1]+1)||((nr[i-1]>2)&&(nr[i]==1))) {  #switch between runs in row (conversion=1 to 0)
        alpI[i]=alpLim
        dadtI[i]=0 }
      else {
        h=ft[i]-ft[i-1]
        cw=ifelse(flag[i]=="d2",cw2,cw1)
        if (method==1) {      # crude Euler method, delete after all!
          alpI[i]=alpI[i-1]+h*dadtI[i-1]
          dadtI[i]=ifelse(alpI[i]<=cw,ePT(lnA1,Ea1,n1,m1,q1,fT1[i],alpI[i]),ePT(lnA2,Ea2,n2,m2,q2,fT1[i],alpI[i])) 
        }
        if (method==2) {      # Improved Euler method, a second order!
          approx_al=alpI[i-1]+h*dadtI[i-1]
          appr_dadt=ifelse(approx_al<=cw,ePT(lnA1,Ea1,n1,m1,q1,fT1[i],approx_al),ePT(lnA2,Ea2,n2,m2,q2,fT1[i],approx_al)) 
          alpI[i]=alpI[i-1]+h*(dadtI[i-1]+appr_dadt)/2
          dadtI[i]=ifelse(alpI[i]<=cw,ePT(lnA1,Ea1,n1,m1,q1,fT1[i],alpI[i]),ePT(lnA2,Ea2,n2,m2,q2,fT1[i],alpI[i]))  }
        ###
        
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
        #if(is.na(alpI[i])) cat(lnA1,Ea1,n1,m1,q1,q2,i,alpI[i-1],dadtI[i-1],alpI[i],fT1[i],dadtI[i],"\n",file="aaaa.txt", append = T)
      }
    }
    #    write.table(cbind(fT1,dadtI1,dadtI2,dadtI), file="SB+SB_diffStepsTGA.txt", eol="\n")
    if (Dtype==3) res=dadtI*WghDA
    else  if(Dtype==1) res=dadtI
    else res=alpI
    res
  }

  funK7=function(lnA1,Ea1,n1,m1,lnA2,Ea2,n2,m2,swA,fT1,ft,fa,nr,Dtype,WghDA,method,q1,q2,amin,amax,Nruns,TempF,supY) {  
    # 1(2) stages: ePT switched at certain swA conversion to another ePT 
    alpLim=amin
    for (i in 1:Nruns) {
      if (Ea2==0) Ea2=Ea1
      v7Func <- function(t, x, parms) {
          with(as.list(parms), {
            x=ifelse(x<MyUnity,x,MyUnity)
            lnA1s=ifelse(x<=swA,lnA1_ode,lnA2_ode)
            Ea1s=ifelse(x<=swA,Ea1_ode,Ea2_ode)
            n1s=ifelse(x<=swA,n1_ode,n2_ode)
            m1s=ifelse(x<=swA,m1_ode,m2_ode)
            q1s=ifelse(x<=swA,q1_ode,q2_ode)
            dx = exp(lnA1s-Ea1s*1000/8.31446/TempF[[i]](t))*((1-x)^n1s)*(1-q1s*(1-x))^m1s
            list(dx)
          }) }
        parms=c(lnA1_ode=lnA1, Ea1_ode=Ea1, n1_ode=n1, m1_ode=m1, q1_ode=q1,
                lnA2_ode=lnA2, Ea2_ode=Ea2, n2_ode=n2, m2_ode=m2, q2_ode=q2)
      x=c(a = alpLim)
      tI=ft[nr==i]
      if (method==1) mD_N = "euler" #---------------------------Selection of the solver (see details in deSolve package)
      if (method==2) mD_N = "rk2"
      if (method==4) mD_N = "rk4"
      if ((method==2)||(method==4)) out=as.data.frame(rk(x, tI, v7Func, parms,method = mD_N))
      if (method==5) out=as.data.frame(ode(x, tI, v7Func, parms, method = rkMethod("rk45ck", densetype = NULL, nknots = 0)))
      if (method==6) out=as.data.frame(ode(x, tI, v7Func, parms, method = rkMethod("rk45dp7", densetype = NULL, nknots = 0)))
      if (method==9) out=as.data.frame(lsoda(x, tI, v7Func, parms, atol = 1e-10, rtol = 1e-10))
      #these solvers only for initial guesses!
      if (method==7) out=as.data.frame(ode(x, tI, v7Func, parms, method = rkMethod("rk45ck")))
      if (method==8) out=as.data.frame(ode(x, tI, v7Func, parms, method = rkMethod("rk45dp7")))
      alp=out$a
      if (supY==T) { #-------------------------------------Whether additionally suppress NaN's and outliers
        alp[is.na(alp)]=2
        unity=0
        for (j in 2:length(alp)) {
          if (unity==0) {
            if (alp[j]<0) alp[j]=0 
            else {
              if (alp[j]>amax) 
                if (alp[j-1]>=amax-0.2) { 
                  alp[j]=amax
                  unity=1 }
            } }
          else alp[j]=amax }   }
      if (Dtype!=2) {
        alph=suppressWarnings(splinefun(tI,alp,method = "fmm"))  #----computes da/dt
        dadt=alph(tI,deriv=1) 
        if(i==1) dadtI=dadt
        else dadtI=c(dadtI,dadt) }
      if(i==1) alpI=alp
      else alpI=c(alpI,alp)
    }
    if (Dtype==3) res=dadtI*WghDA
    else  if(Dtype==1) res=dadtI
    else res=alpI
    res
  }
  
  funK8s=function(lnA1,Ea1,mu1,lnA2,Ea2,lnA3,Ea3,n3,m3,cd1,cd2,fT1,ft,fa,flag,nr,Dtype,WghDA,method,q3,amin,amax,A_trigg) {
    alpLim=amin
    Nal=length(fa)
    alpI=numeric(Nal)
    dadtI=numeric(Nal)
    alpI1=numeric(Nal)
    dadtI1=numeric(Nal)
    alpI2=numeric(Nal)
    dadtI2=numeric(Nal)
    alpI[1]=alpLim
    dadtI[1]=0
    alpI1[1]=alpLim
    dadtI1[1]=0
    alpI2[1]=alpLim
    dadtI2[1]=0
    if (Ea2==0) Ea2=Ea1
    if (Ea3==0) Ea3=Ea1
    for (i in 2:Nal) {
      if ((nr[i]==nr[i-1]+1)||((nr[i-1]>2)&&(nr[i]==1))) {  #switch between runs in row (conversion=1 to 0)
        alpI[i]=alpLim
        dadtI[i]=0
        alpI1[i]=alpLim
        dadtI1[i]=0
        alpI2[i]=alpLim
        dadtI2[i]=0 }
      else {
        h=ft[i]-ft[i-1]
        if (method==1) {      # crude Euler method, delete after all(?)
          alpI1[i]=alpI1[i-1]+h*dadtI1[i-1]
          if (alpI1[i]>=amax) {
            alpI1[i]=1
            dadtI1[i]=0 } 
          else dadtI1[i]=DMM58(lnA1+A_trigg[i],Ea1,mu1,lnA2+A_trigg[i],Ea2,fT1[i],alpI1[i])
          alpI2[i]=alpI2[i-1]+h*dadtI2[i-1]
          if (alpI2[i]>=amax) {
            alpI2[i]=1
            dadtI2[i]=0 } 
          else {
            if (alpI2[i]>alpI1[i]) {
              alpI2[i]=alpI1[i] #ifelse(alpI2[i-1]>alpI1[i],alpI2[i-1],alpI1[i]) 
              dadtI2[i]=0 }
            else dadtI2[i]=ePT(lnA3,Ea3,n3,m3,q3,fT1[i],1-alpI1[i]+alpI2[i]) }
        }
        if (is.na(alpI1[i])) alpI1[i]=0
        if (is.na(alpI2[i])) alpI2[i]=0
        alpI[i]=alpI[i-1]+h*dadtI[i-1] 
        if (alpI[i]>=amax) {
          alpI[i]=1
          dadtI[i]=0 } 
        else {
          if (alpI[i]>alpI1[i]) { #----to prevent error when the second stage kinetics is exceptionally high(!???)
            alpI[i]=alpI1[i]
            dadtI[i]=0 }
          else {
            if (flag[i]=="d2") dadtI[i]=cd2*dadtI1[i]+(1-cd2)*dadtI2[i]
            else dadtI[i]=cd1*dadtI1[i]+(1-cd1)*dadtI2[i]
          }
        }
      }
    }
    #    write.table(cbind(fT1,dadtI1,dadtI2,dadtI), file="SB+SB_diffStepsTGA.txt", eol="\n")
    if (Dtype==3) res=dadtI*WghDA
    else  if(Dtype==1) res=dadtI
    else res=alpI
    res
  }
  
  funK8=function(lnA1,Ea1,mu1,lnA2,Ea2,lnA3,Ea3,n3,m3,cd1,cd2,fT1,ft,fa,flag,nr,Dtype,WghDA,method,q3,amin,amax,Nruns,TempF,supY,ti_trigg) {  
# 3(2) stages: consecutive general autocatalytic (DMM) followed by ePT reaction
    alpLim=amin
    for (i in 1:Nruns) {
      v8Func <- function(t, x, parms) {
        with(as.list(c(parms, x)),{
          lnAt1=ifelse(t>ti_trigg[i], lnA1_ode, lnA1_ode-log(Kt))
          lnAt2=ifelse(t>ti_trigg[i], lnA2_ode, lnA2_ode-log(Kt))
          A=ifelse(A>MyNull,A,MyNull)
          B=ifelse(B<amax,B,amax)
          C=ifelse(C<amax,C,amax)
          dA = -exp(lnAt1-Ea1_ode*1000/8.31446/TempF[[i]](t))*A -
            exp(lnAt2-Ea2_ode*1000/8.31446/TempF[[i]](t))*(1-mu_ode)*B*A/(1-mu_ode*B) #---note that here A = 1 - a !
          dB = exp(lnAt1-Ea1_ode*1000/8.31446/TempF[[i]](t))*A +
                exp(lnAt2-Ea2_ode*1000/8.31446/TempF[[i]](t))*(1-mu_ode)*B*A/(1-mu_ode*B) -
                exp(lnA3_ode-Ea3_ode*1000/8.31446/TempF[[i]](t))*B^n3_ode*(1-q3_ode*(1-C))^m3_ode
          dC = exp(lnA3_ode-Ea3_ode*1000/8.31446/TempF[[i]](t))*B^n3_ode*(1-q3_ode*(1-C))^m3_ode #1-A-B
          #C=1-A-B
          return(list(c(dA,dB,dC)))
        }) }
      x=c(A=1-alpLim, B=alpLim, C=alpLim)
      if (Ea3==0) Ea3=Ea1 #what for??
      if (Ea2==0) Ea2=Ea1 #what for??
      parms=c(lnA1_ode=lnA1, Ea1_ode=Ea1, mu_ode=mu1, lnA2_ode=lnA2, Ea2_ode=Ea2, 
              lnA3_ode=lnA3, Ea3_ode=Ea3, n3_ode=n3, m3_ode=m3, q3_ode=q3)
      tI=ft[nr==i]
      if (method==1) mD_N = "euler" #---------------------------Selection of the solver (see details in deSolve package)
      if (method==2) mD_N = "rk2"
      if (method==4) mD_N = "rk4"
      if ((method==2)||(method==4)) out=as.data.frame(rk(x, tI, v8Func, parms, method = mD_N))
      #if (method==5) out=as.data.frame(ode(x, tI, v8Func, parms, method = rkMethod("rk45ck", densetype = NULL, nknots = 0)))
      #if (method==6) out=as.data.frame(ode(x, tI, v8Func, parms, method = rkMethod("rk45dp7", densetype = NULL, nknots = 0)))
      if (method==9) out=as.data.frame(lsoda(x, tI, v8Func, parms, atol = 1e-10, rtol = 1e-10))
      #if (method==7) out=as.data.frame(ode(x, tI, v8Func, parms, method = rkMethod("rk45ck")))
      if (method==8) out=as.data.frame(ode(x, tI, v8Func, parms, method = rkMethod("rk45dp7")))
      
      fraction=ifelse(flag[nr==i][1] == "d2", cd2, cd1)
      alp1=1-out$A
      alp2=out$C
      if (supY==T) { #-------------------------------------Whether additionally suppress NaN's and outliers
        alp1[is.na(alp1)]=2
        unity=0
        for (j in 2:length(alp1)) {
          if (unity==0) {
            if (alp1[j]<0) alp1[j]=0 
            else {
              if (alp1[j]>amax) 
                if (alp1[j-1]>=amax-0.2) { 
                  alp1[j]=amax
                  unity=1 }
            } }
          else alp1[j]=amax }
        alp2[is.na(alp2)]=2
        unity=0
        for (j in 2:length(alp2)) {
          if (unity==0) {
            if (alp2[j]<0) alp2[j]=0 
            else {
              if (alp2[j]>amax) 
                if (alp2[j-1]>=amax-0.2) { 
                  alp2[j]=amax
                  unity=1 }
            } }
          else alp2[j]=amax }
      }
      alp=fraction*alp1+(1-fraction)*alp2
      if (Dtype!=2) {
        alph=suppressWarnings(splinefun(tI,alp,method = "fmm"))  #----computes da/dt
        dadt=alph(tI,deriv=1) 
        if(i==1) dadtI=dadt
        else dadtI=c(dadtI,dadt) }
      if(i==1) alpI=alp
      else alpI=c(alpI,alp)
    }
    if (Dtype==3) res=dadtI*WghDA
    else  if(Dtype==1) res=dadtI
    else res=alpI
    res
  }

  funK9=function(lnA1,Ea1,n1,m1,lnA2,Ea2,n2,m2,lnA3,Ea3,n3,m3,cd1,cd2,cr21,cr22,fT1,ft,fa,flag,nr,Dtype,WghDA,method,q1,q2,q3,amin,amax,Nruns,TempF,supY,ti_trigg) {  
    # 3 stages: consecutive ePT followed by two ePT reactions
    alpLim=amin
    alpLim_max=amax
    for (i in 1:Nruns) {
      if (Ea3==0) Ea3=Ea1 #what for??
      if (Ea2==0) Ea2=Ea1 #what for??
      if (m1!=10) {
        v9Func <- function(t, x, parms) {
          with(as.list(c(parms, x)),{
            lnAt=ifelse(t>ti_trigg[i], lnA1_ode, lnA1_ode-log(Kt))
            A=ifelse(A>MyNull,A,MyNull)
            B=ifelse(B>MyNull,B,MyNull)
            C=ifelse(C>MyNull,C,MyNull)
            D=ifelse(D>MyNull,D,MyNull)
            dA = -exp(lnAt-Ea1_ode*1000/8.31446/TempF[[i]](t))*(A^n1_ode)*(1-q1_ode*(1-B))^m1_ode #---note that here A = 1 - a !
            dB = exp(lnAt-Ea1_ode*1000/8.31446/TempF[[i]](t))*A^n1_ode*(1-q1_ode*(1-B))^m1_ode - 
              exp(lnA2_ode-Ea2_ode*1000/8.31446/TempF[[i]](t))*B^n2_ode*(1-q2_ode*(1-C))^m2_ode -
              exp(lnA3_ode-Ea3_ode*1000/8.31446/TempF[[i]](t))*B^n3_ode*(1-q3_ode*(1-D))^m3_ode 
            dC = exp(lnA2_ode-Ea2_ode*1000/8.31446/TempF[[i]](t))*B^n2_ode*(1-q2_ode*(1-C))^m2_ode 
            dD = exp(lnA3_ode-Ea3_ode*1000/8.31446/TempF[[i]](t))*B^n3_ode*(1-q3_ode*(1-D))^m3_ode 
            #D=1-A-B-C
            return(list(c(dA,dB,dC,dD)))
          }) }
        parms=c(lnA1_ode=lnA1, Ea1_ode=Ea1, n1_ode=n1, m1_ode=m1, q1_ode=q1,
                lnA2_ode=lnA2, Ea2_ode=Ea2, n2_ode=n2, m2_ode=m2, q2_ode=q2,
                lnA3_ode=lnA3, Ea3_ode=Ea3, n3_ode=n3, m3_ode=m3, q3_ode=q3)
      }
      else {
        v9Func <- function(t, x, parms) {
          with(as.list(c(parms, x)),{
            lnAt=ifelse(t>ti_trigg[i], lnA1_ode, lnA1_ode-log(Kt))
            A_corr=ifelse(A>=MyNull_for_Log,A,MyNull_for_Log)
            B_corr=ifelse(B>=MyNull,B,MyNull)
            dA = -exp(lnAt-Ea1_ode*1000/8.31446/TempF[[i]](t))*n1_ode*A_corr*(-log(A_corr))^((n1_ode-1)/n1_ode) #---note that here A = 1 - a !
            dB = exp(lnAt-Ea1_ode*1000/8.31446/TempF[[i]](t))*n1_ode*A_corr*(-log(A_corr))^((n1_ode-1)/n1_ode) - 
              exp(lnA2_ode-Ea2_ode*1000/8.31446/TempF[[i]](t))*B_corr^n2_ode*(1-q2_ode*(1-C))^m2_ode -
              exp(lnA3_ode-Ea3_ode*1000/8.31446/TempF[[i]](t))*B_corr^n3_ode*(1-q3_ode*(1-D))^m3_ode 
            dC = exp(lnA2_ode-Ea2_ode*1000/8.31446/TempF[[i]](t))*B_corr^n2_ode*(1-q2_ode*(1-C))^m2_ode 
            dD = exp(lnA3_ode-Ea3_ode*1000/8.31446/TempF[[i]](t))*B_corr^n3_ode*(1-q3_ode*(1-D))^m3_ode 
            #D=1-A-B-C
            return(list(c(dA,dB,dC,dD),A_corr=A_corr,B_corr=B_corr))
          }) }
        parms=c(lnA1_ode=lnA1, Ea1_ode=Ea1, n1_ode=n1, 
                lnA2_ode=lnA2, Ea2_ode=Ea2, n2_ode=n2, m2_ode=m2, q2_ode=q2,
                lnA3_ode=lnA3, Ea3_ode=Ea3, n3_ode=n3, m3_ode=m3, q3_ode=q3)
      }      
      x=c(A=alpLim_max, B=alpLim, C=alpLim, D=alpLim)
      tI=ft[nr==i]
      if (method==1) mD_N = "euler" #---------------------------Selection of the solver (see details in deSolve package)
      if (method==2) mD_N = "rk2"
      if (method==4) mD_N = "rk4"
      if (method<5) out=as.data.frame(rk(x, tI, v9Func, parms, method = mD_N))
      #if (method==5) out=as.data.frame(ode(x, tI, v8Func, parms, method = rkMethod("rk45ck", densetype = NULL, nknots = 0)))
      #if (method==6) out=as.data.frame(ode(x, tI, v8Func, parms, method = rkMethod("rk45dp7", densetype = NULL, nknots = 0)))
      if (method==9) out=as.data.frame(lsoda(x, tI, v9Func, parms, atol = 1e-10, rtol = 1e-10))
      #if (method==7) out=as.data.frame(ode(x, tI, v8Func, parms, method = rkMethod("rk45ck")))
      if (method==8) out=as.data.frame(ode(x, tI, v9Func, parms, method = rkMethod("rk45dp7")))
      
      fraction1=ifelse(flag[nr==i][1] == "d2", cd2, cd1)
      outA=ifelse(m1!=10,out$A,out$A_corr)
      alp=fraction1*(1-outA)+(1-fraction1)*(cr21*out$C+cr22*out$D)
      if (supY==T) { #-------------------------------------Whether additionally suppress NaN's and outliers
        alp[is.na(alp)]=2
        unity=0
        for (j in 2:length(alp)) {
          if (unity==0) {
            if (alp[j]<0) alp[j]=0 
            else {
              if (alp[j]>1) 
                if (alp[j-1]>=amax-0.1) { 
                  alp[j]=1
                  unity=1 }
            } }
          else alp[j]=1 }
      }
      if (Dtype!=2) {
        alph=suppressWarnings(splinefun(tI,alp,method = "fmm"))  #----computes da/dt
        dadt=alph(tI,deriv=1) 
        #dadt=fraction1*out$dA+(1-fraction1)*(cr21*out$dC+cr22*out$dD)
        if(i==1) dadtI=dadt
        else dadtI=c(dadtI,dadt) }
      if(i==1) alpI=alp
      else alpI=c(alpI,alp)
    }
    if (Dtype==3) res=dadtI*WghDA
    else  if(Dtype==1) res=dadtI
    else res=alpI
    res
  }
  
  funK10=function(lnA1,Ea1,n1,m1,lnA2,Ea2,n2,m2,lnA3,Ea3,n3,m3,cd1,cd2,swA,fT1,ft,fa,flag,nr,Dtype,WghDA,method,q1,q2,q3,amin,amax,Nruns,TempF,supY) {  
    # 2(3) stages: consecutive ePT switched at certain swA conversion to another ePT followed by ePT reaction
    swA=1-swA
    alpLim=amin
    alpLim_max=amax
    for (i in 1:Nruns) {
      if (Ea3==0) Ea3=Ea1 #what for??
      if (Ea2==0) Ea2=Ea1 
      #if (lnA2<=1) lnA2=lnA1*lnA2 #an experimental feature, no need in it
    v10Func <- function(t, x, parms) {
      with(as.list(c(parms, x)),{
        A=ifelse(A>MyNull,A,MyNull)
        B=ifelse(B>MyNull,B,MyNull)
        C=ifelse(C>MyNull,C,MyNull)
        lnA1s=ifelse(A>=swA,lnA1_ode,lnA2_ode)
        Ea1s=ifelse(A>=swA,Ea1_ode,Ea2_ode)
        n1s=ifelse(A>=swA,n1_ode,n2_ode)
        m1s=ifelse(A>=swA,m1_ode,m2_ode)
        q1s=ifelse(A>=swA,q1_ode,q2_ode)
        dA = -exp(lnA1s-Ea1s*1000/8.31446/TempF[[i]](t))*(A^n1s)*(1-q1s*(1-B))^m1s #---note that here A = 1 - a !
        dB = exp(lnA1s-Ea1s*1000/8.31446/TempF[[i]](t))*(A^n1s)*(1-q1s*(1-B))^m1s - exp(lnA3_ode-Ea3_ode*1000/8.31446/TempF[[i]](t))*B^n3_ode*(1-q3_ode*(1-C))^m3_ode
        dC = exp(lnA3_ode-Ea3_ode*1000/8.31446/TempF[[i]](t))*B^n3_ode*(1-q3_ode*(1-C))^m3_ode #1-A-B
        return(list(c(dA,dB,dC)))
      }) }
      parms=c(lnA1_ode=lnA1, Ea1_ode=Ea1, n1_ode=n1, m1_ode=m1, q1_ode=q1, swA_ode=swA,
              lnA2_ode=lnA2, Ea2_ode=Ea2, n2_ode=n2, m2_ode=m2, q2_ode=q2,
              lnA3_ode=lnA3, Ea3_ode=Ea3, n3_ode=n3, m3_ode=m3, q3_ode=q3)
      
      x=c(A=alpLim_max, B=alpLim, C=alpLim)
      tI=ft[nr==i]
      if (method==1) mD_N = "euler" #---------------------------Selection of the solver (see details in deSolve package)
      if (method==2) mD_N = "rk2"
      if (method==4) mD_N = "rk4"
      if (method<5) out=as.data.frame(rk(x, tI, v10Func, parms, method = mD_N))
      #if (method==5) out=as.data.frame(ode(x, tI, v8Func, parms, method = rkMethod("rk45ck", densetype = NULL, nknots = 0)))
      #if (method==6) out=as.data.frame(ode(x, tI, v8Func, parms, method = rkMethod("rk45dp7", densetype = NULL, nknots = 0)))
      if (method==9) out=as.data.frame(lsoda(x, tI, v10Func, parms, atol = 1e-10, rtol = 1e-10))
      #if (method==7) out=as.data.frame(ode(x, tI, v8Func, parms, method = rkMethod("rk45ck")))
      if (method==8) out=as.data.frame(ode(x, tI, v10Func, parms, method = rkMethod("rk45dp7")))
      
      fraction=ifelse(flag[nr==i][1] == "d2", cd2, cd1)
      alp1=1-out$A
      alp2=1-out$A-out$B
      if (supY==T) { #-------------------------------------Whether additionally suppress NaN's and outliers
        alp1[is.na(alp1)]=2
        unity=0
        for (j in 2:length(alp1)) {
          if (unity==0) {
            if (alp1[j]<0) alp1[j]=0 
            else {
              if (alp1[j]>amax) 
                if (alp1[j-1]>=amax-0.2) { 
                  alp1[j]=amax
                  unity=1 }
            } }
          else alp1[j]=amax }
        alp2[is.na(alp2)]=2
        unity=0
        for (j in 2:length(alp2)) {
          if (unity==0) {
            if (alp2[j]<0) alp2[j]=0 
            else {
              if (alp2[j]>amax) 
                if (alp2[j-1]>=amax-0.2) { 
                  alp2[j]=amax
                  unity=1 }
            } }
          else alp2[j]=amax }
      }
      alp=fraction*alp1+(1-fraction)*alp2
      if (Dtype!=2) {
        alph=suppressWarnings(splinefun(tI,alp,method = "fmm"))  #----computes da/dt
        dadt=alph(tI,deriv=1) 
        if(i==1) dadtI=dadt
        else dadtI=c(dadtI,dadt) }
      if(i==1) alpI=alp
      else alpI=c(alpI,alp)
    }
    if (Dtype==3) res=dadtI*WghDA
    else  if(Dtype==1) res=dadtI
    else res=alpI
    res
  }
  
  treatData <- function(file) { 
    data=read.table(file,header = TRUE, check.names = FALSE)
    Nruns=ncol(data)/4
    Npts=nrow(data)
    TempF=list()
    ti=data[[2]]
    Te=data[[1]]+273.15
    al=suppressWarnings(as.numeric(as.character(data[[3]])))
    #al[al<1E-6]=0
    da=data[[4]]
    wghDA=da*0+1/(max(da)-min(da)) #to make it a vector
    nr=1+0*data[[4]]
    TempF[[1]]=suppressWarnings(splinefun(ti,Te,method = "fmm")) 
    ti_trigg=numeric(Nruns)
    if (input$TrigTf==T) {
      ind=which(Te>input$TrigTfixed)
      ti_trigg[1]=ti[ind[1]]
      A_trigg=0*da-log(Kt)
      A_trigg[ind]=0 }
    else {
      ti_trigg[1]=ti[1]
      A_trigg=0*da }
    if (Nruns>1) {
      for (i in 2:Nruns) {        
        tim=data[[4*(i-1)+2]]
        ti=c(ti,tim)
        Tem=data[[4*(i-1)+1]]+273.15
        Te=c(Te,Tem)
        al=c(al,data[[4*(i-1)+3]])
        dam=data[[4*(i-1)+4]]
        da=c(da,dam)
        nrm=i+dam*0
        nr=c(nr,nrm)
        TempF[[i]]=suppressWarnings(splinefun(data[[4*(i-1)+2]],Tem,method = "fmm")) 
        if (input$TrigTf==T) {
          ind=which(Tem>input$TrigTfixed)
          ti_trigg[i]=tim[ind[1]]
          A_triggm=0*dam-log(Kt)
          A_triggm[ind]=0 }
        else {
          ti_trigg[i]=tim[1]
          A_triggm=0*dam }
        A_trigg=c(A_trigg,A_triggm)
        wghDAm=dam*0+1/(max(dam)-min(dam))
        wghDA=c(wghDA,wghDAm) } }
    
    Te1=1/Te/8.31446
    thin_al_min=isolate(input$Trunc_al_min)
    thin_al_max=isolate(input$Trunc_al_max)  
    ti_thin=ti[al>thin_al_min&al<thin_al_max]
    al_thin=al[al>thin_al_min&al<thin_al_max]
    da_thin=da[al>thin_al_min&al<thin_al_max]
    Te_thin=Te[al>thin_al_min&al<thin_al_max]
    Te1_thin=Te1[al>thin_al_min&al<thin_al_max]
    wghDA_thin=wghDA[al>thin_al_min&al<thin_al_max]
    A_trigg_thin=A_trigg[al>thin_al_min&al<thin_al_max]
    nr_thin=nr[al>thin_al_min&al<thin_al_max]
    
    #write.table(cbind(ti_thin,al_thin,da_thin,Te_thin,Te1_thin,wghDA_thin,nr_thin), file="aaa.txt")
    
    list(ti,al,da,Te,Te1,TempF,Nruns,Npts,wghDA,nr,
         ti_thin,al_thin,da_thin,Te_thin,Te1_thin,wghDA_thin,nr_thin,ti_trigg,A_trigg,A_trigg_thin) #add column with N_run, prepare thin-out data!
  }
  
  output$Plot <- renderPlot(  height= function(){ ifelse(input$Compact==T, 350, input$BaseHeight)   },
    { make_plot1() })
  
  output$download1 <- downloadHandler(
    filename = function() {'plot1.emf'},
    content = function(file) {
      emf(file, emfPlus = FALSE, height= 4.67, width = 6.23)
      print(make_plot1())
      dev.off()
    })
  
  make_plot1 <- function(){    
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    res=treatData(inFile$datapath) #returns list (ti,al,da,Te,Te1,0,Nruns,Npts,WghDA)
    Nall=length(res[[1]])
    Nruns=res[[7]]
    Npts=res[[8]]
    nr=res[[10]]
    nr_thin=res[[17]]
    flag=character(Nall)
    flag[1:Nall]="d"
    if (input$dataTp==2) yExp=res[[2]]
    else yExp=res[[3]]
    if ((input$normC==T)&&(input$dataTp!=2)) yToPlot=yExp*res[[9]]
    else yToPlot=yExp
    if(input$drawX==2) { 
        xf=res[[4]]
        xL=expression('Temperature, K') 
        xf_thin=res[[14]] }
    if(input$drawX==1) { 
      xf=res[[4]]-273.15
      xL=expression('Temperature, °C') 
      xf_thin=res[[14]]-273.15 }
    if(input$drawX==3) { 
      xf=res[[1]]/60
      xL=expression('Time, min') 
      xf_thin=res[[11]]/60 }
    if(input$drawX==4) { 
      xf=res[[1]]
      xL=expression('Time, s') 
      xf_thin=res[[11]] }
     method=input$methodIni   #-----------------------used only for initual guess plotting!
     amin=input$al_min  #-----------------------------previously al_min_global,al_max_global were used here!
     amax=input$al_max
    if (input$taskType==1) {
      Ea1=ifelse(input$Ea1f==T,input$Ea1fixed,input$Ea1)
      n1=ifelse(input$nf1==T,input$n1fixed,input$n1)
      m1=ifelse(input$mf1==T,input$m1fixed,input$m1)
      if (method==1)
          Yf=funK1s(input$lnA1,Ea1,n1,m1,
               res[[5]],res[[1]],res[[2]],nr,input$dataTp,res[[9]],method,input$q,
               amin,amax,res[[19]])
      else Yf=funK1(input$lnA1,Ea1,n1,m1,
               res[[5]],res[[1]],res[[2]],nr,input$dataTp,res[[9]],method,input$q,
               amin,amax, res[[7]],res[[6]],input$suppression,res[[18]])
      }
    if (input$taskType==2) {
      Ea1=ifelse(input$Ea1f==T,input$Ea1fixed,input$Ea1)
      if (method==1)
          Yf=funK2s(input$lnA1,Ea1,input$reMod,
               res[[5]],res[[1]],res[[2]],nr,input$dataTp,res[[9]],method,amin,amax,res[[19]])
      else Yf=funK2(input$lnA1,Ea1,input$reMod,
                     res[[5]],res[[1]],res[[2]],nr,input$dataTp,res[[9]],method,amin,amax, res[[7]],res[[6]],input$suppression,res[[18]])
      }
    if (input$taskType==3) {
      Ea1=ifelse(input$Ea1f==T,input$Ea1fixed,input$Ea1)
      lnA22=ifelse(input$lnA22f==T,input$lnA22,0)
      Ea2=ifelse(input$Ea2f==T,input$Ea2fixed,input$Ea2)
      mu=ifelse(input$muf==T,input$mufixed,input$mu)
      if (method==1)
          Yf=funK3s(input$lnA1,Ea1,mu,input$lnA2,lnA22,Ea2,
               res[[5]],res[[1]],res[[2]],flag,nr,input$dataTp,res[[9]],method,amin,amax,res[[19]]) 
      else Yf=funK3(input$lnA1,Ea1,mu,input$lnA2,lnA22,Ea2,
                     res[[5]],res[[1]],res[[2]],flag,nr,input$dataTp,res[[9]],method,amin,amax, res[[7]],res[[6]],input$suppression,res[[18]])
      }
    if (input$taskType >= 4) {
        Ea1=ifelse(input$Ea1f==T,input$Ea1fixed,input$Ea1)
        Ea2=ifelse(input$Ea2f==T,input$Ea2fixed,input$Ea2)
        n1=ifelse(input$nf1==T,input$n1fixed,input$n1)
        m1=ifelse(input$mf1==T,input$m1fixed,input$m1)
        n2=ifelse(input$nf2==T,input$n2fixed,input$n2)
        m2=ifelse(input$mf2==T,input$m2fixed,input$m2)
        q2=ifelse(input$q2f==T,input$q2,input$q)
        if (input$taskType == 4) {
          if (method==1)
            Yf=funK4s(input$lnA1,Ea1,n1,m1,input$lnA2,Ea2,n2,m2,
               res[[5]],res[[1]],res[[2]],nr,input$dataTp,res[[9]],method,input$q,q2,amin,amax,res[[19]])
          else Yf=funK4(input$lnA1,Ea1,n1,m1,input$lnA2,Ea2,n2,m2,
                         res[[5]],res[[1]],res[[2]],nr,input$dataTp,res[[9]],method,input$q,q2,amin,amax, res[[7]],res[[6]],input$suppression,res[[18]])
        }
        if (input$taskType==5) {
          cd1=ifelse(input$cd1f==T,input$cd1fixed,input$cd1)
          cd2=ifelse(input$datasetsCnt==2,input$cd2,0)
          if (method==1)
              Yf=funK5s(input$lnA1,Ea1,n1,m1,input$lnA2,Ea2,n2,m2,cd1,cd2,
                   res[[5]],res[[1]],res[[2]],flag,nr,input$dataTp,res[[9]],method,input$q,q2,amin,amax,res[[19]]) 
          else Yf=funK5(input$lnA1,Ea1,n1,m1,input$lnA2,Ea2,n2,m2,cd1,cd2,
                         res[[5]],res[[1]],res[[2]],flag,nr,input$dataTp,res[[9]],method,input$q,q2,amin,amax, res[[7]],res[[6]],input$suppression,res[[18]])
          }
        if (input$taskType==6) {
            cd1=ifelse(input$cd1f==T,input$cd1fixed,input$cd1)
            cd2=ifelse(input$datasetsCnt==2,input$cd2,0)
            if (method==1)
                Yf=funK6s(input$lnA1,Ea1,n1,m1,input$lnA2,Ea2,n2,m2,cd1,cd2,
                    res[[5]],res[[1]],res[[2]],flag,nr,input$dataTp,res[[9]],method,input$q,q2,amin,amax,res[[19]]) 
            else Yf=funK6(input$lnA1,Ea1,n1,m1,input$lnA2,Ea2,n2,m2,cd1,cd2,
                           res[[5]],res[[1]],res[[2]],flag,nr,input$dataTp,res[[9]],method,input$q,q2,amin,amax, res[[7]],res[[6]],input$suppression,res[[18]])
            }
        if (input$taskType==7) {
            swA=ifelse(input$cd1f==T,input$cd1fixed,input$cd1)
            if (method==1)
              Yf=funK7s(input$lnA1,Ea1,n1,m1,input$lnA2,Ea2,n2,m2,swA,
                        res[[5]],res[[1]],res[[2]],flag,nr,input$dataTp,res[[9]],method,input$q,q2,amin,amax)
            else Yf=funK7(input$lnA1,Ea1,n1,m1,input$lnA2,Ea2,n2,m2,swA,
                   res[[5]],res[[1]],res[[2]],nr,input$dataTp,res[[9]],method,input$q,q2,amin,amax, res[[7]],res[[6]],input$suppression)
            }
        }
     if (input$taskType==8) {
       Ea1=ifelse(input$Ea1f==T,input$Ea1fixed,input$Ea1)
       Ea2=ifelse(input$Ea2f==T,input$Ea2fixed,input$Ea2)
       Ea3=ifelse(input$Ea3f==T,input$Ea3fixed,input$Ea3)
       lnA3=ifelse(input$lnA3f==T,input$lnA3fixed,input$lnA3)
       mu=ifelse(input$muf==T,input$mufixed,input$mu)
       n3=ifelse(input$nf3==T,input$n3fixed,input$n3)
       m3=ifelse(input$mf3==T,input$m3fixed,input$m3)
       cd1=ifelse(input$cd1f==T,input$cd1fixed,input$cd1)
       cd2=ifelse(input$datasetsCnt==2,input$cd2,0)
       if (method==1)
            Yf=funK8s(input$lnA1,Ea1,mu,input$lnA2,Ea2,lnA3,Ea3,n3,m3,cd1,cd2,
                   res[[5]],res[[1]],res[[2]],flag,nr,input$dataTp,res[[9]],method,input$q3,amin,amax,res[[19]]) 
       else Yf=funK8(input$lnA1,Ea1,mu,input$lnA2,Ea2,lnA3,Ea3,n3,m3,cd1,cd2,
                     res[[5]],res[[1]],res[[2]],flag,nr,input$dataTp,res[[9]],method,input$q3,amin,amax, res[[7]],res[[6]],input$suppression,res[[18]])
     }
     if (input$taskType==9) {
       Ea1=ifelse(input$Ea1f==T,input$Ea1fixed,input$Ea1)
       Ea2=ifelse(input$Ea2f==T,input$Ea2fixed,input$Ea2)
       Ea3=ifelse(input$Ea3f==T,input$Ea3fixed,input$Ea3)
       lnA1=ifelse(input$lnA1f==T,input$lnA1fixed,input$lnA1)
       lnA2=ifelse(input$lnA2f==T,input$lnA2fixed,input$lnA2)
       lnA3=ifelse(input$lnA3f==T,input$lnA3fixed,input$lnA3)
       n1=ifelse(input$nf1==T,input$n1fixed,input$n1)
       n2=ifelse(input$nf2==T,input$n2fixed,input$n2)
       m1=ifelse(input$mf1==T,input$m1fixed,input$m1)
       m2=ifelse(input$mf2==T,input$m2fixed,input$m2)
       n3=ifelse(input$nf3==T,input$n3fixed,input$n3)
       m3=ifelse(input$mf3==T,input$m3fixed,input$m3)
       cd1=ifelse(input$cd1f==T,input$cd1fixed,input$cd1)
       cd2=ifelse(input$datasetsCnt==2,input$cd2,0)
       q2=ifelse(input$q2f==T,input$q2,input$q)
       q3=ifelse(input$q3f==T,input$q3,input$q)
       Yf=funK9(lnA1,Ea1,n1,m1,lnA2,Ea2,n2,m2,lnA3,Ea3,n3,m3,cd1,cd2,input$cr21,input$cr22,
                     res[[5]],res[[1]],res[[2]],flag,nr,input$dataTp,res[[9]],method,input$q,q2,q3,amin,amax, res[[7]],res[[6]],input$suppression,res[[18]])
     }
     if (input$taskType==10) {
       Ea1=ifelse(input$Ea1f==T,input$Ea1fixed,input$Ea1)
       Ea2=ifelse(input$Ea2f==T,input$Ea2fixed,input$Ea2)
       Ea3=ifelse(input$Ea3f==T,input$Ea3fixed,input$Ea3)
       lnA1=ifelse(input$lnA1f==T,input$lnA1fixed,input$lnA1)
       lnA2=ifelse(input$lnA2f==T,input$lnA2fixed,input$lnA2)
       lnA3=ifelse(input$lnA3f==T,input$lnA3fixed,input$lnA3)
       n1=ifelse(input$nf1==T,input$n1fixed,input$n1)
       n2=ifelse(input$nf2==T,input$n2fixed,input$n2)
       m1=ifelse(input$mf1==T,input$m1fixed,input$m1)
       m2=ifelse(input$mf2==T,input$m2fixed,input$m2)
       n3=ifelse(input$nf3==T,input$n3fixed,input$n3)
       m3=ifelse(input$mf3==T,input$m3fixed,input$m3)
       cd1=ifelse(input$cd1f==T,input$cd1fixed,input$cd1)
       cd2=ifelse(input$datasetsCnt==2,input$cd2,0)
       q2=ifelse(input$q2f==T,input$q2,input$q)
       q3=ifelse(input$q3f==T,input$q3,input$q)
       Yf=funK10(lnA1,Ea1,n1,m1,lnA2,Ea2,n2,m2,lnA3,Ea3,n3,m3,cd1,cd2,input$swA,
                res[[5]],res[[1]],res[[2]],flag,nr,input$dataTp,res[[9]],method,input$q,q2,q3,amin,amax, res[[7]],res[[6]],input$suppression)
     }
    if ((input$normC==T)&&(input$dataTp!=2)) Yf=Yf*res[[9]]
    if (isolate(input$dataTp)!=3) { 
        maxYf=max(Yf)
        minYf=min(Yf) }
    else {
        maxYf=max(Yf/res[[9]])
        minYf=max(Yf/res[[9]]) }
    par (mar=c(4.1+input$marginBottom,5.1+input$marginLeft,3.1,1.1),
         mgp=c(3+input$marginTitle,1,0)) #par(oma=c(0,1,0,0))
    if (input$PointThin==1) xToPlot=xf
    else {
      if (ceiling(length(xf)/input$PointThin)>50) {
        xToPlot = xf[seq(1, length(xf), input$PointThin)]
        yToPlot = yToPlot[seq(1, length(yToPlot), input$PointThin)]
      }
    }
    plot(xToPlot,yToPlot, cex.lab=input$LegendSize, ylim = c(min(minYf,yToPlot), ifelse(input$maxAl==1,max(maxYf,yToPlot),input$maxAl)), xlab=xL, 
         ylab=ifelse(input$dataTp==2,expression('Conversion '*alpha),expression('d'*alpha*'/dt, 1/s')),
         cex.lab=input$AxisSize,cex.main=input$TitleSize, las=1, cex.axis=input$LegendSize,
         cex=input$PointSize, pch=input$PointType)
    #yTitle=ifelse(input$dataTp==2,expression('Conversion '*alpha),expression('d'*alpha*'/dt, 1/s'))
    #mtext(yTitle,
    #      side=2, line=3, cex=input$AxisSize)
    
    if (input$TrigTf == T) {
        if(input$drawX==1) abline (v=input$TrigTfixed-273.15, col="gray", lwd=2, lty=2)
        else abline (v=input$TrigTfixed, col="gray", lwd=2, lty=2)
        }
    if (input$iG==T) {
      if (isolate(input$dataTp)!=3) 
          for (i in 1:Nruns) lines(xf[nr==i],Yf[nr==i])
      else for (i in 1:Nruns) lines(xf[nr==i],Yf[nr==i]/res[[9]][nr==i]) }
    
    #write.table(cbind(xf,nr,yToPlot), file="initial_data1.txt") #to export the initial data of the model
    
    if (input$go!=0) {
      #nlrRes=reggr()[[1]]       #predict on full-data ??
      nlf=reggr()[[3]]  #predict(nlrRes)#[(Nruns*Npts+1):(Nruns*Npts*2)] #result - tga then dsc
      if (nlf[1] == "ERROR") return(NULL)
      if (isolate(input$datasetsCnt)==2) nlf=nlf[1:length(nr_thin)]
      if ((input$normC==T)&&(input$dataTp!=2)) nlf=nlf*res[[16]]
      #use array with indices for shift!
      if (isolate(input$dataTp)!=3)  
          for (i in 1:Nruns) lines(xf_thin[nr_thin==i],nlf[nr_thin==i],col='red', lwd=input$nlrLwd) 
                            #lines(xf[((i-1)*Npts+1):(i*Npts)],nlf[((i-1)*Npts+1):(i*Npts)],col='red', lwd=3) 
      else for (i in 1:Nruns) lines(xf_thin[nr_thin==i],nlf[nr_thin==i]/res[[16]][nr_thin==i],
                                    col='red', lwd=input$nlrLwd) 
      
      if (input$iG==T) legend('topleft', inset=c(0,-0.15), horiz=T, bty='n', xpd=NA, legend=c("experiment","initial guess","regression result"),
                              lty=c(NA,1,1), lwd=c(NA,1,input$nlrLwd),
                              pch=c(input$PointType,NA,NA), col=c('black','black','red'), cex=1.1) #input$LegendSize)
      else legend('topleft', inset=c(0,-0.15), bty='n', xpd=NA, horiz=T,
                  legend=c("experiment","regression result"),lty=c(NA,1), lwd=c(NA,input$nlrLwd),
                  pch=c(input$PointType,NA), col=c('black','red'), cex=1.1) #input$LegendSize)
      if (input$normC==T) nlf_out=nlf/res[[16]]
      else nlf_out=nlf
      #write.table(cbind(xf_thin,nlf_out,nr_thin), file="fk_result1.txt") #to export the optimized model data if needed
    }
    else legend('topleft', inset=c(0,-0.15), bty='n', xpd=NA, legend=c("experiment","initial guess"),lty=c(NA,1), 
                lwd=c(NA,1),pch=c(input$PointType,NA), col=c('black','black'), cex=1.1, horiz=T)
  }

  
output$Plot2 <- renderPlot(  height= function(){ ifelse(input$Compact==T, 350, input$BaseHeight)   },
                              { make_plot2() })
                                
  output$download2 <- downloadHandler(
      filename = function() {'plot2.emf'},
      content = function(file) {
      emf(file, emfPlus = FALSE, height= 4.67, width = 6.23)
      print(make_plot2())
      dev.off()
      })
                                
make_plot2 <- function(){
  inFile <- input$file2
  if (is.null(inFile)) return(NULL)
  res=treatData(inFile$datapath) #returns list (ti,al,da,Te,Te1,0,Nruns,Npts,WghDA)
  Nall=length(res[[1]])
  Nruns=res[[7]]
  Npts=res[[8]]
  nr=res[[10]]
  nr_thin=res[[17]]
  flag=character(Nall)
  flag[1:Nall]="d2"
  if (input$dataTp==2) yExp=res[[2]]
  else yExp=res[[3]]
  if ((input$normC==T)&&(input$dataTp!=2)) yToPlot=yExp*res[[9]]
  else yToPlot=yExp
  if(input$drawX2==2) { 
    xf=res[[4]]
    xL=expression('Temperature, K') 
    xf_thin=res[[14]] }
  else if(input$drawX2==1) { 
    xf=res[[4]]-273.15
    xL=expression('Temperature, °C') 
    xf_thin=res[[14]]-273.15 }
  else if(input$drawX2==3) { 
    xf=res[[1]]/60
    xL=expression('Time, min') 
    xf_thin=res[[11]]/60 }
  else if(input$drawX2==4) { 
    xf=res[[1]]
    xL=expression('Time, s') 
    xf_thin=res[[11]] }
  method=input$methodIni   #-----------------------used only for initual guess plotting!
  amin=input$al_min  #-----------------------------previously al_min_global,al_max_global were used here!
  amax=input$al_max
  if (input$taskType==1) {
    Ea1=ifelse(input$Ea1f==T,input$Ea1fixed,input$Ea1)
    n1=ifelse(input$nf1==T,input$n1fixed,input$n1)
    m1=ifelse(input$mf1==T,input$m1fixed,input$m1)
    if (method==1)
      Yf=funK1s(input$lnA1,Ea1,n1,m1,
                res[[5]],res[[1]],res[[2]],nr,input$dataTp,res[[9]],method,input$q,
                amin,amax,res[[19]])
    else Yf=funK1(input$lnA1,Ea1,n1,m1,
                  res[[5]],res[[1]],res[[2]],nr,input$dataTp,res[[9]],method,input$q,
                  amin,amax, res[[7]],res[[6]],input$suppression,res[[18]])
  }
  if (input$taskType==2) {
    Ea1=ifelse(input$Ea1f==T,input$Ea1fixed,input$Ea1)
    if (method==1)
      Yf=funK2s(input$lnA1,Ea1,input$reMod,
                res[[5]],res[[1]],res[[2]],nr,input$dataTp,res[[9]],method,amin,amax,res[[19]])
    else Yf=funK2(input$lnA1,Ea1,input$reMod,
                  res[[5]],res[[1]],res[[2]],nr,input$dataTp,res[[9]],method,amin,amax, res[[7]],res[[6]],input$suppression,res[[18]])
  }
  if (input$taskType==3) {
    Ea1=ifelse(input$Ea1f==T,input$Ea1fixed,input$Ea1)
    lnA22=ifelse(input$lnA22f==T,input$lnA22,0)
    Ea2=ifelse(input$Ea2f==T,input$Ea2fixed,input$Ea2)
    mu=ifelse(input$muf==T,input$mufixed,input$mu)
    if (method==1)
      Yf=funK3s(input$lnA1,Ea1,mu,input$lnA2,lnA22,Ea2,
                res[[5]],res[[1]],res[[2]],flag,nr,input$dataTp,res[[9]],method,amin,amax,res[[19]]) 
    else Yf=funK3(input$lnA1,Ea1,mu,input$lnA2,lnA22,Ea2,
                  res[[5]],res[[1]],res[[2]],flag,nr,input$dataTp,res[[9]],method,amin,amax, res[[7]],res[[6]],input$suppression,res[[18]])
  }
  if (input$taskType >= 4) {
    Ea1=ifelse(input$Ea1f==T,input$Ea1fixed,input$Ea1)
    Ea2=ifelse(input$Ea2f==T,input$Ea2fixed,input$Ea2)
    n1=ifelse(input$nf1==T,input$n1fixed,input$n1)
    m1=ifelse(input$mf1==T,input$m1fixed,input$m1)
    n2=ifelse(input$nf2==T,input$n2fixed,input$n2)
    m2=ifelse(input$mf2==T,input$m2fixed,input$m2)
    q2=ifelse(input$q2f==T,input$q2,input$q)
    if (input$taskType == 4) {
      if (method==1)
        Yf=funK4s(input$lnA1,Ea1,n1,m1,input$lnA2,Ea2,n2,m2,
                  res[[5]],res[[1]],res[[2]],nr,input$dataTp,res[[9]],method,input$q,q2,amin,amax,res[[19]])
      else Yf=funK4(input$lnA1,Ea1,n1,m1,input$lnA2,Ea2,n2,m2,
                    res[[5]],res[[1]],res[[2]],nr,input$dataTp,res[[9]],method,input$q,q2,amin,amax, res[[7]],res[[6]],input$suppression,res[[18]])
    }
    if (input$taskType==5) {
      cd1=ifelse(input$cd1f==T,input$cd1fixed,input$cd1)
      cd2=ifelse(input$datasetsCnt==2,input$cd2,0)
      if (method==1)
        Yf=funK5s(input$lnA1,Ea1,n1,m1,input$lnA2,Ea2,n2,m2,cd1,cd2,
                  res[[5]],res[[1]],res[[2]],flag,nr,input$dataTp,res[[9]],method,input$q,q2,amin,amax,res[[19]]) 
      else Yf=funK5(input$lnA1,Ea1,n1,m1,input$lnA2,Ea2,n2,m2,cd1,cd2,
                    res[[5]],res[[1]],res[[2]],flag,nr,input$dataTp,res[[9]],method,input$q,q2,amin,amax, res[[7]],res[[6]],input$suppression, res[[18]])
    }
    if (input$taskType==6) {
      cd1=ifelse(input$cd1f==T,input$cd1fixed,input$cd1)
      cd2=ifelse(input$datasetsCnt==2,input$cd2,0)
      if (method==1)
        Yf=funK6s(input$lnA1,Ea1,n1,m1,input$lnA2,Ea2,n2,m2,cd1,cd2,
                  res[[5]],res[[1]],res[[2]],flag,nr,input$dataTp,res[[9]],method,input$q,q2,amin,amax,res[[19]]) 
      else Yf=funK6(input$lnA1,Ea1,n1,m1,input$lnA2,Ea2,n2,m2,cd1,cd2,
                    res[[5]],res[[1]],res[[2]],flag,nr,input$dataTp,res[[9]],method,input$q,q2,amin,amax, res[[7]],res[[6]],input$suppression,res[[18]])
    }
    if (input$taskType==7) {
      swA=ifelse(input$cd1f==T,input$cd1fixed,input$cd1)
      if (method==1)
        Yf=funK7s(input$lnA1,Ea1,n1,m1,input$lnA2,Ea2,n2,m2,swA,
                  res[[5]],res[[1]],res[[2]],flag,nr,input$dataTp,res[[9]],method,input$q,q2,amin,amax)
      else Yf=funK7(input$lnA1,Ea1,n1,m1,input$lnA2,Ea2,n2,m2,swA,
                    res[[5]],res[[1]],res[[2]],nr,input$dataTp,res[[9]],method,input$q,q2,amin,amax, res[[7]],res[[6]],input$suppression)
      }
  }
  if (input$taskType==8) {
    Ea1=ifelse(input$Ea1f==T,input$Ea1fixed,input$Ea1)
    Ea2=ifelse(input$Ea2f==T,input$Ea2fixed,input$Ea2)
    Ea3=ifelse(input$Ea3f==T,input$Ea3fixed,input$Ea3)
    lnA3=ifelse(input$lnA3f==T,input$lnA3fixed,input$lnA3)
    mu=ifelse(input$muf==T,input$mufixed,input$mu)
    n3=ifelse(input$nf3==T,input$n3fixed,input$n3)
    m3=ifelse(input$mf3==T,input$m3fixed,input$m3)
    cd1=ifelse(input$cd1f==T,input$cd1fixed,input$cd1)
    cd2=ifelse(input$datasetsCnt==2,input$cd2,0)
    if (method==1)
      Yf=funK8s(input$lnA1,Ea1,mu,input$lnA2,Ea2,lnA3,Ea3,n3,m3,cd1,cd2,
                res[[5]],res[[1]],res[[2]],flag,nr,input$dataTp,res[[9]],method,input$q3,amin,amax,res[[19]]) 
    else Yf=funK8(input$lnA1,Ea1,mu,input$lnA2,Ea2,lnA3,Ea3,n3,m3,cd1,cd2,
                  res[[5]],res[[1]],res[[2]],flag,nr,input$dataTp,res[[9]],method,input$q3,amin,amax, res[[7]],res[[6]],input$suppression,res[[18]])
  }
  if (input$taskType==9) {
    Ea1=ifelse(input$Ea1f==T,input$Ea1fixed,input$Ea1)
    Ea2=ifelse(input$Ea2f==T,input$Ea2fixed,input$Ea2)
    Ea3=ifelse(input$Ea3f==T,input$Ea3fixed,input$Ea3)
    lnA1=ifelse(input$lnA1f==T,input$lnA1fixed,input$lnA1)
    lnA2=ifelse(input$lnA2f==T,input$lnA2fixed,input$lnA2)
    lnA3=ifelse(input$lnA3f==T,input$lnA3fixed,input$lnA3)
    n1=ifelse(input$nf1==T,input$n1fixed,input$n1)
    n2=ifelse(input$nf2==T,input$n2fixed,input$n2)
    m1=ifelse(input$mf1==T,input$m1fixed,input$m1)
    m2=ifelse(input$mf2==T,input$m2fixed,input$m2)
    n3=ifelse(input$nf3==T,input$n3fixed,input$n3)
    m3=ifelse(input$mf3==T,input$m3fixed,input$m3)
    cd1=ifelse(input$cd1f==T,input$cd1fixed,input$cd1)
    cd2=ifelse(input$datasetsCnt==2,input$cd2,0)
    q2=ifelse(input$q2f==T,input$q2,input$q)
    q3=ifelse(input$q3f==T,input$q3,input$q)
    Yf=funK9(lnA1,Ea1,n1,m1,lnA2,Ea2,n2,m2,lnA3,Ea3,n3,m3,cd1,cd2,input$cr21,input$cr22,
             res[[5]],res[[1]],res[[2]],flag,nr,input$dataTp,res[[9]],method,input$q,q2,q3,amin,amax, res[[7]],res[[6]],input$suppression,res[[18]])
  }
  if (input$taskType==10) {
    Ea1=ifelse(input$Ea1f==T,input$Ea1fixed,input$Ea1)
    Ea2=ifelse(input$Ea2f==T,input$Ea2fixed,input$Ea2)
    Ea3=ifelse(input$Ea3f==T,input$Ea3fixed,input$Ea3)
    lnA1=ifelse(input$lnA1f==T,input$lnA1fixed,input$lnA1)
    lnA2=ifelse(input$lnA2f==T,input$lnA2fixed,input$lnA2)
    lnA3=ifelse(input$lnA3f==T,input$lnA3fixed,input$lnA3)
    n1=ifelse(input$nf1==T,input$n1fixed,input$n1)
    n2=ifelse(input$nf2==T,input$n2fixed,input$n2)
    m1=ifelse(input$mf1==T,input$m1fixed,input$m1)
    m2=ifelse(input$mf2==T,input$m2fixed,input$m2)
    n3=ifelse(input$nf3==T,input$n3fixed,input$n3)
    m3=ifelse(input$mf3==T,input$m3fixed,input$m3)
    cd1=ifelse(input$cd1f==T,input$cd1fixed,input$cd1)
    cd2=ifelse(input$datasetsCnt==2,input$cd2,0)
    q2=ifelse(input$q2f==T,input$q2,input$q)
    q3=ifelse(input$q3f==T,input$q3,input$q)
    Yf=funK10(lnA1,Ea1,n1,m1,lnA2,Ea2,n2,m2,lnA3,Ea3,n3,m3,cd1,cd2,input$swA,
              res[[5]],res[[1]],res[[2]],flag,nr,input$dataTp,res[[9]],method,input$q,q2,q3,amin,amax, res[[7]],res[[6]],input$suppression)
  }
  if ((input$normC==T)&&(input$dataTp!=2)) Yf=Yf*res[[9]]
  if (isolate(input$dataTp)!=3) { 
    maxYf=max(Yf)
    minYf=min(Yf) }
  else {
    maxYf=max(Yf/res[[9]])
    minYf=max(Yf/res[[9]]) }

  par (mar=c(4.1+input$marginBottom,5.1+input$marginLeft,3.1,1.1),
       mgp=c(3+input$marginTitle,1,0)) #par(oma=c(0,1,0,0))
  if (input$PointThin==1) xToPlot=xf
  else {
    if (ceiling(length(xf)/input$PointThin)>50) {
      xToPlot = xf[seq(1, length(xf), input$PointThin)]
      yToPlot = yToPlot[seq(1, length(yToPlot), input$PointThin)]
    }
  }
  plot(xToPlot,yToPlot, cex.lab=input$LegendSize, ylim = c(min(minYf,yToPlot), max(maxYf,yToPlot)), xlab=xL, 
       ylab=ifelse(input$dataTp==2,expression('Conversion '*alpha),expression('d'*alpha*'/dt, 1/s')),
       cex.lab=input$AxisSize,cex.main=input$TitleSize, las=1, cex.axis=input$LegendSize,
       cex=input$PointSize, pch=input$PointType)
  
  if (input$TrigTf == T) {
    if(input$drawX==1) abline (v=input$TrigTfixed-273.15, col="gray", lwd=2, lty=2)
    else abline (v=input$TrigTfixed, col="gray", lwd=2, lty=2)
  }
  if (input$iG==T) {
    if (isolate(input$dataTp)!=3) 
      for (i in 1:Nruns) lines(xf[nr==i],Yf[nr==i])
    else for (i in 1:Nruns) lines(xf[nr==i],Yf[nr==i]/res[[9]][nr==i]) }
  
  #write.table(cbind(xf,nr,yToPlot), file="initial_data2.txt") #to export the initial data of the model
  
  if (input$go!=0) {
    #nlrRes=reggr()[[1]]       #predict on full-data ??
    nlf=reggr()[[3]]  #predict(nlrRes)#[(Nruns*Npts+1):(Nruns*Npts*2)] #result - tga then dsc
    if (nlf[1] == "ERROR") return(NULL)
    nlf=nlf[(length(nlf)-length(nr_thin)+1):length(nlf)]  #---different from plot1!
    if ((input$normC==T)&&(input$dataTp!=2)) nlf=nlf*res[[16]]
    #use array with indices for shift!
    if (isolate(input$dataTp)!=3)  
      for (i in 1:Nruns) lines(xf_thin[nr_thin==i],nlf[nr_thin==i],col='red', lwd=input$nlrLwd) 
    #lines(xf[((i-1)*Npts+1):(i*Npts)],nlf[((i-1)*Npts+1):(i*Npts)],col='red', lwd=3) 
    else for (i in 1:Nruns) lines(xf_thin[nr_thin==i],nlf[nr_thin==i]/res[[16]][nr_thin==i],
                                  col='red', lwd=input$nlrLwd) 
    
    if (input$iG==T) legend('topleft', inset=c(0,-0.15), horiz=T, bty='n', xpd=NA, legend=c("experiment","initial guess","regression result"),
                            lty=c(NA,1,1), lwd=c(NA,1,input$nlrLwd),
                            pch=c(input$PointType,NA,NA), col=c('black','black','red'), cex=1.1) #input$LegendSize)
    else legend('topleft', inset=c(0,-0.15), bty='n', xpd=NA, horiz=T,
                legend=c("experiment","regression result"),lty=c(NA,1), lwd=c(NA,input$nlrLwd),
                pch=c(input$PointType,NA), col=c('black','red'), cex=1.1)
    
    #write.table(cbind(xf_thin,nlf/res[[16]],nr_thin), file="fk_result2.txt") #to export the optimized model data if needed
  }
  else legend('topleft', inset=c(0,-0.15), bty='n', xpd=NA, legend=c("experiment","initial guess"),lty=c(NA,1), 
              lwd=c(NA,1),pch=c(input$PointType,NA), col=c('black','black'), cex=1.1, horiz=T)
}
    
  reggr=eventReactive(input$go, {
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    #input$go
    shinyjs::disable("go")
    resD=treatData(inFile$datapath) #returns list (ti,al,da,Te,Te1,TempF,Nruns,Npts,wghDA,nrm,
                                    #                ti_thin,al_thin,da_thin,Te_thin,Te1_thin,wghDA_thin,nr_thin)
    Nall=length(resD[[11]])
    dat=resD[[13]]
    T1t=resD[[15]]
    at=resD[[12]]
    tt=resD[[11]]
    Tt=resD[[14]]
    wgh0=numeric(Nall)
    wghDA=resD[[16]] # or first resT[[9]]
    flag=character(Nall)
    flag[1:Nall]="d"
    nr=resD[[17]] 
    ti_trigg=resD[[18]]
    A_trigg=resD[[20]]
    TempF=resD[[6]]
    Nruns=resD[[7]]
    beginning = Sys.time()
    if (isolate(input$datasetsCnt)==2) {
      inFile2 <- input$file2
      res2=treatData(inFile2$datapath) 
      Nall2=length(res2[[11]])
      dat=c(dat,res2[[13]])
      T1t=c(T1t,res2[[15]])
      at=c(at,res2[[12]])
      tt=c(tt,res2[[11]])
      Tt=c(Tt,res2[[14]])
      wgh02=numeric(Nall2)
      wghDA=c(wghDA,res2[[16]]) # or first resT[[9]]
      flag2=character(Nall2)
      flag2[1:Nall2]="d2"
      flag=c(flag,flag2)
      nr2=res2[[17]]+nr[Nall]
      ti_trigg2=res2[[18]]
      A_trigg2=res2[[20]]
      nr=c(nr,nr2)
      ti_trigg=c(ti_trigg,ti_trigg2)
      A_trigg=c(A_trigg,A_trigg2)
      TempF=c(TempF,res2[[6]])
      Nruns=Nruns+res2[[7]]
    }
    Ea_limits=c(50,400) #-----------global limits for the kinetic parametes (if change it, change also input limits)
    lnA_limits=c(5,120)
    n_limits=c(0,3)
    if (isolate(input$m1fixed)==10) n_limits=c(0.5,6) #-------specific limits for KJMAE equation
    else n_limits=c(0,3)
    n_ePT_limits=c(0,3)
    m_limits=c(-1,3)
    Dtype=isolate(input$dataTp)
    method=isolate(input$methodDiff)
    optimType=isolate(input$optimType)
    suppression=isolate(input$suppression)
    q=isolate(input$q)
    amin=isolate(input$al_min)
    amax=isolate(input$al_max)
    if (Dtype==3) nlY=dat*wghDA
    else if (Dtype==1) nlY=dat
    else nlY=at

    if (isolate(input$taskType)==1) {
      startList=isolate(list(lnA1=input$lnA1,Ea1=input$Ea1,n1=input$n1,m1=input$m1))
      lowerList=c(lnA_limits[1],Ea_limits[1],n_limits[1],m_limits[1])
      upperList=c(lnA_limits[2],Ea_limits[2],n_limits[2],m_limits[2])
      nInd=0
      if (isolate(input$Ea1f)==T) {
        Ea1=isolate(input$Ea1fixed)
        indF=2
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$nf1)==T) {
        n1=isolate(input$n1fixed)
        indF=3-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$mf1)==T) {
        m1=isolate(input$m1fixed)
        indF=4-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] }
      if (method>1) {
        if (optimType==1) 
          nlrRes=try(nlsLM(nlY ~ funK1(lnA1,Ea1,n1,m1,T1t,tt,at,nr,Dtype,wghDA,method,q,amin,amax,Nruns,TempF,suppression,ti_trigg), 
                     start=startList,lower=lowerList, upper=upperList,
                     control= nls.lm.control(maxiter=1000)))
        else if (optimType==2) 
          nlrRes=try(nls(nlY ~ funK1(lnA1,Ea1,n1,m1,T1t,tt,at,nr,Dtype,wghDA,method,q,amin,amax,Nruns,TempF,suppression,ti_trigg), 
                         start=startList,algorithm = "port",lower=lowerList, upper=upperList,
                         control= nls.control(maxiter=1000, warnOnly=T)))  }
      else nlrRes=try(nlsLM(nlY ~ funK1s(lnA1,Ea1,n1,m1,T1t,tt,at,nr,Dtype,wghDA,method,q,amin,amax,A_trigg), 
                        start=startList,lower=lowerList, upper=upperList,
                        control= nls.lm.control(maxiter=1000)))
    }
    if (isolate(input$taskType)==2) {
      startList=isolate(list(lnA1=input$lnA1,Ea1=input$Ea1))
      RType=isolate(input$reMod)
      lowerList=c(lnA_limits[1],Ea_limits[1])
      upperList=c(lnA_limits[2],Ea_limits[2])
      nInd=0
      if (isolate(input$Ea1f)==T) {
        Ea1=isolate(input$Ea1fixed)
        indF=2
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (method>1) {
        if (optimType==1) 
          nlrRes=try(nlsLM(nlY ~ funK2(lnA1,Ea1,RType,T1t,tt,at,nr,Dtype,wghDA,method,amin,amax,Nruns,TempF,suppression,ti_trigg), 
                       start=startList,lower=lowerList, upper=upperList,
                       control= nls.lm.control(maxiter=1000)))
        else if (optimType==2)
          nlrRes=try(nls(nlY ~ funK2(lnA1,Ea1,RType,T1t,tt,at,nr,Dtype,wghDA,method,amin,amax,Nruns,TempF,suppression,ti_trigg), 
                       start=startList,algorithm = "port",lower=lowerList, upper=upperList,
                       control= nls.control(maxiter=1000, warnOnly=T))) }
      else nlrRes=try(nlsLM(nlY ~ funK2s(lnA1,Ea1,RType,T1t,tt,at,nr,Dtype,wghDA,method,amin,amax,A_trigg), 
                        start=startList,lower=lowerList, upper=upperList,
                        control= nls.lm.control(maxiter=1000)))
    }
    if (isolate(input$taskType)==3) { #lnA1,Ea1,mu,lnA2,lnA22,Ea2,fT1,ft,fa,flag,nr,Dtype,WghDA,method,amin,amax
      lnA22val=ifelse(isolate(input$lnA22f)==T,isolate(input$lnA22),0)
      startList=isolate(list(lnA1=input$lnA1,Ea1=input$Ea1,mu=input$mu,lnA2=input$lnA2,lnA22=lnA22val,Ea2=input$Ea2))
      lowerList=c(lnA_limits[1],Ea_limits[1],0,lnA_limits[1],lnA_limits[1],Ea_limits[1])
      upperList=c(lnA_limits[2],Ea_limits[2],1,lnA_limits[2],lnA_limits[2],Ea_limits[2])
      nInd=0
      if (isolate(input$Ea1f)==T) {
        Ea1=isolate(input$Ea1fixed)
        indF=2
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$muf)==T) {
        mu=isolate(input$mufixed)
        indF=3-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF]
        nInd=nInd+1 }
      if (isolate(input$lnA22f)==F) {
        lnA22=0
        indF=5-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] }
      if (isolate(input$Ea2f)==T) {
        Ea2=isolate(input$Ea2fixed)
        indF=6-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] }
      
      if (method>1) {
        if (optimType==1) 
            nlrRes=try(nlsLM(nlY ~ funK3(lnA1,Ea1,mu,lnA2,lnA22,Ea2,T1t,tt,at,flag,nr,Dtype,wghDA,method,amin,amax,Nruns,TempF,suppression,ti_trigg), 
                     start=startList,lower=lowerList, upper=upperList,
                     control= nls.lm.control(maxiter=1000)))
    else if (optimType==2) 
            nlrRes=try(nls(nlY ~ funK3(lnA1,Ea1,mu,lnA2,lnA22,Ea2,T1t,tt,at,flag,nr,Dtype,wghDA,method,amin,amax,Nruns,TempF,suppression,ti_trigg), 
                          start=startList,algorithm = "port",lower=lowerList, upper=upperList,
                          control= nls.control(maxiter=1000, warnOnly=T)))
#    else if (optimType==3) {
#            myfun="nlY ~ funK3(lnA1,Ea1,mu,lnA2,Ea2,T1t,tt,at,nr,Dtype,wghDA,method,amin,amax,Nruns,TempF,suppression)"
#            myres=model2resfun(myfun, startList)
#            nlrRes=nlfb(startList, resfn=myres, lower=lowerList, upper=upperList) }
#    else if (optimType==4) #{
#            nlrRes=nlxb(nlY ~ funK3(lnA1,Ea1,mu,lnA2,Ea2,T1t,tt,at,nr,Dtype,wghDA,method,amin,amax,Nruns,TempF,suppression), 
#                    start=startList, lower=lowerList, upper=upperList)
#          
#          min.RSS <- function(lnA1,Ea1,mu,lnA2,Ea2) {
#            with(as.list(c(lnA1,Ea1,mu,lnA2,Ea2)),{
#              sse = sum((nlY - funK3(lnA1,Ea1,mu,lnA2,Ea2,T1t,tt,at,nr,Dtype,wghDA,method,amin,amax,Nruns,TempF,suppression))^2)
#              return(sse)
#            }) }
#            nlrRes=nlsLM(min.RSS(lnA1,Ea1,mu,lnA2,Ea2), 
#                       start=startList,lower=lowerList, upper=upperList,
#                       control= nls.lm.control(maxiter=1000)) 
#            strFunk="funK3(lnA1,Ea1,mu,lnA2,Ea2,T1t,tt,at,nr,Dtype,wghDA,method,amin,amax,Nruns,TempF,suppression)"
#            for (i in 1:length(coef(nlrRes))) {
#                   strFunk=gsub(names(coef(nlrRes))[i],coef(nlrRes)[i],strFunk) }
#            predictY=eval(parse(text=strFunk))
#        } else if (optimType==4) {
#          nlrRes=nls(sum((nlY - funK3(lnA1,Ea1,mu,lnA2,Ea2,T1t,tt,at,nr,Dtype,wghDA,method,amin,amax,Nruns,TempF,suppression))^2), 
#                                          start=startList,algorithm = "port",lower=lowerList, upper=upperList,
#                                          control= nls.lm.control(maxiter=1000)) 
        }
      else 
        nlrRes=try(nlsLM(nlY ~ funK3s(lnA1,Ea1,mu,lnA2,lnA22,Ea2,T1t,tt,at,flag,nr,Dtype,wghDA,method,amin,amax,A_trigg), 
                     start=startList,lower=lowerList, upper=upperList,
                     control= nls.lm.control(maxiter=1000)))
    }
    if (isolate(input$taskType)==4) {
      startList=isolate(list(lnA1=input$lnA1,Ea1=input$Ea1,n1=input$n1,m1=input$m1,
                             lnA2=input$lnA2,Ea2=input$Ea2,n2=input$n2,m2=input$m2))
      lowerList=c(lnA_limits[1],Ea_limits[1],n_limits[1],m_limits[1],lnA_limits[1],Ea_limits[1],n_limits[1],m_limits[1]) 
      upperList=c(lnA_limits[2],Ea_limits[2],n_limits[2],m_limits[2],lnA_limits[2],Ea_limits[2],n_limits[2],m_limits[2]) 
      nInd=0
      if (isolate(input$Ea1f)==T) {
        Ea1=isolate(input$Ea1fixed)
        indF=2-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$nf1)==T) {
        n1=isolate(input$n1fixed)
        indF=3-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$mf1)==T) {
        m1=isolate(input$m1fixed)
        indF=4-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$Ea2f)==T) {
        Ea2=isolate(input$Ea2fixed)
        indF=6-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$nf2)==T) {
        n2=isolate(input$n2fixed)
        indF=7-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$mf2)==T) {
        m2=isolate(input$m2fixed)
        indF=8-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      q2=ifelse(isolate(input$q2f)==T,isolate(input$q2),isolate(input$q)) 
      if (method>1) {
        if (optimType==1) 
          nlrRes=try(nlsLM(nlY ~ funK4(lnA1,Ea1,n1,m1,lnA2,Ea2,n2,m2,T1t,tt,at,nr,Dtype,wghDA,method,q,q2,amin,amax,Nruns,TempF,suppression,ti_trigg), 
                     start=startList,lower=lowerList, upper=upperList,control= nls.lm.control(maxiter=1000,maxfev=10000)))
        else if (optimType==2) 
          nlrRes=try(nls(nlY ~ funK4(lnA1,Ea1,n1,m1,lnA2,Ea2,n2,m2,T1t,tt,at,nr,Dtype,wghDA,method,q,q2,amin,amax,Nruns,TempF,suppression,ti_trigg), 
                         start=startList,algorithm = "port",lower=lowerList, upper=upperList,
                         control= nls.control(maxiter=1000, warnOnly=T))) }
      else
          nlrRes=try(nlsLM(nlY ~ funK4s(lnA1,Ea1,n1,m1,lnA2,Ea2,n2,m2,T1t,tt,at,nr,Dtype,wghDA,method,q,q2,amin,amax,A_trigg), 
                     start=startList,lower=lowerList, upper=upperList,control= nls.lm.control(maxiter=1000,maxfev=10000)))
    }

    if (isolate(input$taskType)>=5) {
      cd2val=ifelse(isolate(input$datasetsCnt)==2,isolate(input$cd2),0)
      startList=isolate(list(lnA1=input$lnA1,Ea1=input$Ea1,n1=input$n1,m1=input$m1,
                             lnA2=input$lnA2,Ea2=input$Ea2,n2=input$n2,m2=input$m2,cd1=input$cd1,cd2=cd2val))
      cd1min=isolate(input$cd1Lim[1])
      cd1max=isolate(input$cd1Lim[2])
      lowerList=c(lnA_limits[1],Ea_limits[1],n_limits[1],-1,lnA_limits[1],Ea_limits[1],0,-1,cd1min,0.01)
      upperList=c(lnA_limits[2],Ea_limits[2],n_limits[2],3,lnA_limits[2],Ea_limits[2],3,3,cd1max,0.99)
      nInd=0
      if (isolate(input$Ea1f)==T) {
        Ea1=isolate(input$Ea1fixed)
        indF=2-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$nf1)==T) {
        n1=isolate(input$n1fixed)
        indF=3-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$mf1)==T) {
        m1=isolate(input$m1fixed)
        indF=4-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$Ea2f)==T) {
        Ea2=isolate(input$Ea2fixed)
        indF=6-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$nf2)==T) {
        n2=isolate(input$n2fixed)
        indF=7-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$mf2)==T) {
        m2=isolate(input$m2fixed)
        indF=8-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$cd1f)==T) {
        cd1=isolate(input$cd1fixed)
        indF=9-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$datasetsCnt)==1) {
        cd2=0
        indF=10-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] }
      q2=ifelse(isolate(input$q2f)==T,isolate(input$q2),isolate(input$q))
      if (isolate(input$taskType)==5) {
        if (method>1) {
          if (optimType==1) 
            nlrRes=try(nlsLM(nlY ~ funK5(lnA1,Ea1,n1,m1,lnA2,Ea2,n2,m2,cd1,cd2,T1t,tt,at,flag,nr,Dtype,wghDA,method,q,q2,amin,amax,Nruns,TempF,suppression,ti_trigg), 
                             start=startList,lower=lowerList, upper=upperList,control= nls.lm.control(maxiter=1000,maxfev=10000)))
          else if (optimType==2) 
            nlrRes=try(nls(nlY ~ funK5(lnA1,Ea1,n1,m1,lnA2,Ea2,n2,m2,cd1,cd2,T1t,tt,at,flag,nr,Dtype,wghDA,method,q,q2,amin,amax,Nruns,TempF,suppression,ti_trigg), 
                           start=startList,algorithm = "port",lower=lowerList, upper=upperList,
                           control= nls.control(maxiter=1000, warnOnly=T))) }
        else
          nlrRes=try(nlsLM(nlY ~ funK5s(lnA1,Ea1,n1,m1,lnA2,Ea2,n2,m2,cd1,cd2,T1t,tt,at,flag,nr,Dtype,wghDA,method,q,q2,amin,amax,A_trigg), 
                           start=startList,lower=lowerList, upper=upperList,control= nls.lm.control(maxiter=1000,maxfev=10000)))
      }
      if (isolate(input$taskType)==6) {
        if (method>1) {
          if (optimType==1) 
            nlrRes=try(nlsLM(nlY ~ funK6(lnA1,Ea1,n1,m1,lnA2,Ea2,n2,m2,cd1,cd2,T1t,tt,at,flag,nr,Dtype,wghDA,method,q,q2,amin,amax,Nruns,TempF,suppression,ti_trigg), 
                             start=startList,lower=lowerList, upper=upperList,control= nls.lm.control(maxiter=1000,maxfev=10000)))
          else if (optimType==2) 
            nlrRes=try(nls(nlY ~ funK6(lnA1,Ea1,n1,m1,lnA2,Ea2,n2,m2,cd1,cd2,T1t,tt,at,flag,nr,Dtype,wghDA,method,q,q2,amin,amax,Nruns,TempF,suppression,ti_trigg), 
                           start=startList,algorithm = "port",lower=lowerList, upper=upperList,
                           control= nls.control(maxiter=1000, warnOnly=T))) }
        else
          nlrRes=try(nlsLM(nlY ~ funK6s(lnA1,Ea1,n1,m1,lnA2,Ea2,n2,m2,cd1,cd2,T1t,tt,at,flag,nr,Dtype,wghDA,method,q,q2,amin,amax,A_trigg), 
                           start=startList,lower=lowerList, upper=upperList,control= nls.lm.control(maxiter=1000,maxfev=10000)))
      }
    }
    
    
    if (isolate(input$taskType)==7) {
      startList=isolate(list(lnA1=input$lnA1,Ea1=input$Ea1,n1=input$n1,m1=input$m1,
                             lnA2=input$lnA2,Ea2=input$Ea2,n2=input$n2,m2=input$m2,swA=input$cd1))
      cd1min=isolate(input$cd1Lim[1])
      cd1max=isolate(input$cd1Lim[2])
      lowerList=c(lnA_limits[1],Ea_limits[1],n_limits[1],m_limits[1],lnA_limits[1],Ea_limits[1],n_limits[1],m_limits[1],cd1min)
      upperList=c(lnA_limits[2],Ea_limits[2],n_limits[2],m_limits[2],lnA_limits[2],Ea_limits[2],n_limits[2],m_limits[2],cd1max)
      nInd=0
      if (isolate(input$Ea1f)==T) {
        Ea1=isolate(input$Ea1fixed)
        indF=2-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$nf1)==T) {
        n1=isolate(input$n1fixed)
        indF=3-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$mf1)==T) {
        m1=isolate(input$m1fixed)
        indF=4-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$Ea2f)==T) {
        Ea2=isolate(input$Ea2fixed)
        indF=6-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$nf2)==T) {
        n2=isolate(input$n2fixed)
        indF=7-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$mf2)==T) {
        m2=isolate(input$m2fixed)
        indF=8-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$cd1f)==T) {
        swA=isolate(input$cd1fixed)
        indF=9-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      q2=ifelse(isolate(input$q2f)==T,isolate(input$q2),isolate(input$q))
      if (method>1) {
        if (optimType==1) 
          nlrRes=try(nlsLM(nlY ~ funK7(lnA1,Ea1,n1,m1,lnA2,Ea2,n2,m2,swA,T1t,tt,at,nr,Dtype,wghDA,method,q,q2,amin,amax,Nruns,TempF,suppression),
                           start=startList,lower=lowerList, upper=upperList,control= nls.lm.control(maxiter=1000,maxfev=10000)))
        else if (optimType==2) 
          nlrRes=try(nls(nlY ~ funK7(lnA1,Ea1,n1,m1,lnA2,Ea2,n2,m2,swA,T1t,tt,at,nr,Dtype,wghDA,method,q,q2,amin,amax,Nruns,TempF,suppression),
                         start=startList,algorithm = "port",lower=lowerList, upper=upperList,
                         control= nls.control(maxiter=1000, warnOnly=T))) }
      else
        nlrRes=try(nlsLM(nlY ~ funK7s(lnA1,Ea1,n1,m1,lnA2,Ea2,n2,m2,swA,T1t,tt,at,flag,nr,Dtype,wghDA,method,q,q2,amin,amax), 
                         start=startList,lower=lowerList, upper=upperList,control= nls.lm.control(maxiter=1000,maxfev=10000)))
    }
    
    if (isolate(input$taskType)==8) {
      cd2val=ifelse(isolate(input$datasetsCnt)==2,isolate(input$cd2),0)
      startList=isolate(list(lnA1=input$lnA1,Ea1=input$Ea1,mu=input$mu,
                             lnA2=input$lnA2,Ea2=input$Ea2, lnA3=input$lnA3,Ea3=input$Ea3,
                             n3=input$n3,m3=input$m3,cd1=input$cd1,cd2=cd2val))
      cd1min=isolate(input$cd1Lim[1])
      cd1max=isolate(input$cd1Lim[2])
      #lnA_limits[1],Ea_limits[1],n_limits[1],m_limits[1]
      lowerList=c(lnA_limits[1],Ea_limits[1],0,lnA_limits[1],Ea_limits[1],lnA_limits[1],Ea_limits[1],
                  n_limits[1],m_limits[1],cd1min,0.01)
      upperList=c(lnA_limits[2],Ea_limits[2],1,lnA_limits[2],Ea_limits[2],lnA_limits[2],Ea_limits[2],
                  n_limits[2],m_limits[2],cd1max,0.99)
      nInd=0
      if (isolate(input$Ea1f)==T) {
        Ea1=isolate(input$Ea1fixed)
        indF=2-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$muf)==T) {
        mu=isolate(input$mufixed)
        indF=3-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$Ea2f)==T) {
        Ea2=isolate(input$Ea2fixed)
        indF=5-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$lnA3f)==T) {
        lnA3=isolate(input$lnA3fixed)
        indF=6-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$Ea3f)==T) {
        Ea3=isolate(input$Ea3fixed)
        indF=7-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$nf3)==T) {
        n3=isolate(input$n3fixed)
        indF=8-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$mf3)==T) {
        m3=isolate(input$m3fixed)
        indF=9-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$cd1f)==T) {
        cd1=isolate(input$cd1fixed)
        indF=10-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$datasetsCnt)==1) {
        cd2=0
        indF=11-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] }
      q3=isolate(input$q3)
      
        if (method>1) {
          if (optimType==1) 
            nlrRes=try(nlsLM(nlY ~ funK8(lnA1,Ea1,mu,lnA2,Ea2,lnA3,Ea3,n3,m3,cd1,cd2,T1t,tt,at,flag,nr,Dtype,wghDA,method,q3,amin,amax,Nruns,TempF,suppression,ti_trigg), 
                             start=startList,lower=lowerList, upper=upperList,control= nls.lm.control(maxiter=1000,maxfev=10000)))
          else if (optimType==2) 
            nlrRes=try(nls(nlY ~ funK8(lnA1,Ea1,mu,lnA2,Ea2,lnA3,Ea3,n3,m3,cd1,cd2,T1t,tt,at,flag,nr,Dtype,wghDA,method,q3,amin,amax,Nruns,TempF,suppression,ti_trigg), 
                           start=startList,algorithm = "port",lower=lowerList, upper=upperList,
                           control= nls.control(maxiter=1000, warnOnly=T))) }
        else
          nlrRes=try(nlsLM(nlY ~ funK8s(lnA1,Ea1,mu,lnA2,Ea2,lnA3,Ea3,n3,m3,cd1,cd2,T1t,tt,at,flag,nr,Dtype,wghDA,method,q3,amin,amax,A_trigg), 
                           start=startList,lower=lowerList, upper=upperList,control= nls.lm.control(maxiter=1000,maxfev=10000)))
      
    }

    if (isolate(input$taskType)==9) {
      cd2val=ifelse(isolate(input$datasetsCnt)==2,isolate(input$cd2),0)
      startList=isolate(list(lnA1=input$lnA1,Ea1=input$Ea1,n1=input$n1,m1=input$m1,
                             lnA2=input$lnA2,Ea2=input$Ea2, n2=input$n2,m2=input$m2,
                             lnA3=input$lnA3,Ea3=input$Ea3, n3=input$n3,m3=input$m3,
                             cd1=input$cd1,cd2=cd2val, cr21=input$cr21, cr22=input$cr22))
      cd1min=isolate(input$cd1Lim[1])
      cd1max=isolate(input$cd1Lim[2])
      cr_limits=c(0.1,10)
      lowerList=c(lnA_limits[1],Ea_limits[1],n_limits[1],m_limits[1],
                  lnA_limits[1],Ea_limits[1],n_ePT_limits[1],m_limits[1],
                  lnA_limits[1],Ea_limits[1],n_ePT_limits[1],m_limits[1],
                  cd1min,cd1min,cr_limits[1],cr_limits[1])
      upperList=c(lnA_limits[2],Ea_limits[2],n_limits[2],m_limits[2],
                  lnA_limits[2],Ea_limits[2],n_ePT_limits[2],m_limits[2],
                  lnA_limits[2],Ea_limits[2],n_ePT_limits[2],m_limits[2],
                  cd1max,cd1max,cr_limits[2],cr_limits[2])
      nInd=0
      if (isolate(input$lnA1f)==T) {
        lnA1=isolate(input$lnA1fixed)
        indF=1-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$Ea1f)==T) {
        Ea1=isolate(input$Ea1fixed)
        indF=2-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$nf1)==T) {
        n1=isolate(input$n1fixed)
        indF=3-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$mf1)==T) {
        m1=isolate(input$m1fixed)
        indF=4-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$lnA2f)==T) {
        lnA2=isolate(input$lnA2fixed)
        indF=5-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$Ea2f)==T) {
        Ea2=isolate(input$Ea2fixed)
        indF=6-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$nf2)==T) {
        n2=isolate(input$n2fixed)
        indF=7-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$mf2)==T) {
        m2=isolate(input$m2fixed)
        indF=8-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$lnA3f)==T) {
        lnA3=isolate(input$lnA3fixed)
        indF=9-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$Ea3f)==T) {
        Ea3=isolate(input$Ea3fixed)
        indF=10-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$nf3)==T) {
        n3=isolate(input$n3fixed)
        indF=11-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$mf3)==T) {
        m3=isolate(input$m3fixed)
        indF=12-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$cd1f)==T) {
        cd1=isolate(input$cd1fixed)
        indF=13-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$datasetsCnt)==1) {
        cd2=0
        indF=14-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] }
      q2=ifelse(isolate(input$q2f)==T,isolate(input$q2),isolate(input$q))
      q3=ifelse(isolate(input$q3f)==T,isolate(input$q3),isolate(input$q))
      if (optimType==1) 
          nlrRes=try(nlsLM(nlY ~ funK9(lnA1,Ea1,n1,m1,lnA2,Ea2,n2,m2,lnA3,Ea3,n3,m3,cd1,cd2,cr21,cr22,T1t,tt,at,flag,nr,Dtype,
                                       wghDA,method,q,q2,q3,amin,amax,Nruns,TempF,suppression,ti_trigg), 
                           start=startList,lower=lowerList, upper=upperList,control= nls.lm.control(maxiter=1000,maxfev=10000)))
      else if (optimType==2) 
          nlrRes=try(nls(nlY ~ funK9(lnA1,Ea1,n1,m1,lnA2,Ea2,n2,m2,lnA3,Ea3,n3,m3,cd1,cd2,cr21,cr22,T1t,tt,at,flag,nr,Dtype,
                                     wghDA,method,q,q2,q3,amin,amax,Nruns,TempF,suppression,ti_trigg), 
                         start=startList,algorithm = "port",lower=lowerList, upper=upperList,
                         control= nls.control(maxiter=1000, warnOnly=T))) 
    }    

    if (isolate(input$taskType)==10) {
      cd2val=ifelse(isolate(input$datasetsCnt)==2,isolate(input$cd2),0)
      startList=isolate(list(lnA1=input$lnA1,Ea1=input$Ea1,n1=input$n1,m1=input$m1,
                             lnA2=input$lnA2,Ea2=input$Ea2, n2=input$n2,m2=input$m2,
                             lnA3=input$lnA3,Ea3=input$Ea3, n3=input$n3,m3=input$m3,
                             cd1=input$cd1,cd2=cd2val))
      cd1min=isolate(input$cd1Lim[1])
      cd1max=isolate(input$cd1Lim[2])
      lowerList=c(lnA_limits[1],Ea_limits[1],n_limits[1],m_limits[1],
                  lnA_limits[1],Ea_limits[1],n_ePT_limits[1],m_limits[1],
                  lnA_limits[1],Ea_limits[1],n_ePT_limits[1],m_limits[1],
                  cd1min,cd1min)
      upperList=c(lnA_limits[2],Ea_limits[2],n_limits[2],m_limits[2],
                  lnA_limits[2],Ea_limits[2],n_ePT_limits[2],m_limits[2],
                  lnA_limits[2],Ea_limits[2],n_ePT_limits[2],m_limits[2],
                  cd1max,cd1max)
      nInd=0
      if (isolate(input$lnA1f)==T) {
        lnA1=isolate(input$lnA1fixed)
        indF=1-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$Ea1f)==T) {
        Ea1=isolate(input$Ea1fixed)
        indF=2-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$nf1)==T) {
        n1=isolate(input$n1fixed)
        indF=3-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$mf1)==T) {
        m1=isolate(input$m1fixed)
        indF=4-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$lnA2f)==T) {
        lnA2=isolate(input$lnA2fixed)
        indF=5-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$Ea2f)==T) {
        Ea2=isolate(input$Ea2fixed)
        indF=6-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$nf2)==T) {
        n2=isolate(input$n2fixed)
        indF=7-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$mf2)==T) {
        m2=isolate(input$m2fixed)
        indF=8-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$lnA3f)==T) {
        lnA3=isolate(input$lnA3fixed)
        indF=9-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$Ea3f)==T) {
        Ea3=isolate(input$Ea3fixed)
        indF=10-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$nf3)==T) {
        n3=isolate(input$n3fixed)
        indF=11-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$mf3)==T) {
        m3=isolate(input$m3fixed)
        indF=12-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$cd1f)==T) {
        cd1=isolate(input$cd1fixed)
        indF=13-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] 
        nInd=nInd+1 }
      if (isolate(input$datasetsCnt)==1) {
        cd2=0
        indF=14-nInd
        startList[[indF]]=NULL
        lowerList=lowerList[-indF]
        upperList=upperList[-indF] }
      q2=ifelse(isolate(input$q2f)==T,isolate(input$q2),isolate(input$q))
      q3=ifelse(isolate(input$q3f)==T,isolate(input$q3),isolate(input$q))
      swA=input$swA
      if (optimType==1) 
        nlrRes=try(nlsLM(nlY ~ funK10(lnA1,Ea1,n1,m1,lnA2,Ea2,n2,m2,lnA3,Ea3,n3,m3,cd1,cd2,swA,T1t,tt,at,flag,nr,Dtype,
                                     wghDA,method,q,q2,q3,amin,amax,Nruns,TempF,suppression), 
                         start=startList,lower=lowerList, upper=upperList,control= nls.lm.control(maxiter=1000,maxfev=10000)))
      else if (optimType==2) 
        nlrRes=try(nls(nlY ~ funK10(lnA1,Ea1,n1,m1,lnA2,Ea2,n2,m2,lnA3,Ea3,n3,m3,cd1,cd2,swA,T1t,tt,at,flag,nr,Dtype,
                                   wghDA,method,q,q2,q3,amin,amax,Nruns,TempF,suppression), 
                       start=startList,algorithm = "port",lower=lowerList, upper=upperList,
                       control= nls.control(maxiter=1000, warnOnly=T))) 
    }  
    
    end = Sys.time()
    #cat(file=stderr(), end-beginning, "\n")
    #Dat=cbind(Tt,at,predict(nlrRes))
    #write.table(Dat, file = "logfilename2.txt", sep=" ", row.names=FALSE, quote = FALSE) 
    #Dat=cbind(Tt,at,funK1(45.9386,169.5096,1,T1t,tt,at))#predict(nlrRes1),predict(nlrRes2),predict(nlrRes3),predict(nlrRes4))
    #write.table(Dat, file = "logfilename.txt", sep=" ", row.names=FALSE, quote = FALSE)
    if (inherits(nlrRes, "try-error"))   predNLR="ERROR" 
    else predNLR = predict(nlrRes)
    shinyjs::enable("go")
    list(nlrRes,difftime(end,beginning,units='secs'),predNLR)
  })
  
  output$reggrSummary <- renderPrint({
    if (input$go==0) cat("Regression summary will shown be here")
    else {
      res=reggr()
      if (res[[3]][1] == "ERROR") cat(res[[1]][1])
      else summary(res[[1]]) }
  })

  output$reggrSummaryTime <- renderPrint({
    if (input$go==0) cat("The time of NLR optimization will shown be here")
    else {
      res=reggr()
      if (res[[3]][1] == "ERROR") return(NULL)
      else {
        Duration=reggr()[[2]]
        cat("Running time of NLR: ", round(Duration, digits=1), "seconds") } }
  })
    
  output$reggrSummary2 <- renderPrint({
    if (input$go==0) cat("The statistical quality of model fit will shown be here")
    else {
      res=reggr()
      if (res[[3]][1] == "ERROR") return(NULL)
      nlrRes=res[[1]]
      resOut=array(dim=c(1,0))
      if(match(1,isolate(input$statM),nomatch = 0)) {
        RSS=sum(residuals(nlrRes)^2)
        resOut=cbind(resOut,round(as.numeric(RSS)*10)/10)
        colnames(resOut)[ncol(resOut)]="RSS"
      }
      if(match(3,isolate(input$statM),nomatch = 0)) {
        AICval=AIC(nlrRes)
        resOut=cbind(resOut,round(as.numeric(AICval)*10)/10)
        colnames(resOut)[ncol(resOut)]="AIC"
      }
      if(match(4,isolate(input$statM),nomatch = 0)) {
        BICval=BIC(nlrRes)
        resOut=cbind(resOut,round(as.numeric(BICval)*10)/10)
        colnames(resOut)[ncol(resOut)]="BIC"
      }
      noquote(resOut) }
  })
  
  
})
