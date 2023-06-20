#THINKS thermokinetic software, version #01.07.20 Website: thinks.chemphys.ru. To cite: 10.3390/molecules24122298
#This routine performs the isoconversional analysis.
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
  
  treatData <- reactive({ #creates the input data for isoconversional methods by alpha in a(min)..a(max) with defined step
    inFile <- input$file1
    data=read.table(inFile$datapath,header = TRUE, check.names = FALSE)
    Nruns=ncol(data)/4
    Npts=nrow(data)
    alpha=suppressWarnings(as.numeric(as.character(data[[3]])))
    alpha[alpha<1E-6]=0   #-----------------------------------------------------add this check to Import routine!?
    Temp=suppressWarnings(as.numeric(as.character(data[[1]]+273.15)))
    tFunc=suppressWarnings(splinefun(alpha,data[[2]],method = "fmm"))
    TFunc=suppressWarnings(splinefun(alpha,Temp,method = "fmm"))
    daFunc=suppressWarnings(splinefun(alpha,data[[4]],method = "fmm"))
    
    alExtr=double(Nruns)    # calculation by Kissinger method (to draw its results on the Ea(alpha) plot)
    TExtr=double(Nruns)
    HRs=double(Nruns)
    HRs[1]=line(data[[2]],data[[1]]+273.15)$coefficients[2]
    if (HRs[1]>(0.05/60)) {  # finds maximum of da/dt for non-isothermal runs
        alExtr[1]=optimize(daFunc,c(0,1), maximum = TRUE)$maximum  
        TExtr[1]=TFunc(alExtr[1]) }
    else {
      alExtr[1]=0
      TExtr[1]=0 }
    AStep=ifelse(input$aStep>=0.005,input$aStep,0.01)  #-------thinning-out of data
    alNew=seq(input$aMin,input$aMax,by=AStep)
    Nal=length(alNew)
    ti=matrix(nrow=Nal,ncol=Nruns)
    Te=matrix(nrow=Nal,ncol=Nruns)
    da=matrix(nrow=Nal,ncol=Nruns)
    lda=matrix(nrow=Nal,ncol=Nruns)
    ti[,1]=tFunc(alNew)
    Te[,1]=TFunc(alNew)
    da[,1]=daFunc(alNew)
    for (k in 1:Nal) lda[k,1]=ifelse(da[k,1]<=0,0,log(da[k,1]))
    for (i in 2:Nruns) {  #-----------------------------------------------all the same for next runs
      alpha=suppressWarnings(as.numeric(as.character(data[[4*(i-1)+3]])))
      Temp=suppressWarnings(as.numeric(as.character(data[[4*(i-1)+1]]+273.15)))
      alpha[alpha<1E-6]=0
      tFunc=suppressWarnings(splinefun(alpha,data[[4*(i-1)+2]],method = "fmm"))
      TFunc=suppressWarnings(splinefun(alpha,Temp,method = "fmm"))
      daFunc=suppressWarnings(splinefun(alpha,data[[4*(i-1)+4]],method = "fmm"))
      
      HRs[i]=line(data[[4*(i-1)+2]],data[[4*(i-1)+1]]+273.15)$coefficients[2]
      if (HRs[i]>(0.05/60)) {   # finds maximum of da/dt for non-isothermal runs
          alExtr[i]=optimize(daFunc,c(0,1), maximum = TRUE)$maximum
          TExtr[i]=TFunc(alExtr[i]) }
      else {
          alExtr[i]=0
          TExtr[i]=0 }
      ti[,i]=tFunc(alNew)
      Te[,i]=TFunc(alNew)
      da[,i]=daFunc(alNew)
      for (k in 1:Nal) lda[k,i]=ifelse(da[k,i]<=0,0,log(da[k,i])) }
    HRs=HRs[HRs>(0.05/60)]
    TExtr=TExtr[TExtr!=0]
    alExtr=alExtr[alExtr!=0]
    Te1=1/Te
    list(ti,alNew,da,Te,Te1,-Te1/8.314,lda,Nruns,Nal,alExtr,
         TExtr,HRs) 
  })

  Kissinger=reactive({
    resK=treatData()
    yKiss=log(resK[[12]]/(resK[[11]]^2))     # Kissinger method block -add preexponent?
    xKiss=1/resK[[11]]
    lmRes=lm(yKiss~xKiss)
    lmRes_Summary=summary(lmRes)
    std_Err=qt(0.95,df=lmRes$df.residual)*lmRes_Summary$coefficients[,2]
    lgA=log10(-exp(lmRes$coefficients[1])*lmRes$coefficients[2])
    dlgA=((std_Err[2]/exp(lmRes$coefficients[1])/lmRes$coefficients[2])^2+std_Err[1]^2)^0.5/log(10)
    list(lmRes, -lmRes$coefficients[2]*8.314/1000, std_Err[2]*8.314/1000, lgA, dlgA,
         lmRes_Summary$r.squared)
  })
  
  output$textKiss <- renderPrint({ #returns full regression output for Kissinger
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    if (input$showKiss==F) return(NULL)
    else {
      reggr=Kissinger()
      Ea=Kissinger()[[2]]
      dEa=Kissinger()[[3]]
      lgA=Kissinger()[[4]]
      dlgA=Kissinger()[[5]]
      R2=Kissinger()[[6]]
      cat(c("","/Summary of the Kissinger analysis/","\n",
            "Ea = ",round(Ea,digits=2),"±",round(dEa,digits=2),"kJ/mol", "\n",
            "lg (A, 1/s) = ",round(lgA,digits=2),"±",round(dlgA,digits=2),"\n",
            "R2 = ", round(R2,digits=3)
      ),sep=" ") }
  })  
  
  Friedman=reactive({ #Friedman method
    res=treatData()
    Nruns=res[[8]]
    Nal=res[[9]]
    al=res[[2]]
    ldat=res[[7]]
    T1t=res[[6]]
    Ea=double(Nal)
    dEa=double(Nal)
    lnf=double(Nal)
    dlnf=double(Nal)
    R2=double(Nal)
    input$go
    for (i in 1:Nal) { 
      lmRes=lm(ldat[i,] ~ T1t[i,])
      Ea[i]=coef(lmRes)[2]
      lnf[i]=coef(lmRes)[1]
      StuCoef=qt(1-.05/2, df=df.residual(lmRes))
      stdErrors=coef(summary(lmRes))[,2]#*StuCoef
      dEa[i]=stdErrors[2]
      dlnf[i]=stdErrors[1]
      R2[i]=summary(lmRes)$r.squared
    }
    #write.table(cbind(al,Ea/1000,dEa/1000,lnf,dlnf,R2), file = "Friedman_results.txt", sep=" ", row.names=FALSE, quote = FALSE)
    list(al,Ea/1000,dEa/1000,lnf,dlnf,R2)
  }) 
  
  Vyazovkin=reactive({ #advanced Vyazovkin method [Journal of Computational Chemistry, Vol.22, No.2, 178-183 (2001)]
    res=treatData() #-------------list(ti,alNew,da,Te,Te1,-Te1/8.314,lda,Nruns,Nal,alExtr, TExtr,HRs) 
    Nruns=res[[8]]
    Nal=res[[9]]
    al_old=res[[2]]
    d_al=al_old[2]-al_old[1]
    ldat=res[[7]]
    T1t=res[[6]]
    ti=res[[1]]
    al=double(Nal-1)
    Ea=double(Nal-1)
    dEa=double(Nal-1)
    lnf=double(Nal-1)
    dlnf=double(Nal-1)
    R2=double(Nal-1)
    input$go
    for (i in 2:Nal) { 
      E1=10000 #-----minimal possible activation energy = 10 kJ/mol
      E2=600000 #--maximal possible activation energy = 600 kJ/mol
      E3=(E1+E2)/2
      F_EE <- function(EE) {
        ress=0
        for (j in 1:Nruns) {
          for (k in 1:Nruns) {  #-check this use the trapezoid rule
            if(k != j) ress=ress+ (ti[i,j]-ti[i-1,j])*(exp(EE*T1t[i,j])+exp(EE*T1t[i-1,j]))  /
                                    ((ti[i,k]-ti[i-1,k])*(exp(EE*T1t[i,k])+exp(EE*T1t[i-1,k])))   # -1)^2 
          }
        }
        return(ress/((Nruns-1)*Nruns))
      }
#      vY=matrix(c(F_EE(E1),F_EE(E2),F_EE(E3)))
#      vX=matrix(c(E1^2,E2^2,E3^2, E1,E2,E3, 1,1,1), nrow=3, ncol=3)
#      v=vX^-1 %*% vY
#      E0=-v[2]/2/v[1]
#      EaP=E0
#      delta_Ea=1000 #---the accuracy of the output activation energy = 1 kJ/mol
#      delta=delta_Ea+10
      
#      while (delta>delta_Ea) {
#        if (F_EE(E1)>F_EE(E1+delta)) {
#          EaP=E1
#          while (F_EE(EaP+delta)<F_EE(EaP)) EaP=EaP+delta
#          E1=EaP-delta
#          E2=EaP+delta
#          delta=(E2-E1)/5
#        }
#        else if (F_EE(0.95*E1)<F_EE(E1)) E1=0.95*E1
#      }
      EaP=optimize(F_EE, interval=c(E1,E2), maximum=FALSE)$minimum
      al[i-1]=al_old[1]+(i-1)*d_al
      Ea[i-1]=EaP
      j=floor(Nruns/2)-1 
      A=2*d_al/(ti[i,j]-ti[i-1,j])/(exp(EaP*T1t[i,j])+exp(EaP*T1t[i-1,j]))
      lnf[i-1]=ifelse(A<=0,0,log(A))
      dEa[i-1]=0
      dlnf[i-1]=0
      R2[i-1]=F_EE(EaP)*((Nruns-1)*Nruns) #---it is not R^2 but Ω (optimization indicator [DOI: 10.1016/j.tca.2016.05.004]) 
    }
    list(al,Ea/1000,dEa/1000,lnf,dlnf,R2)
  })   
  
  output$Plot1 <- renderPlot(height= function(){ ifelse(input$Compact==T, ifelse(input$showKiss==T,300, 350), input$BaseHeight)   },
                             { #--------------------main plot with the isoconversional data
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    if (input$Method==1) res=Friedman()
    else res=Vyazovkin()
    par (fig=c(0,0.8,0,1.0),
         mar=c(5.1+input$marginBottom,5.1+input$marginLeft,3.1,0),
         mgp=c(3+input$marginTitle,1,0)) #par(oma=c(0,1,0,0))
    #par (fig=c(0,0.8,0,1.0),mar=c(5.1,5.1,3.1,0)) #par (fig=c(0,0.8,0,1.0),mar=c(6.1,4.9,4.1,0),mgp=c(2.6,0.7,0),tcl=-0.5)
    EaStep=max((max(res[[2]])-min(res[[2]]))/10,50)
    EaStep2=min((max(res[[2]])-min(res[[2]]))/10,20)
    plot(res[[1]],res[[2]], type="n", xlim=c(0, 1), ylim=c(min(res[[2]])-EaStep, max(res[[2]])+EaStep),
           xlab=NA, ylab=NA,pch=21, cex=2.5, bg="white", col="red",las=1, cex.axis=input$LegendSize)
    title(main="Isoconversional plot of activation energy",
          xlab=expression('Conversion '~alpha),ylab=expression(italic('E')[a]*', kJ mol'^-1),
          cex.lab=input$AxisSize,cex.main=input$LegendSize)
    arrows(res[[1]],res[[2]]-res[[3]],res[[1]],res[[2]]+res[[3]],length=0,angle=90,code=3,col="red",lwd=0.5)
   # if(input$KinFT==1) EaFunc=splinefun(res[[1]],res[[2]],method = "fmm")
  #  if(input$KinFT==2) EaFunc=splinefun(res[[1]],res[[2]],method = "natural")
  #  if(input$KinFT==3) EaFunc=approxfun(res[[1]],res[[2]],method = "linear",rule=2)
    EaFunc=suppressWarnings(splinefun(res[[1]],res[[2]],method = "natural"))
    alpSeq=seq(input$aMin,input$aMax,by=0.02)
    lines(alpSeq,EaFunc(alpSeq), col = "red",lwd=2) 
    #lines(lowess(res[[1]],res[[2]],f=.1), col = "red",lwd=1.5) 
    points(res[[1]],res[[2]], pch=21, cex=2.5, bg="white", col="red")
    if (input$showBoth==T) {
      if (input$Method==1) res2=Vyazovkin()
      else res2=Friedman()
      lines(res2[[1]],res2[[2]], col = "blue",lwd=2) }
    if (input$showKiss==T) {
      fitKiss=Kissinger()[[1]]
      EaKiss=-coef(fitKiss)[2]*8.314/1000
      #StuCoef=qt(1-.05/2, df=df.residual(fitKiss))
      #stdErrors=coef(summary(fitKiss))[,2]*StuCoef
      dEaKiss=coef(summary(fitKiss))[2,2]*8.314/1000                 #check the error evaluation!
      resK=treatData()
      radiusX=(max(resK[[10]])-min(resK[[10]]))/2
      radiusY=dEaKiss
      theta=seq(0, 2 * pi, length = 200)
      lines(x = mean(resK[[10]])+radiusX*cos(theta), y = EaKiss+radiusY*sin(theta))
    }
    par (fig=c(0.85,1.0,0,1.0),
         mar=c(5.1+input$marginBottom,0,3.1,1.0),new=TRUE) 
    #par (fig=c(0.85,1.0,0,1.0), mar=c(5.1,0,3.1,1.5),new=TRUE) #par (fig=c(0.85,1.0,0,1.0), mar=c(6.1,0,4.1,1.5),new=TRUE)
    hs=hist(res[[2]], plot=FALSE, breaks=seq(min(res[[2]])-EaStep, max(res[[2]])+EaStep, EaStep2)) #breaks="FD"
    barplot(hs$counts, space=0, horiz=T, col = "red",border="white", axes=F)
  })
  
  output$Plot2 <- renderPlot(height= function(){ ifelse(input$Compact==T, ifelse(input$showKiss==T,230, 270), input$BaseHeight)   },
  { #------------------------------------------------Other isoconversional data
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    if (input$Method==1) res=Friedman()
    else res=Vyazovkin()
    par (fig=c(0,0.8,0,1.0),
         mar=c(5.1+input$marginBottom,5.1+input$marginLeft,3.1,0),
         mgp=c(3+input$marginTitle,1,0))
    #par (fig=c(0,0.8,0,1.0),mar=c(4.1,5.1,3.1,0)) #par (fig=c(0,0.8,0,1.0),mar=c(6.1,4.9,4.1,0),mgp=c(2.6,0.7,0),tcl=-0.5)
    lnAfStep=max((max(res[[4]])-min(res[[4]]))/10,50)
    lnAfStep2=min((max(res[[4]])-min(res[[4]]))/10,10)
    plot(res[[1]],res[[4]], type="n", xlim=c(0, 1), ylim=c(min(res[[4]])-lnAfStep, max(res[[4]])+lnAfStep),
         xlab=NA, ylab=NA,pch=21, cex.axis=input$LegendSize, bg="white", col="blue",las=1)
    title(main="Isoconversional results",
      xlab=expression('Conversion'~alpha),ylab=expression('ln('*italic('A*f')*')'),
          cex.lab=input$AxisSize,cex.main=input$LegendSize)
    arrows(res[[1]],res[[4]]-res[[5]],res[[1]],res[[4]]+res[[5]],length=0,angle=90,code=3,col="blue",lwd=0.5)
    axis(side = 2, col="blue", labels=F)
  #  if(input$KinFT==1) lnAfFunc=splinefun(res[[1]],res[[4]],method = "fmm")
  #  if(input$KinFT==2) lnAfFunc=splinefun(res[[1]],res[[4]],method = "natural")
  #  if(input$KinFT==3) lnAfFunc=approxfun(res[[1]],res[[4]],method = "linear",rule=2)
    lnAfFunc=suppressWarnings(splinefun(res[[1]],res[[4]],method = "natural"))
    alpSeq=seq(input$aMin,input$aMax,by=0.02)
    lines(alpSeq,lnAfFunc(alpSeq), col = "red",lwd=2) 
    #lines(lowess(res[[1]],res[[4]],f=.1), col = "blue",lwd=1.5) 
    points(res[[1]],res[[4]], pch=21, cex=2.5, bg="white", col="red")
    if (input$showBoth==T) {
      if (input$Method==1) res2=Vyazovkin()
      else res2=Friedman()
      lines(res2[[1]],res2[[4]], col = "blue",lwd=2) }
    if (input$showR2==T) {
        par(new = T)
        plot(res[[1]],res[[6]], type="l", col="gray36", lwd=1.5, axes=F, xlab=NA, ylab=NA, cex.axis=input$LegendSize)
        mtext(ifelse(input$Method==1,expression(italic('R')^2),"Ω"),
              side=4,col="gray36",line=4,cex=input$AxisSize)
        axis(side = 4, col="gray36",col.axis="gray36",las=1, cex.axis=input$LegendSize) }
  })
  
  output$Plot4 <- renderPlot(height= function(){ ifelse(input$Compact==T, 250, input$BaseHeight)   },
                             { #-------Friedman plot
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    res=treatData()
    al=res[[2]]
    ldat=res[[7]]
    T1t=res[[5]]*1000
    resF=Friedman()
    Ea=resF[[2]]
    lnAf=resF[[4]]
    #windowsFonts(A = windowsFont("Times New Roman"))
    #windowsFonts(B = windowsFont("Arial"))
    daX=ifelse((input$daX>=1)&&(input$daX<=10),input$daX,3) 
  #  AStep=ifelse(input$aStep>=0.005,input$aStep,0.01)
  #  al=seq(input$aMin, input$aMax, by=AStep)
   # alNewX=seq(input$aMin+AStep*daX, input$aMax-AStep*daX, by=AStep*daX)
    alNewX=seq(min(al),max(al), by = daX*(al[2]-al[1]))
    #alNewX=al[indAl]
    indAl=which(al%in%alNewX)
    T1t=T1t[indAl,]
    ldat=ldat[indAl,]
    Ea=Ea[indAl]
    lnAf=lnAf[indAl]
    NalX=length(indAl)
    par (#fig=c(0,0.8,0,1.0),
         mar=c(5.1+input$marginBottom,6.1+input$marginLeft,3.1,6.7),
         mgp=c(3+input$marginTitle,1,0))
    #par (fig=c(0,0.8,0,1.0),mar=c(5.1,4.9,3.1,4.7), las=1) #par (mar=c(6.1,4.9,4.1,0),mgp=c(2.6,0.7,0),tcl=-0.5)
    plot(T1t,ldat, type="n", xlim = c((min(T1t)-(max(T1t)-min(T1t))/10), (max(T1t)+(max(T1t)-min(T1t))/10)),
         ylim=c((min(ldat)-(max(ldat)-min(ldat))/10), (max(ldat)+(max(ldat)-min(ldat))/10)), xlab=NA, ylab=NA,
         cex.axis=input$LegendSize, las=1)
    title(main="Friedman plot",
      xlab=expression(italic('T')^-1*', (kK)'^-1),ylab=expression('ln[d'~alpha*'/d'*italic('t')*', s'^-1*']'),
                    cex.lab=input$AxisSize,cex.main=input$LegendSize)
    if (NalX>5) col=rainbow(NalX)
    else col=c('red','blue','magenta','green','grey')
    legA=character(NalX)
    for (k in 1:NalX) {
      points(T1t[k,],ldat[k,], pch=21, cex=2.5, bg="white", col=col[k])  
      T1tfri=seq((min(T1t[k,])-(max(T1t[k,])-min(T1t[k,]))/10), (max(T1t[k,])+(max(T1t[k,])-min(T1t[k,]))/10), length.out=3)
      lines(T1tfri, (lnAf[k]-Ea[k]*T1tfri/8.314), col = col[k],lwd=1.5) 
      legA[k]=alNewX[k] #paste0("α = ",alNewX[k])
      }
    legend("topright", legend=legA, pch=21, col=col, bty = "n", cex=input$LegendSize, xpd=T, inset=c(input$LegendInset,0))
    
#    if (input$FriD==2) {  #3D-plot?!
#        x1=alNewX
#        x2=T1t
#        y=seq(0, 1, by=0.1)
#        x=seq(floor(min(x2)*20)/20-0.05, ceiling(max(x2)*20)/20+0.05, by=0.05)
#        z=outer(x,y,function(x,y){x*0+y*0+min(ldat)})
#        par (mar=c(3, 3, 2, 1) + 0.1, mgp=c(3, 1, 0),tcl=-0.5)
#        p=persp(x,y,z, zlim=c(min(ldat),max(ldat)), xlim=range(x), ylim=range(y),
#            xlab=expression(italic('T')^-1*' / (kK)'^-1), ylab=expression('Conversion'~alpha),
#            zlab=expression('ln[d'~alpha*'/d'*italic('t')*' / s'^-1*']'),
#            family="A",font.main=1,cex.lab=2,cex.main=2,
#            theta=45, phi=10, col="white", shade=0.01, ticktype="detailed", d=2)
#        yObs=alNewX[1]+ldat[1,]*0
#        obs=trans3d(T1t[1,], yObs,ldat[1,], p)
#        points(obs, col = "red", pch = 1, cex=2.5)  }
    
  })
  
  output$Plot5 <- renderPlot(height= function(){ ifelse(input$Compact==T, 250, input$BaseHeight)   },
    { #-------Kissinger plot
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    resK=treatData()
    yKiss=log(resK[[12]]/(resK[[11]]^2))    
    xKiss=1000/resK[[11]]
    par (#fig=c(0,0.8,0,1.0),
         mar=c(5.1+input$marginBottom,6.1+input$marginLeft,3.1,3.0),
         mgp=c(3.5+input$marginTitle,1,0))
    #par (fig=c(0,0.8,0,1.0),mar=c(5.1,4.9,3.1,4.7), las=1) 
    plot(xKiss,yKiss, type="p", xlab=NA, ylab=NA,pch=21, cex=2.5,
         cex.axis=input$LegendSize, las=1)
    title(main="Kissinger plot",
          xlab=expression(italic('T')[p]^-1*', (kK)'^-1),ylab=expression('ln('~beta*'/'*italic('T')[p]^2*')'),
                                     cex.lab=input$AxisSize,cex.main=input$LegendSize)
    abline(lm(yKiss ~ xKiss), col="red")
    })
  
  output$downloadData1 <- downloadHandler( #generates file with kinetic data 
    filename = function() {
      paste("IsoKin-", Sys.Date(), ".txt", sep="") },
    content = function(file) {
      if (input$Method==1) resF=Friedman()
      else resF=Vyazovkin()
      IsoKinData=cbind(resF[[1]],resF[[2]],resF[[3]],resF[[4]],resF[[5]],resF[[6]])
      last_col_name=ifelse(input$Method==1,"R^2","Omega")
      write.table(IsoKinData, file, sep=" ", eol = "\r\n", col.names=c("Conversion","Ea[kJ/mol]","dEa[kJ/mol]",
                                                                       "ln(A*f)","dln(A*f)",last_col_name),
                  row.names=FALSE, quote = FALSE) 
    })

  
})
