#THINKS thermokinetic software, version #01.07.20 Website: thinks.chemphys.ru. To cite: https://doi.org/10.3390/molecules24122298
#This routine performs the model-fitting linear regression analysis.
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
    
  treatData <- reactive({ #from input kinetic data creates single rows of conversion, temperature etc. and cuts data at a<a(min) and a>1-a(min)
    inFile <- input$file1
    data=read.table(inFile$datapath,header = TRUE, check.names = FALSE)
    Nruns=ncol(data)/4
    Npts=nrow(data)
    ti=data[[2]]
    Te=data[[1]]+273.15
    al=data[[3]]
    da=data[[4]]
    for (i in 1:Nruns) {        
      ti=c(ti,data[[4*(i-1)+2]])
      Te=c(Te,data[[4*(i-1)+1]]+273.15)
      al=c(al,data[[4*(i-1)+3]])
      da=c(da,data[[4*(i-1)+4]]) }
    AMin=ifelse(((input$aMin<=0.3)&&(input$aMin>=0.01)),input$aMin,0.05)
    AMax=ifelse(((input$aMax<=0.99)&&(input$aMax>=0.6)),input$aMax,0.95)
    RightInd=which((al>=AMin)&(al<=AMax))
    al=al[RightInd]
    da=da[RightInd]
    Te=Te[RightInd]
    lowerInd=numeric(Nruns)
    higherInd=numeric(Nruns)
    for (i in 1:Nruns) {
      insideInd=which((RightInd>=((i-1)*Npts+1))&(RightInd<=i*Npts))
      lowerInd[i]=min(insideInd)
      higherInd[i]=max(insideInd) }
    Te1=1/Te
    ldat=ifelse(da<=0,0,log(da))
    l1at=ifelse(al>=1,0,log(1-al))
    lat=ifelse(al<=0,0,log(al))
    list(ti,al,da,Te,Te1,Te1/8.314,-1000*Te1/8.314,ldat,l1at,lat,lowerInd,higherInd,Nruns) 
  })

  reggr=reactive({ #linear regression (simple case) or nonlinear regression (with constraints defined)
    res=treatData()
    ldat=res[[8]]
    T1t=res[[7]]
    l1at=res[[9]]
    lat=res[[10]]
    input$go
    if(isolate(input$DL)==TRUE) {
      startList=isolate(list(lncA=(input$lA[1]+input$lA[2])/2, Ea=(input$Ea[1]+input$Ea[2])/2,
                             n=(input$n1[1]+input$n1[2])/2, m=(input$m1[1]+input$m1[2])/2))
      lowerList=isolate(list(input$lA[1],input$Ea[1],input$n1[1],input$m1[1]))
      upperList=isolate(list(input$lA[2],input$Ea[2],input$n1[2],input$m1[2]))
      lmRes=nls(ldat ~ lncA+Ea*T1t+n*l1at+m*lat, algorithm = "port",start=startList,lower=lowerList, upper=upperList) 
    }
    else { 
      lmRes=lm(ldat ~ T1t+l1at+lat) #-------change names?
      names(lmRes$coefficients) = c('lncA','Ea','n', 'm') }
    lmRes
  })
  
  output$distPlot <- renderPlot( height= function(){ ifelse(input$Compact==T, 350, input$BaseHeight)   },
    { #plot of ln(da/dt)-ln((1-a)^na^m) against 1/Temperature
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    res=treatData()
    NrunsD=res[[13]]
    if(!input$go) return(NULL)
    lmRes=reggr()
    StuCoef=qt(1-.05/2, df=df.residual(lmRes))
    stdErrors=coef(summary(lmRes))[,2]*StuCoef
    xLin=seq(min(res[[5]]),max(res[[5]]),length.out=3)
    yLin=coef(lmRes)[1]-1000*coef(lmRes)[2]*xLin/8.314
    par (fig=c(0,0.8,0,1.0),mar=c(5.1,5.1,4.1,0)) #par(mar = c(5,7,4,2) + 0.1)
    plot(xLin,yLin,type="l",lwd=2, xlab=NA, ylab=NA, las=1, cex.axis=input$LegendSize) 
    title(xlab=expression('Temperature'^-1*', K'^-1), main="Combined kinetic analysis plot",
          ylab=expression('ln(d'*alpha*'/d'*italic('t')*') - ln[(1-'*alpha*')'^italic('n')*alpha^italic('m')*']'),
          font.main=input$TitleFont,cex.lab=input$AxisSize,cex.main=input$TitleSize) 
    yF=res[[8]]-coef(lmRes)[3]*res[[9]]-coef(lmRes)[4]*res[[10]]
    xF=res[[5]] 
    #if one needs output data:
    #write.table(cbind(xLin,yLin), file = "line.txt", sep=" ", row.names=FALSE, quote = FALSE)
    #write.table(cbind(xF,yF), file = "points.txt", sep=" ", row.names=FALSE, quote = FALSE)
    col = rainbow(NrunsD)
    lowerInd=res[[11]]
    higherInd=res[[12]]
    legA=character(NrunsD)
    legA[1]="linear fit"
    lwdA=numeric(NrunsD)
    lwdA[1]=2
    colA=character(NrunsD)
    colA[1]="black"
    pchA=numeric(NrunsD)
    pchA[1]=NA
    for (i in 1:NrunsD) {
        lI=lowerInd[i]
        hI=higherInd[i]
        points(xF[lI:hI],yF[lI:hI], col=col[i], pch=19) 
        legA[i+1]=paste0("run ",i)
        lwdA[i+1]=NA
        colA[i+1]=col[i]
        pchA[i+1]=19 }    
    usr=par("usr")
    legend(usr[2], usr[4], bty='n', 
           legend=legA, lwd=lwdA, pch=pchA, col=colA, cex=input$LegendSize,  xpd=NA)
    
    xCoord=(usr[2]-usr[1])*29/30+usr[1]
    yCoordmin=0.07
    yCoordStep=0.07   
    mess10=paste(round(coef(lmRes)[2],digits=1),"±",round(stdErrors[2],digits=1)," kJ/mol")
    mess1=bquote(italic('E')[a]*plain(" = ")*.(mess10))
    text(xCoord,(usr[3]-usr[4])*yCoordmin+usr[4], mess1, pos=2, cex=input$LegendSize)
    mess20=paste(round(coef(lmRes)[1],digits=1),"±",round(stdErrors[1],digits=1)," 1/s")
    mess2=bquote(plain("ln(")*italic('cA')*plain(") = ")*.(mess20))    
    text(xCoord,(usr[3]-usr[4])*(yCoordmin+yCoordStep)+usr[4], mess2, pos=2, cex=input$LegendSize)
    mess3=bquote(italic('n')*plain(" = ")*.(round(coef(lmRes)[3],digits=2))*plain(" ± ")*.(round(stdErrors[3],digits=2)))    
    text(xCoord,(usr[3]-usr[4])*(yCoordmin+yCoordStep*2)+usr[4], mess3, pos=2, cex=input$LegendSize)
    mess4=bquote(italic('m')*plain(" = ")*.(round(coef(lmRes)[4],digits=2))*plain(" ± ")*.(round(stdErrors[3],digits=2)))    
    text(xCoord,(usr[3]-usr[4])*(yCoordmin+yCoordStep*3)+usr[4], mess4, pos=2, cex=input$LegendSize)
  })

  output$distPlot2 <- renderPlot( height= function(){ ifelse(input$Compact==T, 300, input$BaseHeight)   },
    { #plot with reduced reaction model types f(a)/f(0.5)
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    res=treatData()
    if(!input$go) return(NULL)
    lmRes=reggr()
    #par(mar = c(5,7,4,2) + 0.1)
    #par(oma=c(0, 0, 0, 5))
    par (fig=c(0,0.8,0,1.0),mar=c(5.1,5.1,2.5,0)) 
    AMin=ifelse(((isolate(input$aMin)<=0.3)&&(isolate(input$aMin)>=0.01)),isolate(input$aMin),0.05)
    AMax=ifelse(((isolate(input$aMax)<=0.99)&&(isolate(input$aMax)>=0.6)),isolate(input$aMax),0.95)
    alp=seq(AMin,AMax,length.out=20)
    ff=((1-alp)^coef(lmRes)[3]*alp^coef(lmRes)[4])/(0.5^coef(lmRes)[3]*0.5^coef(lmRes)[4])
    plot(alp,ff,type="p",cex=input$LegendSize, cex.lab=input$AxisSize, pch=19, xlim=c(0,1), ylim=c(0,2.5), 
         xlab=expression('Conversion '*alpha), main="Master plot", font.main=input$TitleFont, cex.main=input$TitleSize,
         ylab=expression(italic('f')*'('*alpha*') / '*italic('f')*'(0.5)'))
    Nmodels=length(input$reMod)
    if (Nmodels>5) col=rainbow(Nmodels)
    else col=c('red','blue','magenta','green','grey')
    legA=character(Nmodels)
    legA[1]="Model"
    lwdA=numeric(Nmodels)
    lwdA[1]=NA
    colA=character(Nmodels)
    colA[1]="black"
    pchA=numeric(Nmodels)
    pchA[1]=19
    alp2=seq(0.05,0.95,length.out=50)
    for (i in 1:Nmodels) {
      ffRM=reactionModel(input$reMod[i],alp2)/reactionModel(input$reMod[i],0.5)
      lines(alp2,ffRM, col=col[i], lwd=1) 
      legA[i+1]=input$reMod[i]
      lwdA[i+1]=1
      colA[i+1]=col[i]
      pchA[i+1]=NA }      
    legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA, legend=legA, lwd=lwdA, pch=pchA, col=colA,
           cex=input$LegendSize)
    })
  
  output$text <- renderPrint({ #returns full regression output if selected
    if ((input$go==0)||(is.null(input$file1)))  cat("Here will be the regression summary")
    else {
      lmRes=reggr()
      summary(lmRes)
    }
  })

  reactionModel=function(Rtype,alpha) {  #idealized reaction model types. Five group of models from:
  #S. Vyazovkin, A. K. Burnham, J. M. Criado, et al, ICTAC Kinetics Committee recommendations for performing kinetic computations on thermal analysis data, Thermochimica Acta, Vol. 520 (1–2), 2011, pp. 1–19.
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
    if (Rtype=="A2") res=2*(1-alpha)*(-log(1-alpha))^0.5
    if (Rtype=="A3") res=3*(1-alpha)*(-log(1-alpha))^(2/3)
    if (Rtype=="A4") res=4*(1-alpha)*(-log(1-alpha))^(3/4)
    if (Rtype=="P2") res=2*alpha^0.5
    if (Rtype=="P3") res=3*alpha^(2/3)
    if (Rtype=="P4") res=4*alpha^(3/4)
    if (Rtype=="P23") res=2/3*alpha^(-0.5)
    if (Rtype=="D1") res=0.5/alpha
    if (Rtype=="D2") res=-1/log(1-alpha)
    if (Rtype=="D3") res=1.5*(1-alpha)^(2/3)/(1-(1-alpha)^(1/3))
    if (Rtype=="D4") res=1.5/((1-alpha)^(-1/3)-1)
    if (Rtype=="L2") res=2*(alpha^0.5-alpha)
    if (Rtype=="B1") res=(1-alpha)*alpha
    res }
  
})
