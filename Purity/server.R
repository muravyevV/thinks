#THINKS thermokinetic software, version #01.07.20 Website: thinks.chemphys.ru. To cite: 10.3390/molecules24122298
#This routine performs the purity analysis using the DSC data on melting of the respective substance.
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

shinyServer(function(input, output, session) {
  
  session$onSessionEnded(function() {
    observe({ 
      stopApp() })
  })
  
  observeEvent(input$go, {
    updateTabsetPanel(session, "MainTabset",
                      selected = "panel1")
  }, once=TRUE)
  
  simpson <- function(fun,a,b,n) { #calculation of the integral using Simpson equation
    n=ifelse(n<10,10,n)
    h=(b-a)/n
    x=seq(a, b, by=h)
    if (a==b) return (0)
    if (n==2) s=fun(x[1])+4*fun(x[2])+fun(x[3])  
    else s=fun(x[1])+fun(x[n+1])+2*sum(fun(x[seq(2,n,by=2)]))+4 *sum(fun(x[seq(3,n-1, by=2)])) 
    s*h/3 }

  readData <- reactive({#----------------------------------------------------Reading file defined in "file1" field------
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    InputData=read.table(inFile$datapath,sep=input$sep1, skip=input$BegSkipL)         #isolate and add "refresh"(action button)?
    removeLines=nrow(InputData)-input$EndSkipL
    InputData=InputData[1:removeLines,]
    x1=InputData[,as.numeric(input$tempCol)]-273.15*(as.numeric(input$tempUnits)-1)   #Temperature in [°C]
    if (as.numeric(input$dscUnits)>2)                                                 #DSC signal in W/g or μV/mg=mV/g
      y1=InputData[,as.numeric(input$dscCol)]/input$mass1
    else
      y1=InputData[,as.numeric(input$dscCol)]
    t1=InputData[,as.numeric(input$timeCol)]*as.numeric(input$timeUnits)              #Time in [s]
#    yPts=length(y1)                                                                  #Thinning out, if too many points?                                              
#    if (yPts>500) {
#        splD=splinefun(x1,y1,method = "fmm")
#        splT=splinefun(x1,t1,method = "fmm")
#        xN=seq(x1[1],x1[yPts],length.out = 500)                            
#        x1=xN
#        y1=splD(x1)
#        t1=splT(x1) }
    headF=head(as.data.frame(InputData), n=3L, addrownums = TRUE)
    list(x1,y1,t1,headF)  #Output: (T[°C],DSC[W/g or mV/g],t[s],first three lines from the file loaded)
    })                                   
  
  output$headFile <- renderPrint({#--------------------Returns first three lines from the file loaded as text output------
    if (is.null(input$file1)) cat("No data loaded") 
    else {
      headF=readData()[[4]]
      cat(c("First three lines of the file loaded:\n","[1]",t(headF[1,]),"\n","[2]",t(headF[2,]),"\n","[3]",t(headF[3,]),"\n"),sep=" ") }
    })
    
  output$BLpointsBeg <- renderUI({#--------------------Introduces slider with 10% shift from initial T in data--------------------
    if (is.null(input$file1)) return(NULL)
    res=readData()
    Tx=res[[1]]
    endTx=length(Tx) 
    Tmin=round(Tx[1]*10)/10
    Tmax=round(Tx[endTx]*10)/10
    Tval=round(((Tx[endTx]-Tx[1])*0.1+Tx[1])*10)/10
    sliderInput("sta", "Starting T :",min = Tmin, max = Tmax, value = Tval, step=0.01)
  })
  
  output$BLpointsEnd <- renderUI({#--------------------Introduces slider with 10% shift from final T in data--------------------
    if (is.null(input$file1)) return(NULL)
    res=readData()
    Tx=res[[1]]
    endTx=length(Tx) 
    Tmin=round(Tx[1]*10)/10
    Tmax=round(Tx[endTx]*10)/10
    Tval=round(((Tx[endTx]-Tx[1])*0.9+Tx[1])*10)/10
    sliderInput("end", "Final T :",min = Tmin, max = Tmax, value = Tval, step=0.01)
  })
  
  output$thinOut <- renderUI({#--------------------Slider for number of points reduction (thinning out)--------------------
    if (is.null(input$file1)) return(NULL)
    res=readData()
    Npts=length(res[[1]])
    IniValPts=floor(Npts/50)*50
    sliderInput("thinOutTo", "Reduce N pts to:",min = 50, max = Npts, value = IniValPts, step=10)
  })
  
  treatData <- reactive({#----------------------------Treats the imported data (this experiment data)---------------------
    if (input$go == 0) return()
    res=readData()
    x1=res[[1]]
    y1=res[[2]]
    t1=res[[3]]
    if ((input$ThinData==T)&&(!is.null(input$thinOutTo))) {               #thin-out of the data
      t1new=seq(min(t1),max(t1),length.out=input$thinOutTo)
      y1func=suppressWarnings(splinefun(t1,y1,method = "fmm"))
      y1=y1func(t1new)
      x1func=suppressWarnings(splinefun(t1,x1,method = "fmm"))
      x1=x1func(t1new)
      t1=t1new }
    endX=length(x1) 
    xLeft=ifelse(is.null(input$sta),(x1[endX]-x1[1])*0.1+x1[1],input$sta)
    xRight=ifelse(is.null(input$end),(x1[endX]-x1[1])*0.9+x1[1],input$end)
    BLind=which((x1<xLeft)|(x1>xRight))
    xBL=x1[BLind]
    if (input$BLsmooth==T) yBL=predict(loess(y1[BLind]~xBL, span = 0.1))  #smoothing of data used for baseline construction
    else yBL=y1[BLind]
    if ((all(is.na(yBL)))&&(all(is.na(xBL)))) BLline=0
    else {
      if (input$BLtype=="line") {                                         #linear baseline substraction
        BLfit=lm(yBL~xBL)
        BLline=BLfit$coefficients[1]+BLfit$coefficients[2]*x1 }                          
      if (input$BLtype=="spln") {                                         #spline baseline substraction
        BLspl= suppressWarnings(splinefun(xBL, yBL,method = "natural"))
        BLline=BLspl(x1) }     
      if (input$BLtype=="linst") {
        BLline=y1
        BLind.left=which(x1<xLeft)
        BLfit=lm(y1[BLind.left]~x1[BLind.left])
        BLline[BLind.left]=BLfit$coefficients[1]+BLfit$coefficients[2]*x1[BLind.left]
        BLind.right=which(x1>xRight)
        BLfit=lm(y1[BLind.right]~x1[BLind.right])
        BLline[BLind.right]=BLfit$coefficients[1]+BLfit$coefficients[2]*x1[BLind.right]
        k.line=(BLline[min(BLind.right)]-BLline[max(BLind.left)])/(x1[min(BLind.right)]-x1[max(BLind.left)])
        BLline[-BLind]=k.line*x1[-BLind]+(BLline[min(BLind.right)]-k.line*x1[min(BLind.right)]) } 
      if (input$BLtype=="spln2") {
        BLline=y1
        BLind.left=which(x1<xLeft)
        BLfit=lm(y1[BLind.left]~x1[BLind.left])
        BLline[BLind.left]=BLfit$coefficients[1]+BLfit$coefficients[2]*x1[BLind.left]
        BLind.right=which(x1>xRight)
        BLfit=lm(y1[BLind.right]~x1[BLind.right])
        BLline[BLind.right]=BLfit$coefficients[1]+BLfit$coefficients[2]*x1[BLind.right]
        BLspl= suppressWarnings(splinefun(x1[BLind], BLline[BLind],method = "natural"))
        BLline=BLspl(x1) } 
      }
    spl0=suppressWarnings(splinefun(x1,y1,method = "fmm"))
    y1=y1-BLline
    yPts=length(y1)
    y1min=min(y1)
    posMin=which.min(y1)
    spl1=suppressWarnings(splinefun(x1,y1,method = "fmm"))
    bl=y1[1]
    dy1=spl1(x1[1:posMin],deriv=1)
    pos=which.min(dy1) 
    HRind=which((x1>xLeft)&(x1<xRight))
    Ttfit=lm(x1[HRind]~t1[HRind])
    HR=Ttfit$coefficients[2]                                               #heating rate in [°C/s]
    dH=simpson(spl1,x1[1],x1[yPts],input$N_Simp)            #heat effect in W*°C/g or mV*°C/g
    linInd=which((y1<(y1min-bl)*0.05)&(y1>y1min*0.95))
    linInd=linInd[linInd<posMin]
    nearPercent=ifelse(input$TRtype=="near15",0.15,0.05)                   #adjust this number if needed
    deltaPos=ceiling(length(linInd)*nearPercent)                           
    posDeltamin=ifelse((pos-deltaPos)>1,pos-deltaPos,1)
    posDeltamax=ifelse((pos+deltaPos)>length(dy1),length(dy1),pos+deltaPos)
    if(input$TRtype=="min") dy1select=dy1[pos]                             #thermal resistance in point where gradient minimal
    if((input$TRtype=="near5")||(input$TRtype=="near15")) 
        dy1select=mean(dy1[posDeltamin:posDeltamax])                       #thermal resistance near point where gradient minimal
    if(input$TRtype=="ave5") {
      TRfit=lm(y1[linInd]~x1[linInd])
      dy1select=TRfit$coefficients[2] }               #thermal resistance in 5-95% region 
    if(input$TRtype=="ave10") {
      linInd=which((y1<(y1min-bl)*0.1)&(y1>y1min*0.9))
      linInd=linInd[linInd<posMin]
      TRfit=lm(y1[linInd]~x1[linInd])
      dy1select=TRfit$coefficients[2] }               #thermal resistance in 10-90% region 
    if(input$TRtype=="us") #thermal resistance user-specified
      dy1select=-1/(input$usTR*input$mass1)
    
    TRselect=1/(abs(dy1select)*input$mass1)
    meanPos=round((linInd[1]+linInd[length(linInd)])/2)
    x1pos=ifelse(((input$TRtype=="near5")||(input$TRtype=="near15")||(input$TRtype=="min")),x1[pos],x1[meanPos])
    y1pos=ifelse(((input$TRtype=="near5")||(input$TRtype=="near15")||(input$TRtype=="min")),y1[pos],y1[meanPos])
    yint=spl1(x1pos)-(dy1select*x1pos)
    xint=(bl-yint)/dy1select                                                #extrapolated onset in [°C]
    yLine=yint+dy1select*x1                                                 #coordinates for line for extrapolated onset
    spl2=function(z){z*spl1(z)}
    dTH=simpson(spl2,x1[1],x1[yPts],input$N_Simp) 
    Tint=dTH/dH
    xCorr=x1-y1/dy1select                                                   #corrected by selected thermal resistance
    spl2=suppressWarnings(splinefun(xCorr,y1,method = "fmm"))
    #block below for time constant, or the DSC signal return to baseline, calculation. Will be finished!
    #xReggr=seq(xCorr[posMin],xCorr[yPts],length.out=100)
    #tReggr=seq(t1[posMin],t1[yPts],length.out=100)
    posRight=which(x1>xRight)[1]
    yReggr=y1[posMin:posRight]
    xReggr=x1[posMin:posRight] #or Xcorr as above?
    tReggr=t1[posMin:posRight]-t1[posMin]
#    if (input$TimeConst==T) {
#        nlrRes=nls(yReggr ~ a1*exp(-tReggr/tau1), start=list(a1=-1,tau1=5))
#        #TBD. or two time constants, as below?
#        #nlrRes=nls(yReggr ~ a1*exp(-tReggr/tau1)+a1*exp(-tReggr/tau2), start=list(a1=coef(nlrRes1)[1],tau1=coef(nlrRes1)[2],tau2=30))
#    }
#    else 
    nlrRes=1 #just not to be undefined
    
    #spl5=splinefun(t1,y1,method = "fmm")
    #dH2=simpson(spl5,t1[1],t1[yPts]) 
    HRind_left=HRind[1]:posMin
      
    list(x1,y1,bl,spl1(x1),x1pos,y1pos,yLine,xint,HR,dH/HR,
         TRselect,HRind_left,Tint,spl1(Tint),xCorr,xReggr,spl2(xReggr),t1[posMin:posRight],spl2(xCorr),nlrRes,
         t1,spl0(input$sta),spl0(input$end),y1+BLline,BLline)
  }) 

  purityAnalysis <- reactive({#----------------------------Computes the purity of the sample based on DSC of melting 
    if (input$go == 0) return()
    res=treatData()
    HRind=res[[12]] #indices of data from left baseline border points to peak maximum
    x1=res[[1]]+273.15     #temperature in K
    xCorr=res[[15]]+273.15 #temperature in K corrected to thermal lag
    y1=res[[2]]     #DSC with baseline substracted
    t1=res[[21]]    #time in seconds
    HE=res[[10]]    #heat effect in J/g
    HR=res[[9]]
    Area=HE*HR
    Tint=res[[13]]
    yS=y1[HRind]
    xS=x1[HRind]
    tS=t1[HRind]
    pts=length(tS)
    A=yS #to define an array
    spl4=suppressWarnings(splinefun(xS,yS,method = "natural"))
    for (m in 1:pts) A[m]=simpson(spl4,min(xS),xS[m],ceiling(m*input$N_Simp/pts))
    #A=suppressWarnings(as.numeric(as.character(A)))
    if(input$PurCorr=="y") xS=xCorr[HRind]
    Fn=A/Area
    #HRind_lim=which((Fn>0.005)&(Fn<0.995))
    #Fn=Fn[HRind_lim]
    #xS=xS[HRind_lim]
    #A=A[HRind_lim]
    #write.table(cbind(tS,xS,yS,Fn), file = "line2.txt", sep=" ", row.names=FALSE, quote = FALSE)
    if (input$task=="pa")  {
      HRind_F=which((Fn>(input$Flim[1]/100))&(Fn<(input$Flim[2]/100)))
      TF=xS[HRind_F]#-273.15
      Hp=A[HRind_F] 
      Fn_sel=Fn[HRind_F]
      xS=xS[HRind_F]
      yS=yS[HRind_F]
      if (input$PurrThin=="y") {
          oFn_thin=seq(min(1/Fn_sel),max(1/Fn_sel),length.out=input$PT_pts)
          TFfunc=suppressWarnings(splinefun(1/Fn_sel,TF,method = "natural"))
          TF=TFfunc(oFn_thin) 
          Hp=Area/oFn_thin
          Fn_sel=1/oFn_thin }
      #write.table(cbind(TF,Area,Hp,1/Fn_sel), file = "line3.txt", sep=" ", row.names=FALSE, quote = FALSE)
      #lmRes=nls(TF ~ To+slope*(Area+c)/(Hp+c), start=list(To=Tint, slope=-0.2,c=0))
      #oFn_corr=(Area+coef(lmRes)[3])/(Hp+coef(lmRes)[3])
      #Area_corr=coef(lmRes)[3]/Area*100
      #To=coef(lmRes)[1]
      #Purity=(1-(coef(lmRes)[2]*(Area+coef(lmRes)[3])/HR*input$MolarMass/8.314/(coef(lmRes)[1]+273.15)^2))*100
      lmRes=lm(TF*Hp ~ Hp+TF)
      oFn_corr=(Area-lmRes$coefficients[3])/(Hp-lmRes$coefficients[3])
      To=lmRes$coefficients[2]#-273.15
      slope=(lmRes$coefficients[1]+lmRes$coefficients[3]*lmRes$coefficients[2])/(Area-lmRes$coefficients[3])
      Area_corr=-lmRes$coefficients[3]/Area*100
      HE_corr=(Area-lmRes$coefficients[3])/HR
      Purity=(1-(slope*HE_corr*input$MolarMass/8.314/(lmRes$coefficients[2])^2))*100
      lmRes_Summary=summary(lmRes)
      std_Err=qt(0.95,df=lmRes$df.residual)*lmRes_Summary$coefficients[,2]
      #dPurity=100*input$MolarMass/8.314/HR*((std_Err[1]/(lmRes$coefficients[2])^2)^2 +
      #                                         (std_Err[2]*(2*lmRes$coefficients[1]/(lmRes$coefficients[2])^3+lmRes$coefficients[3]/(lmRes$coefficients[2])^2))^2 + 
      #                                         (std_Err[3]/lmRes$coefficients[2])^2 )^0.5
      lmRes2=lm(TF-273.15~oFn_corr)
      lmRes2_Summary=summary(lmRes2)
      P2=100-100*lmRes2$coefficients[2]*HE_corr*input$MolarMass/8.314/(To^2)
      dPurity=( (qt(0.95,df=lmRes2$df.residual)*lmRes2_Summary$coefficients[2,2]/lmRes2$coefficients[2])^2 + 
                (qt(0.95,df=lmRes2$df.residual)*lmRes2_Summary$coefficients[1,2]/(lmRes2$coefficients[2]+273.15))^2 )^0.5*(100-Purity)
      To2=lmRes2$coefficients[1]
    }
    #if (input$task=="pb") {
    #  Nrange=11
    #  yFactor=seq(input$Flim[1],input$Flim[2],length.out = Nrange)
    #  tP=numeric(Nrange)
    #  xP=numeric(Nrange)
    #  yP=numeric(Nrange)
    #  Hp=numeric(Nrange)
    #  for (k in 1:Nrange) {
    #    yP[k]=ySmin*yFactor[k]/100
    #    xP[k]=spl3(yP[k])
    #    tP[k]=spl6(yP[k])
    #    Hp[k]=simpson(spl4,tS[1],tP[k],floor(k/Nrange*input$N_Simp))*input$mass1
    #  }
    #  oxP=1/xP
    #  HxP=Hp/xP
    #  lmRes=lm(Hp ~ oxP+HxP)
    #  cat(y1min, file = "line.txt")
    #  write.table(cbind(yFactor,tP,xP,yP,Hp,oxP,HxP), file = "line.txt", sep=" ", row.names=FALSE, quote = FALSE)
    #  Purity=(1-(lmRes$coefficients[3]*lmRes$coefficients[1]-lmRes$coefficients[2])*input$MolarMass/(8.314*lmRes$coefficients[3]^2*input$mass1))*100
    #}

    list (Purity,1/Fn_sel,TF-273.15,oFn_corr,TF-273.15,lmRes_Summary,Area_corr,1/Fn,xS-273.15,To-273.15,
          std_Err,HE_corr,std_Err[3]/HR,dPurity,lmRes2,P2,To2,yS)
  })

  output$textPurity <- renderPrint({ #returns full regression output if selected
    if (is.null(input$file1)) return(cat("Load file for analysis"))
    if (input$go == 0) return(cat("Press \"Start\""))
    if ((input$task!="pb")&&(input$task!="pa")) return(cat("Select the appropriate task"))
    else {
      std_Err=purityAnalysis()[[11]]
      Purity=purityAnalysis()[[1]]
      dPurity=purityAnalysis()[[14]]
      Correction=purityAnalysis()[[7]]
      To=purityAnalysis()[[10]]
      HE=purityAnalysis()[[12]]
      dHE=purityAnalysis()[[13]]
      cat(c("","/Summary of analysis/","\n",
            "Purity",round(Purity,digits=3),"±",formatC(dPurity, format="e",digits=2),"mol.%", #round(purityAnalysis()[[16]],digits=3),
            "with correction",round(Correction,digits=1),"%","\n",
            "Melting temperature",round(To,digits=2),"±",round(std_Err[2],digits=2),"°C","\n", #round(purityAnalysis()[[17]],digits=3),
            "Corrected heat effect",round(HE,digits=2),"±",round(dHE,digits=2),"J/g (",
            round(HE*input$MolarMass/1000,digits=1),"±",round(dHE*input$MolarMass/1000,digits=1),"kJ/mol )"
            ),sep=" ") }
  })
  
  output$distPlot <- renderPlot({
    if (is.null(input$file1)) return(NULL)
    if (input$go == 0) return()
    res=treatData()
    DSCunits=ifelse((as.numeric(input$dscUnits)==1)||(as.numeric(input$dscUnits)==3),"mW/mg","mkV/mg")
    par (mar=c(4.1,5.1,3.1,0.3))
    plot(res[[1]],res[[2]], xlab="Temperature, °C", ylab=paste("DSC,",DSCunits), main="Corrected DSC data",
         cex.lab=1.5) 
    #    title(main="DSC data with substracted baseline (hollow points) / Corrected to thermal resistance (blue line) /
    #              Extrapolated onset (green point) / Point where gradient is maximal (blue point) /
    #              Average integral temperature (red point)", cex.main = 1,   font.main= 1)
    abline(h=res[[3]], col=8)        #horizontal baseline
    lines(res[[1]],res[[4]])
    lines(res[[1]],res[[7]], col=3)
    #block below adds text on the graph with real heating rate inside the analyzed region and heat effect (area of the analyzed peak)
    usr=par("usr")
    mess1=paste("Heating rate:",round(res[[9]]*60,digits=2),"K/min")
    xCoord=(usr[2]-usr[1])*1/30+usr[1]
    yCoordmin=0.55
    yCoordStep=0.05
    text(xCoord,(usr[3]-usr[4])*yCoordmin+usr[4], mess1, pos=4)
    Heunits=ifelse((as.numeric(input$dscUnits)==1)||(as.numeric(input$dscUnits)==3),"J/g","mV*s/g")
    mess2=paste("Heat effect:",round(res[[10]],digits=2),Heunits)
    text(xCoord,(usr[3]-usr[4])*(yCoordmin+yCoordStep)+usr[4], mess2, pos=4)
    #following data is shown inside the legend area
    mess3=paste("Extrapolated onset:",round(res[[8]],digits=2),"°C")
    points(res[[8]],res[[3]], col=3, pch=19, cex=2) #point denoting the extrapolated onset
    TRunits=ifelse((as.numeric(input$dscUnits)==1)||(as.numeric(input$dscUnits)==3),"K/mW","K/mkV")
    mess4=paste("Thermal resistance:",round(res[[11]],digits=3),TRunits)
    points(res[[5]],res[[6]], col=2, pch=19, cex=2) #point where gradient minimal
    mess5=paste("Average integral temperature:",round(res[[13]],digits=2),"°C")
    points(res[[13]],res[[14]], col=4, pch=19, cex=2) #point with average integral temperature
    #text(xCoord,(usr[3]-usr[4])*(yCoordmin+yCoordStep*2)+usr[4], mess3, pos=4)
    #text(xCoord,(usr[3]-usr[4])*(yCoordmin+yCoordStep*3)+usr[4], mess4, pos=4)
    #text(xCoord,(usr[3]-usr[4])*(yCoordmin+yCoordStep*4)+usr[4], mess5, pos=4)  
    lines(res[[15]],res[[2]], col=28) #draw corrected by thermal resistance DSC
    #TBD. Will add the time constant correction...
    #text(xCoord,(usr[3]-usr[4])*(yCoordmin+yCoordStep*5)+usr[4], "Time const", pos=4)    
    #if (input$TimeConst==T) lines(res[[16]],predict(res[[20]]), col=3) #now it only draws the exponential fit, TBD
    legend(xCoord,(usr[3]-usr[4])*(yCoordmin+yCoordStep*1.1)+usr[4], bty="n", #fill = NULL,
           c("DSC with substracted baseline","DSC corrected to thermal resistance",mess3,mess4,mess5), 
           col=c(1,28,3,2,4), lty = c(NA,1,NA,NA,NA), pch = c(1, NA, 19, 19,19))
  })
  
  output$purityPlot2 <- renderPlot({
    if (is.null(input$file1)) return(NULL)
    if (input$go == 0) return()
    if ((input$task!="pb")&&(input$task!="pa")) return(NULL)
    res=treatData()
    res2=purityAnalysis()
    DSCunits=ifelse((as.numeric(input$dscUnits)==1)||(as.numeric(input$dscUnits)==3),"mW/mg","mkV/mg")
    par (mar=c(4.1,5.1,3.1,0.3))
    plot(res[[1]],res[[2]], xlab="Temperature, °C", ylab=paste("DSC,",DSCunits), cex.lab=1.5,
         main="DSC with indicated area range used for purity analysis") 
    
    lines(res[[1]],res[[4]])
    lines(res[[15]],res[[2]], col=28) #draw corrected by thermal resistance DSC
    min_x=min(res2[[9]])
    max_x =max(res2[[9]])
    max_y =max(res2[[18]]) 
    polygon(c(res2[[9]],max_x,min_x,min_x),
            c(res2[[18]],0,0,max_y),col=rgb(1, 0, 0,0.5), border=NA)
    
    #points(res2[[9]],res2[[18]], col=2, pch=19, cex=2) #points used for purity calculation
    
  })
    
  output$textPurity2 <- renderPrint({ #returns full regression output if selected
    if (is.null(input$file1)) return("Load file for analysis")
    if (input$go == 0) return("Press \"Start\"")
    if ((input$task!="pb")&&(input$task!="pa")) return("Select the appropriate task")
    else {
      res=purityAnalysis()[[6]]
      #write.table(as.data.frame(summary(res)), file = "line1.txt")
      res
      }
  })
  
  output$PurityPlot <- renderPlot({
    if (is.null(input$file1)) return(NULL)
    if (input$go == 0) return()
    if ((input$task!="pb")&&(input$task!="pa")) return(NULL)
    res=purityAnalysis()
    oFn=res[[2]]
    TF=res[[3]]
    par (mar=c(5.1,5.1,3.1,0.3))
    plot(oFn,res[[3]], xlab="1/F", ylab="Temperature, °C", xlim=c(0,max(oFn)), ylim=c(min(TF),res[[10]]), 
         main="Purity Analysis according to ASTM E928-08", cex.lab=1.5) 

    #lines(res[[1]],res[[7]], col=3)
   
    points(res[[4]],res[[5]], col=3, pch=19, cex=2)
    abline(res[[15]])

  })
  
  output$distPlot2 <- renderPlot({#-------------------Draw the bottom plot for baseline adjustment-------------------------
    if (is.null(input$file1)) return(NULL)
    if (input$go == 0) return()
    res=treatData()
    dscUnc=res[[24]]
    dPts=length(dscUnc) 
    dMin=min(dscUnc[1],dscUnc[dPts])
    dMax=max(dscUnc[1],dscUnc[dPts])
    yMin=min(res[[24]])
    yMax=max(res[[24]])
    DSCunits=ifelse((as.numeric(input$dscUnits)==1)||(as.numeric(input$dscUnits)==3),"mW/mg","mkV/mg")
    #par(fig = c(0,1,0,1))
    par (mar=c(4.1,5.1,4.1,0.3))
    plot(res[[1]],res[[24]],
         ylim=c(min(dMin,dMin-(yMax-yMin)*0.07),max(dMax,dMax+(yMax-yMin)*0.03)), #specifies the zoomed window for data view
         xlab="Temperature / °C", ylab=paste("DSC",DSCunits), main="Raw DSC data / Baseline adjustment",
         cex.lab=1.5)
    lines(res[[1]],res[[25]], col=4)
    points(input$sta,res[[22]], col=28, pch=19, cex=2) #denoted right side of initial region for baseline construction
    points(input$end,res[[23]], col=28, pch=19, cex=2) #denoted left side of final region for baseline construction
    usr=par("usr")
    xCoord=(usr[2]-usr[1])*1/30+usr[1]
    yCoordmin=0.55
    legend(xCoord,(usr[3]-usr[4])*yCoordmin+usr[4], bty="n", 
           c("raw DSC data","Left/Right limits for baseline","Baseline"), 
           col=c(1,28,4), lty = c(NA,NA,1), pch = c(1, 19,NA))
    
    if (input$ThinData==T) {#-------------------Draw the inset plot for thinning-out adjustment-------------------------
      resIn=readData()
      par(fig = c(0.5,1, 0, 0.6), new = T)
      plot(resIn[[1]],resIn[[2]], type="n",ylab='',xlab='',xaxt="n",yaxt="n")
      rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col=gray(0.9))
      points(resIn[[1]],resIn[[2]])
      lines(res[[1]],res[[24]], col=4)
      usr=par("usr")
      xCoord=(usr[2]-usr[1])*1/30+usr[1]
      yCoordmin=0.55
      legend("bottomleft", bty="n", cex=1,
             c("raw DSC data","thin-out curve"), 
             col=c(1,4), lty = c(NA,1), pch = c(1, NA))
      }
    
    })

  
})
