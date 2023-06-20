#THINKS thermokinetic software, version #1.06.20 Website: thinks.chemphys.ru. To cite: 10.3390/molecules24122298
#This routine is used to import the experimental data into the kinetic project.
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
library(takos)

options(shiny.maxRequestSize=5*1024^2)  #change the file-size limit here----in the server there is 10 MB limit

shinyServer(function(input, output,session) {

  session$onSessionEnded(function() { 
    observe({ 
      stopApp() })
  })
  
#  output$cond2 <- reactive({ #----condition for small data (e.g., with the smoothing errors)
#    inFile <- input$file1
#    if (is.null(inFile)) return(FALSE)
#    res=readData()
#    cat(IQR(res[[2]]),file="iqr.txt",append=TRUE)
#    IQR(res[[2]]) > 0 
#  })
#  outputOptions(output, "cond2", suspendWhenHidden = FALSE)
  
  output$cond <- reactive({
    is.null(input$file1)
  })
  outputOptions(output, "cond", suspendWhenHidden = FALSE)
  
  observeEvent(input$file1, {
    updateTabsetPanel(session, "MainTabset",
                      selected = "panel1")
  }, once=TRUE)
  
  simpson <- function(fun,a,b,n) { #calculation of the integral using Simpson equation, checked some special functions, but this works fact and accurate
      n=ifelse(n<10,10,n)
      h=(b-a)/n
      x=seq(a, b, by=h)
      if (a==b) return (0)
      if (n==2) s=fun(x[1])+4*fun(x[2])+fun(x[3])  
      else s=fun(x[1])+fun(x[n+1])+2*sum(fun(x[seq(2,n,by=2)]))+4 *sum(fun(x[seq(3,n-1, by=2)])) 
      s*h/3 }
  
  readData <- reactive({ #only reads data from file, without further calculations
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    InputData=read.table(inFile$datapath,sep=input$sep1, skip=input$SkipL)              #read file defined in "file1" with separators "sep1" and skip X lines
    removeLines=nrow(InputData)-input$EndSkipL
    InputData=InputData[1:removeLines,]
    Temp=InputData[,as.numeric(input$tempCol)]-273.15*(as.numeric(input$tempUnits)-1)   #temperature in [°C]
    ti=InputData[,as.numeric(input$timeCol)]*as.numeric(input$timeUnits)                #time in [s]
    if (input$signal=="int") y1=InputData[,as.numeric(input$tgaCol)]
    else y1=InputData[,as.numeric(input$dscCol)]
    list(Temp,ti,y1) 
  })
  
  output$BLpointsBeg <- renderUI({ #generates slider for left point selection for baseline generation or mass loss correction
    if ((is.null(input$file1))||(as.numeric(input$tgaUnits)>2)) return(NULL)
    res=readData()
    if (input$tprg=="niso") Tx=res[[1]]
    else Tx=res[[2]]
    endTx=length(Tx) 
    Tmin=round(Tx[1]*10)/10
    Tmax=round(Tx[endTx]*10)/10
    Tval=round(((Tx[endTx]-Tx[1])*0.1+Tx[1])*10)/10
    sliderInput("sta", "Starting point :",min = Tmin, max = Tmax, value = Tval, step=0.05)
  })
  
  output$BLpointsEnd <- renderUI({ #generates slider for right point selection for baseline generation or mass loss correction
    if ((is.null(input$file1))||(as.numeric(input$tgaUnits)>2)) return(NULL)
    res=readData()
    if (input$tprg=="niso") Tx=res[[1]]
    else Tx=res[[2]]
    endTx=length(Tx) 
    Tmin=round(Tx[1]*10)/10
    Tmax=round(Tx[endTx]*10)/10
    Tval=round(((Tx[endTx]-Tx[1])*0.9+Tx[1])*10)/10
    sliderInput("end", "Final point :",min = Tmin, max = Tmax, value = Tval, step=0.05)
  })

  output$BLpointsCntr <- renderUI({ #generates slider for central point selection for baseline generation 
    if ((is.null(input$file1))||(as.numeric(input$tgaUnits)>2)) return(NULL)
    res=readData()
    if (input$tprg=="niso") Tx=res[[1]]
    else Tx=res[[2]]
    endTx=length(Tx) 
    Tmin=round(Tx[1]*10)/10
    Tmax=round(Tx[endTx]*10)/10
    Tval=round(((Tx[endTx]-Tx[1])*0.5+Tx[1])*10)/10
    sliderInput("cntr", "Central T :",min = Tmin, max = Tmax, value = Tval, step=0.05)
  })
  
  output$Yscale <- renderUI({ #generates slider for Y scale
    if ((is.null(input$file1))||(input$signal=="int")||(input$BLtype=="no")) return(NULL)
    res=isolate(treatData())
    dscUnc=res[[12]]
    dPts=length(dscUnc) 
    dMin=min(dscUnc[1],dscUnc[dPts])
    dMax=max(dscUnc[1],dscUnc[dPts])
    yMin=min(res[[12]])
    yMax=max(res[[12]])
    minAppr=signif(min(dMin,min(yMax,dMin)-(yMax-yMin)*0.05), digits=4)
    maxAppr=signif(max(dMax,min(yMax,dMax)+(yMax-yMin)*0.05), digits=4)
    sliderInput("YscaleRange", "Y-scale :", 
                min = min(yMin,minAppr-(maxAppr-minAppr)*0.2), max = max(yMax,maxAppr+(maxAppr-minAppr)*0.5), 
                #min=yMin , max=yMax , #min = minAppr-(maxAppr-minAppr)*0.2, max = maxAppr+(maxAppr-minAppr)*0.5, 
                value = c(minAppr,maxAppr))
  })
  
  treatData <- reactive({ #converts data from file to conversion, conversion rate and its derivatives
    res=readData()
    y1=res[[3]]  # DSC (TGA) in as loaded units
    if (input$tprg=="niso") {
        x1=res[[1]]     # Temperature in °C if nonisothermal run
        t1=res[[2]]     # Time in s
        xS=t1
        TS=x1
        }  
    else {
        x1=res[[2]]     # Time in s if isothermal
        t1=res[[1]] 
        xS=x1
        TS=t1
        }
    endX=length(x1) 
    if (input$smType=="cub") y1=smooth.spline(xS,y1,spar=input$CubSpar)$y                 #two types of smoothing
    if (input$smType=="pol") y1=predict(loess(y1~xS,span=10^input$PolSpan), data.frame(x=xS))
    xLeft=ifelse(is.null(input$sta),(x1[endX]-x1[1])*0.1+x1[1],input$sta)
    xCenter=ifelse(is.null(input$cntr),(x1[endX]-x1[1])*0.5+x1[1],input$cntr)
    xRight=ifelse(is.null(input$end),(x1[endX]-x1[1])*0.9+x1[1],input$end)
    NPTS=ifelse(input$Npts>2000,2000,input$Npts) #---------------------this limit can be extended!
    if (input$signal=="dif") {
        if (input$BLtype!="no") BLind=which((x1<xLeft)|(x1>xRight))
        if (input$BLtype=="horRight") {
          BLind.right=which(x1>xRight)
          BLfit=mean(y1[BLind.right])
          BLline=BLfit+0*x1 }                                                             #horizontal-right baseline substraction
        if (input$BLtype=="line") {
            BLfit=lm(y1[BLind]~x1[BLind])
            BLline=BLfit$coefficients[1]+BLfit$coefficients[2]*x1 }                       #linear baseline substraction
        if (input$BLtype=="spln") {
            BLspl= suppressWarnings(splinefun(x1[BLind], y1[BLind],method = "natural"))
            BLline=BLspl(x1) }                                                            #spline baseline substraction
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
        if (input$BLtype=="spln2c") {
          BLline=y1
          BLind.left=which(x1<xLeft)
          BLfit=lm(y1[BLind.left]~x1[BLind.left])
          BLline[BLind.left]=BLfit$coefficients[1]+BLfit$coefficients[2]*x1[BLind.left]
          BLind.right=which(x1>xRight)
          BLfit=lm(y1[BLind.right]~x1[BLind.right])
          BLline[BLind.right]=BLfit$coefficients[1]+BLfit$coefficients[2]*x1[BLind.right]
          deltaX=x1[2]-x1[1]
          Ampl=x1[endX]-x1[1]
          #N_cntr_pts=ceiling(min(length(BLind.left),length(BLind.right))*input$FNcentrPts/100)
          #cat(as.numeric(input$FNcentrPts),"\n", file="a2.txt", append = T)
          prct=ifelse(is.na(input$FNcentrPts)||is.null(input$FNcentrPts),1,min(input$FNcentrPts,30))
          gap=max(Ampl*prct/100,5*deltaX)
          BLind.center=which((x1<(xCenter+gap))&(x1>(xCenter-gap)))
          BLline[BLind.center]=min(y1[BLind.center])
          BLind=c(BLind.left,BLind.center,BLind.right)
          fit.ns=lm(y~poly(x, degree=3), data=data.frame(x=x1[BLind], y=BLline[BLind]))
          BLline=predict(fit.ns, newdata=list(x=x1)) }
        if (input$BLtype=="tang") {
          BLind.left=which(x1<xLeft)
          BLind.right=which(x1>xRight)
          BLline=TAPPA(x1,y1, interval=min(length(BLind.left),length(BLind.right)))  } 
        
        if (input$BLtype=="glasstr") {
          
          BLind.left=which(x1<xLeft)
          BLfit1=lm(y1[BLind.left]~x1[BLind.left])
          
          #BLline[BLind.left]=BLfit$coefficients[1]+BLfit$coefficients[2]*x1[BLind.left]
          
          BLind.right=which(x1>xRight)
          BLfit2=lm(y1[BLind.right]~x1[BLind.right])
          
          #BLline=BLfit2$coefficients[1]+BLfit2$coefficients[2]*x1
          BLline=y1*0
          } 
               
        if (input$BLtype!="no") {
            spl0=suppressWarnings(splinefun(x1,y1,method = "fmm"))
            y1=y1-BLline }
        if (input$DSCLess0==T) y1[y1<0]=0

        spl1=suppressWarnings(splinefun(xS,y1,method = "fmm"))
        N_SIMP=ifelse(input$N_Simp>2000,2000,input$N_Simp)
        al=double(NPTS)
        ti=seq(min(xS),max(xS),length.out=NPTS)
        Tem=suppressWarnings(splinefun(xS,TS,method = "fmm"))
        Te=Tem(ti)
        
        if (input$tprg=="niso") {
            Ttfit=lm(x1~t1)
            Area=simpson(spl1,xS[1],xS[endX],N_SIMP) } 
        else { 
            Ttfit=lm(t1~x1)
            Area=simpson(spl1,xS[1],xS[endX],N_SIMP) }
        HR=Ttfit$coefficients[2]                                                       #heating rate in [°C/s]

        if (input$BLtype!="glasstr") {
            dH=Area/ifelse(as.numeric(input$dscUnits)==2,input$mass1,1)                #heat effect in J/g, J/g or [-]
            da=spl1(ti)/Area
            for (m in 1:NPTS) al[m]=simpson(spl1,min(ti),ti[m],ceiling(m*N_SIMP/NPTS))/Area*as.numeric(input$afin)
            al=suppressWarnings(as.numeric(as.character(al))) }
        else {
            dH=0
            for (m in 1:NPTS) {
              Fg=BLfit1$coefficients[1]+BLfit1$coefficients[2]*Te[m]
              al[m] = (spl0(Te[m])-Fg)/((BLfit2$coefficients[1]+BLfit2$coefficients[2]*Te[m])-Fg) }
            alphGT=suppressWarnings(splinefun(ti,al,method = "fmm"))
            da=alphGT(ti,deriv=1)
        }
        if (input$aLess0==F) al[al<1E-6]=0
        if (input$aHigher1==F) al[al>(1-1E-6)]=1       
        if (input$BLtype=="no") return(list(x1,y1,t1,y1[1],dH,round(HR*60,digits=2),ti,Te,al,da,da*Area))
        else                    return(list(x1,y1,t1,y1[1],dH,round(HR*60,digits=2),ti,Te,al,da,da*Area,
                                            y1+BLline,BLline,spl0(xLeft),spl0(xRight),spl0(xCenter)))
        }
    else {
      StaInd=min(which(x1>xLeft))
      EndInd=min(which(x1>xRight))
      if (as.numeric(input$tgaUnits)==4) y0=y1/100
      if (as.numeric(input$tgaUnits)==3) y0=y1
      if (as.numeric(input$tgaUnits)<3) {
        ML=(y1[StaInd]-y1[EndInd])/ifelse(as.numeric(input$tgaUnits)==2,input$mass1/100,1)   #mass loss in %, % or [-]
        y0=(y1[StaInd]-y1)/(y1[StaInd]-y1[EndInd]) }
      else ML=-y0[StaInd]+y0[EndInd]
      spl0=suppressWarnings(splinefun(x1,y1,method = "fmm"))
      y0=y0*as.numeric(input$afin)
      if (input$tprg=="niso") {
        Ttfit=lm(x1~t1)
        ti=seq(min(t1),max(t1),length.out=NPTS)
        Tem=suppressWarnings(splinefun(t1,x1,method = "fmm"))
        Te=Tem(ti)
        alph=suppressWarnings(splinefun(t1,y0,method = "fmm"))
        al=alph(ti)
        da=alph(ti,deriv=1) }
      else { 
        Ttfit=lm(t1~x1)
        ti=seq(min(x1),max(x1),length.out=NPTS)
        Tem=suppressWarnings(splinefun(x1,t1,method = "fmm"))
        Te=Tem(ti) 
        alph=suppressWarnings(splinefun(x1,y0,method = "fmm"))
        al=alph(ti)
        da=alph(ti,deriv=1) }
      HR=Ttfit$coefficients[2]                                                            #heating rate in [°C/s]
      al=suppressWarnings(as.numeric(as.character(al)))
      if (input$aLess0==F) al[al<1E-6]=0
      if (input$aHigher1==F) al[al>(1-1E-6)]=1     
      if (as.numeric(input$tgaUnits)>2) 
        if (input$tprg=="niso")
            return(list(x1,y1,t1,"",ML,round(HR*60,digits=2),ti,Te,al,da,spl0(Te)))
        else return(list(t1,y1,x1,"",ML,round(HR*60,digits=2),ti,Te,al,da,spl0(ti)))
      else                              
        if (input$tprg=="niso")
            return(list(x1,y1,t1,"",ML,round(HR*60,digits=2),ti,Te,al,da,spl0(Te),
                                                spl0(input$sta),spl0(input$end),StaInd,y1[StaInd],EndInd,y1[EndInd]))
        else 
            return(list(t1,y1,x1,"",ML,round(HR*60,digits=2),ti,Te,al,da,spl0(ti),
                      spl0(input$sta),spl0(input$end),StaInd,y1[StaInd],EndInd,y1[EndInd]))
      }
  })

  logfilename <- paste0('input',floor(runif(1, 0, 1e+06)),'.txt')               #name of logfile with project kinetic data #originally: 1e+05, 1e+06 - 1
  
  writeKinData <- reactive({ 
    input$go1                                                                           #do only when button pressed!
    inFile <- isolate(input$file1) 
    if (is.null(inFile)) return(NULL)
    res=isolate(treatData())
    Fname=gsub(" ", "", inFile$name)
    FN=ifelse(nchar(Fname)>15,substr(Fname,1,15),Fname)
    mess1=paste("HR",res[[6]],"Kmin", sep="")
    if (isolate(input$signal)=="dif") {
      numb=abs(res[[5]])
      dHshow=round(numb/10^(trunc(log10(numb))-2))*10^(trunc(log10(numb))-2)*ifelse(res[[5]]<0,-1,1)
      dHunits=ifelse(as.numeric(isolate(input$dscUnits))==3,"-","Jg")
      mess2=paste("HE",dHshow,dHunits, sep="") }
    else {
      dMshow=round(res[[5]]*100)/100
      dMunits=ifelse(as.numeric(isolate(input$tgaUnits))>2,"-","%")
      mess2=paste("ML",dMshow,dMunits, sep="") }
    if (!is.null(inFile$name)) {
      if (is.na(file.size(logfilename))) {
      #write data to file (4 columns: T[°C], t [s], alpha [-], da/dt [s])
      #column names: filename(or truncated), heating rate, heat effect or mass loss and "TtaDa" (to be replaced by some useful)
        Dat2<-cbind(format(res[[8]], digits=7),format(res[[7]], digits=7),
                    format(res[[9]], digits=4, scientific=T),format(res[[10]], digits=4, scientific=T))
        Dat<-rbind(c(FN,mess1,mess2,"TtaDa"),Dat2)
        write.table(Dat, file = logfilename, sep=" ", col.names=FALSE, row.names=FALSE, quote = FALSE) 
        return(list(FN,Dat))
      }
      else {
        oldData=read.table(logfilename, header = TRUE, check.names = FALSE)
        Fnlist=character(ncol(oldData)/4)
        for (i in 1:(ncol(oldData)/4)) Fnlist[i]=names(oldData)[4*(i-1)+1]
        newFile=1
        FileInd=match(FN,Fnlist,nomatch=0)
        if (FileInd!=0) {
          FtypeList=character(ncol(oldData)/4)
          for (i in 1:(ncol(oldData)/4)) FtypeList[i]=ifelse(pmatch("HE",names(oldData)[4*(i-1)+3],nomatch=0)==0,"int","dif")
          for (k in 1:length(FileInd)) if (FtypeList[FileInd[k]]==isolate(input$signal)) newFile=0
          }
        if (newFile==1) {
          #write data to file (4 columns) as above
          Dat2<-cbind(format(res[[8]], digits=7),format(res[[7]], digits=7),
                      format(res[[9]], digits=4, scientific=T),format(res[[10]], digits=4, scientific=T))
          Dat2<-rbind(c(FN,mess1,mess2,"TtaDa"),Dat2)
          oldDat<-rbind(names(oldData),oldData)
          Dat<-unname(cbind(oldDat,Dat2))
          write.table(Dat, file = logfilename, sep=" ", col.names=FALSE, row.names=FALSE, quote = FALSE) 
          Fnlist[length(Fnlist)+1]=FN
          return(list(Fnlist,Dat))
        }
        aS=1
        Dat<-cbind(oldData,aS)                                                        #three lines doing rain dance to made Dat suitable for plotting
        Dat$aS=NULL
        return(list(Fnlist,unname(Dat)))
      }
    }
    else return (NULL)
  })
  
  output$downloadData <- downloadHandler( #generates file with kinetic data checked 
    filename = function() {             
      paste("KinData-", Sys.Date(), ".txt", sep="") },
    content = function(file) {
      res=writeKinData()
      KinDat=as.matrix(res[[2]])
      FNlist=res[[1]]
      FLind=match(input$FNchecked,FNlist)
      k=1
      arrInd=integer(length(FLind)*4)
      for (i in 1:length(FLind)) {
        arrInd[k:(k+3)]=(4*(FLind[i]-1)+1):(4*(FLind[i]-1)+4)
        k=k+4 }
      KinDat=KinDat[,arrInd]
      write.table(KinDat[-1,], file, sep=" ", eol = "\r\n", col.names=KinDat[1,], row.names=FALSE, quote = FALSE) 
    })

  output$Plot1 <- renderPlot( height= function(){ ifelse(input$Compact==T, 350, input$BaseHeight)   },
    {  #plot with raw input data (signal vs. temperature/time) and processed (baseline substracted, smoothed) data
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    resOriginal=readData() #temp, t, y
    res=treatData()
    yLabs=ifelse(input$signal=="dif","Differential signal", "Integral signal")
    if (input$tprg=="niso") { 
        x1=resOriginal[[1]]
        XLab=expression('Temperature, °C') }
    else { 
        x1=resOriginal[[2]]
        XLab=expression('Time, s') }
    if (input$signal=="dif") {
      if (input$tprg=="niso") xA=res[[8]]
      else xA=res[[7]]
      yA=res[[11]]
      numb=abs(res[[5]])
      dHshow=round(numb/10^(trunc(log10(numb))-2))*10^(trunc(log10(numb))-2)*ifelse(res[[5]]<0,-1,1)
      dHunits=ifelse(as.numeric(input$dscUnits)==3," ","J/g")
      mess2=paste("Heat effect =",dHshow,dHunits)
    }
    if (input$signal=="int") {
      numb=res[[5]]
      if (input$tprg=="niso") xA=res[[8]]
      else xA=res[[7]]
      if (as.numeric(input$tgaUnits)==4) yA=res[[9]]*100
      if (as.numeric(input$tgaUnits)==3) yA=res[[9]]
      if (as.numeric(input$tgaUnits)<3) yA=res[[15]]-res[[9]]*(res[[15]]-res[[17]])
      dMshow=ifelse(numb>10,round(numb*10)/10,round(numb*100)/100)
      dMunits=ifelse(as.numeric(input$tgaUnits)>2," ","%")
      mess2=paste("Mass loss =",dMshow,dMunits)
      if (as.numeric(input$tgaUnits)<3) {
        points(input$sta,res[[12]], col=28, pch=19, cex=2)                        #denoted right side for TGA normalization
        points(input$end,res[[13]], col=28, pch=19, cex=2) }                      #denoted left side for TGA normalization
      #abline(h=res[[4]], col=8)
    }
    par (fig=c(0,0.8,0,1.0),mar=c(5.1,5.5,4.1,0),las=1)
    plot(x1,resOriginal[[3]], xlab=XLab, ylab=yLabs, cex.lab=input$AxisSize, ylim=c(min(resOriginal[[3]],yA),max(resOriginal[[3]],yA)),
         font.main=input$TitleFont, cex.main=input$TitleSize, type="l", lwd=6, col=rgb(0,0,0,0.3))
    if (input$signal=="dif") abline(h=0, col=1)
    mess1=paste("Heating rate:",res[[6]],"K/min")
    lines(xA,yA, col='red', lwd=1.5)
    mtext(mess2, side = 3, line = 1, adj=0, cex=input$LegendSize)
    mtext(mess1, side = 3, line = 2, adj=0, cex=input$LegendSize)
    #orient=ifelse(input$signal=="dif","topleft","bottomright")
    legend("topright", inset=c(0,-0.24),c("Original data","Modified data"), lwd=c(6,1.5), col=c(rgb(0,0,0,0.3),"red"), 
           bty = "n", cex=input$LegendSize, xpd=T)
  })
  
  output$Plot2 <- renderPlot(height= function(){ ifelse(input$Compact==T, 300, input$BaseHeight) }, 
    { #-------------------------------------------------------------------------draw plot for baseline construction
    if ((is.null(input$file1))||(input$signal=="int")||(input$BLtype=="no")) return(NULL)
    res=treatData()
    dscUnc=res[[12]]
    dPts=length(dscUnc) 
    dMin=min(dscUnc[1],dscUnc[dPts])
    dMax=max(dscUnc[1],dscUnc[dPts])
    yMin=min(res[[12]])
    yMax=max(res[[12]])
    if (input$tprg=="niso") XLab=expression('Temperature, °C') 
    else XLab=expression('Time, s') 
    par (fig=c(0,0.8,0,1.0),mar=c(5.1,5.5,2.1,0),las=1)
    plot(res[[1]],res[[12]],
         main = "Original data zoomed for baseline construction",
         xlab=XLab, ylab="Differential signal", cex.lab=input$AxisSize, font.main=input$TitleFont, cex.main=input$TitleSize,
         ylim=input$YscaleRange)
    lines(res[[1]],res[[13]], col="blue")
    points(input$sta,res[[14]], col=28, pch=19, cex=2)                                    #denoted right side of initial region for baseline construction
    if (input$BLtype=="spln2c") points(input$cntr,res[[16]], col=28, pch=19, cex=2)         #denoted central point for baseline construction
    points(input$end,res[[15]], col=28, pch=19, cex=2)                                    #denoted left side of final region for baseline construction
  })

  output$Plot3 <- renderPlot( height= function(){ ifelse(input$Compact==T, 300, input$BaseHeight)  },
    { #plot with conversion and conversion rate for the current datafile
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    res=treatData()
    par (fig=c(0,0.8,0,1.0),mar=c(4.1,5.5,3.1,0),las=1)
    if (input$tprg=="niso") { 
      x1=res[[8]]
      XLab=expression('Temperature, °C') }
    else { 
      x1=res[[7]]
      XLab=expression('Time, s') }
    plot(x1,res[[9]], 
         main = "Conversion data prepared to load into kinetic project",
         type="l", cex.lab=input$AxisSize, font.main=input$TitleFont, cex.main=input$TitleSize,
         xlab=XLab, ylab=expression('Conversion '*alpha), lwd=2, col='red')
    par(new = TRUE)
    plot(x1,res[[10]], type ="l", axes=FALSE, bty="n", xlab="", ylab="", lwd=2, col='blue')
    axis(side=4, at = pretty(range(res[[10]])))
    mtext(expression('d'*alpha*'d'*italic('t')*', s'^-1), side=4, line=4, cex=input$AxisSize, las=3)
  })
  
  output$Plot4 <- renderPlot( height= function(){ ifelse(input$Compact==T, 350, input$BaseHeight)  },
    { #plot conversion or conversion rates of files loaded in project and selected
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    if (!input$go1) return(NULL)
    if (is.null(input$FNchecked)) return(NULL)
    res=writeKinData()
    KinDat=as.matrix(res[[2]])
    FNlist=res[[1]]
    FLind=match(input$FNchecked,FNlist)
    k=1
    arrInd=integer(length(FLind)*4)
    for (i in 1:length(FLind)) {
      arrInd[k:(k+3)]=(4*(FLind[i]-1)+1):(4*(FLind[i]-1)+4)
      k=k+4 }
    KinDat=KinDat[,arrInd]
    NrunsD=ncol(KinDat)/4
    Tmin=30000
    Tmax=-30000
    dAmin=30000
    dAmax=-30000
    desS=ifelse(input$showA==0,4,3)
    xesS=ifelse(input$showT==1,1,2)
    KinDatNames=KinDat[1,]
    KinDat=KinDat[-1,]
    HRval=character(NrunsD)
    for (i in 1:NrunsD) {
      Tval=as.numeric(unlist(KinDat[,(4*(i-1)+xesS)]))
      Tmin=ifelse(Tval[1]<Tmin,Tval[1],Tmin)
      Tmax=ifelse(Tval[length(Tval)]>Tmax,Tval[length(Tval)],Tmax) 
      dAval=as.numeric(unlist(KinDat[,(4*(i-1)+desS)]))
      minV=min(dAval)
      dAmin=ifelse(minV<dAmin,minV,dAmin)
      maxV=max(dAval)
      dAmax=ifelse(maxV>dAmax,maxV,dAmax) 
      HRval[i]=gsub('HR|Kmin','',KinDatNames[(4*(i-1)+2)])
      }
    if (NrunsD>5) col=rainbow(NrunsD)
    else col=c('red','blue','black','green','grey')
    par (fig=c(0,0.8,0,1.0),mar=c(5.1,5.5,3.1,0),las=1)
    plot(KinDat[,xesS],KinDat[,desS],
         main = "Current kinetic project",
         type="l",lwd=2, cex.lab=input$AxisSize, font.main=input$TitleFont, cex.main=input$TitleSize,
         xlab=ifelse(input$showT==1,'Temperature, °C','Time, s'), 
         ylab=ifelse(input$showA==0,expression('d'*alpha*'d'*italic('t')*', s'^-1),expression('Conversion '*alpha)), 
         xlim=c(Tmin,Tmax), ylim=c(dAmin,dAmax), col=col[1])
    legA=character(NrunsD)
    legA[1]=paste0(HRval[1]," K/min")
    lwdA=numeric(NrunsD)
    lwdA[1]=2
    if (NrunsD>1) {
      for (i in 2:NrunsD) {
        lines(as.numeric(unlist(KinDat[,(4*(i-1)+1)])),as.numeric(unlist(KinDat[,(4*(i-1)+desS)])), col=col[i], lwd=2) 
        legA[i]=paste0(HRval[i]," K/min")
        lwdA[i]=2 } }
    orient=ifelse(input$showA==0,"topleft","bottomright")
    legend(orient, legend=legA, lwd=lwdA, col=col[1:NrunsD],inset = .02, cex=input$LegendSize, bty = "n")
  })

  output$prjData <- renderUI({    #returns list of loaded in project files as group of checkboxes
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    if (!input$go1) return(NULL)
    FNlist=writeKinData()
    checkboxGroupInput("FNchecked", "Loaded to project files:", choices = FNlist[[1]], selected = FNlist[[1]])
  })
  
  output$thirdCol <- renderUI({  #three functions below return selectors for tga/dsc - column number, units, name
    if (input$signal=="dif") return (selectInput("dscCol", label=NULL,#width="75px",
                                         choices = c("1"=1,"2"=2,"3"=3,"4"=4,"5"=5),selected = 3))
    else return(selectInput("tgaCol", label=NULL,#width="75px",
                            choices = c("1"=1,"2"=2,"3"=3,"4"=4,"5"=5),selected = 3))
  })
  output$thirdName <- renderUI({
    if (input$signal=="dif") return ("DSC")
    else return("TGA")
  })  
  output$thirdUnits <- renderUI({
    if (input$signal=="dif") return (selectInput("dscUnits", label=NULL,#width="75px",
                                                 choices = c("W/g"=1, "mW"=2, "1/s"=3),selected = 1))
    else return(selectInput("tgaUnits", label=NULL,#width="75px",
                            choices = c("%"=1, "mg"=2, "0..1"=3, "0..100"=4),selected = 1))
  }) 
    
})
