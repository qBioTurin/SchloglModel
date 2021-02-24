

Error_Distance<-function(Trace1,Trace2,K=100)
{
  
  N=table(Trace1$Time)[1]
  M=table(Trace2$Time)[1]

  max<- max(Trace1$X1,Trace2$X1)
  min<-min(Trace1$X1,Trace2$X1)
  L=(max - min )
  br = seq(min,max,by=L/K)

  h_1<-aggregate(Trace1$X1,by = list(Time = Trace1$Time), FUN= function(SSASol)  {
    freq_SSA   = hist(SSASol, breaks=br, include.lowest=TRUE, plot=FALSE)
    h_SSA<-freq_SSA$counts*1/N
    return(h_SSA)
    })
  h_2<-aggregate(Trace2$X1,by = list(Time = Trace2$Time), FUN= function(SSASol)  {
    freq_SSA   = hist(SSASol, breaks=br, include.lowest=TRUE, plot=FALSE)
    h_SSA<-freq_SSA$counts*1/N
    return(h_SSA)
  })
  
  DErr <- apply(abs(h_1[,-1] - h_2[,-1]),1,sum)
  return(DErr)
}


ErrD<-Error_Distance(Gentrace,trace)


ggplot(SSAdensity , aes(x= x , y= h_SSA )) + geom_line()
