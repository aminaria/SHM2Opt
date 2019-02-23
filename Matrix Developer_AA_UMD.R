#############################################
##########Matrix Developer Module############ Copyright 2018, Amin Aria, All rights reserved.
#############################################

!!!###NOTE: The Minimum Spanning Tree sub-module of the clustering module is included here. 


#####Sample Assingment##########Sample Assingment#####Sample Assingment#####Sample Assingment#####Sample Assingment
#####Sample Assingment##########Sample Assingment#####Sample Assingment#####Sample Assingment#####Sample Assingment  

#outloc is the output of Damage Generator module.You can pass the input in a matrix format similar to the format of 
#Damage generator output

first=1  #first cell of the sample under study
last=11   #Last cell of the sample under study
x=outloc[first:last,1]
t=outloc[first:last,2]
y=outloc[first:last,3]
ty=outloc[first:last,4]



#############INPUTs#################INPUTs#################INPUTs#################INPUTs#################INPUTs####
#############INPUTs#################INPUTs#################INPUTs#################INPUTs#################INPUTs####

#radius of pipeline
r=1;
#number of damages
n=length(x) 
#number of classes
ns=4 
#number of sensors types
nsensor=2 

#####Minimum Spanning Tree library : used for clustering damages of each sample
library(optrees)

##Distance between damages and nodes
dist=matrix(0,nrow=length(x),ncol=length(y))

##Probability of missing a defect
POND=matrix(0,nrow=length(x),ncol=length(y))

#Damage Missing relation: exp(deltax/a +b ), boundry conditions: exp(deltax=0)=0.1, exp(deltax=20 meters)= 1 >>>
#!!!###NoTE: you can change this relation considering the geometry of pipeline and pattern recognition resutls 
  ##of the pipeline. 
  
a1=c(8.3,8,8,7)
b1=c(-2.7,-3,-3.5,-2.7)
LPONDint=matrix(,nrow=0,ncol=length(x))

###parameters of PoD(size) for Acoustic Emission
#(Rabiei and Modarres, 2013) when PoD(depth=1.5mm)=1
aAE=3.07 


###parameters of PoD(size) for HI with Ultrasonic tool
#(Chatterjee and Modarres, 2013) 
aUS=7.32
bUS=0.57

#Considering the size distribution in (Ossai et al., 2015)
#Size distribution location (mean) parameter
sm=-1.12

#Size distribution scale (variance) parameter
#ns is the number of classes
sv=0.78
quan=seq(0,1,1/ns)
yquan=qlnorm(quan,-1.12,0.78)

#Maximum size:
yquan[length(yquan)]=1.5

#number of sample points in each class
nsam=1000;

#PoD(distance) for AE Based on (Pollock, 2007)
#Go to Line 167, PoD as a function of distance, to make changes

##Measurement Error for Acoutsic Emission based on (Rabeie and Modarres, 2013)
MEAE=rnorm(nsam*ns,0,0.2/20)

##ME for HI with Ultra sonic based on (Chatterjee and Modarres, 2013) ##########################
#GO TO LINE 215 to change this Formula

##Human Inspection Coverage radius
#For more info, refer to the journal paper.
Coverage=20





############End of INPUTs#####################End of INPUTs#####################End of INPUTs######################
############End of INPUTs#####################End of INPUTs#####################End of INPUTs######################
############End of INPUTs#####################End of INPUTs#####################End of INPUTs######################
############End of INPUTs#####################End of INPUTs#####################End of INPUTs######################



###############MAIN CODE##############MAIN CODE##########MAIN CODE############MAIN CODE#########MAIN CODE##########
###############MAIN CODE##############MAIN CODE##########MAIN CODE############MAIN CODE#########MAIN CODE##########
###############MAIN CODE##############MAIN CODE##########MAIN CODE############MAIN CODE#########MAIN CODE##########
###############MAIN CODE##############MAIN CODE##########MAIN CODE############MAIN CODE#########MAIN CODE##########


###########defect type-size indicator###########defect type-size indicator###########defect type-size indicator####

type0=outloc[first:last,]
typeind=outloc[first:last,5:(4+ns)]

###########calculating LPOND###########calculating LPOND###########calculating LPOND###########calculating LPOND###

for (k in 1:length(a1))
{
  
  a=a1[k]
  b=b1[k]

POND=matrix(0, nrow= length(x),ncol=length(y))

for (i in 1:length(x))
{
  for (j in 1:length(y))
  {dist[i,j]=sqrt((x[i]-y[j])^2+min(abs((t[i]-ty[j])),6.28*r-abs((t[i]-ty[j])))^2)
  POND[i,j]= exp((dist[i,j]/a)+b)
  if (POND[i,j]>1) {POND[i,j]<- 1}
  }
}

LPONDint=rbind(LPONDint,-log(POND))

}

#final LPOND
LPOND= matrix(0,length(x),length(y))

for (i in 1:length(x)){
  LPOND[i,]= LPONDint[i+10*(which(typeind[i,]>0,arr.ind = TRUE)-1),]
  
}

#>>>LPOND output is in the output section

##################Sample generation from classes' truncated distributions################## 


#Mean PoD for each class for AE sensor using Exponential relation in 

pAE=matrix(0,ncol=ns,nrow=1) # 1-exp(-3.07*xplot)

#Mean PoD for each class for HI with Ultrasonic using PoD2 relation in 
#(Chatterjee and Modarres, 2013) with a(threshold)=0
pUS=matrix(0,ncol=ns,nrow=1)

#To calculate POD for each class
xplot=NULL
yy=NULL


for (i in 1:length(yquan-1))
{  
  xtrunc2=qlnorm((runif(nsam)*(plnorm(yquan[i+1],sm,sv)-plnorm(yquan[i],sm,sv))+plnorm(yquan[i],sm,sv)),sm,sv)
  pAE[1,i]=mean(1-exp(-aAE*xtrunc2[1:nsam]))
  pUS[1,i]=mean(1-((1+exp(-aUS*bUS))/(1+exp(aUS*(xtrunc2[1:nsam]-bUS)))))
  xplot=union(xplot,xtrunc2[1:nsam])
}

yy=dlnorm(xplot,sm,sv)
plot.new()
plot(xplot,yy,xlab = 'pit depth (mm)', ylab='probability density',cex.axis=1.5,cex.lab=1.5)

#class lines
for (i in 2:length(yquan)){ abline(v=yquan[i],col='gray',lty=2,lwd=2) }

xLHS<-xplot

###################POD of AE as a function of distance ( PoD(coverage=0.4)=0.00)##########################
PODdist=1-15.625*(dist**3)
PODdist[PODdist<0]=0
PODdist[PODdist>1]=1


#deltaAE is used for Pexist indicator in optimizer
deltaAE=PODdist
deltaAE[deltaAE>0]=1

#deltaHI is used for Pexist indicator in optimizer

deltaHI=dist
deltaHI[deltaHI>Coverage]=0
deltaHI[deltaHI>0]=1


#PoD(size-type) for AE sensor
PODsity=matrix(0,nsensor,ns)
PODsity[1,]=pAE
PODsize1=t(apply(typeind,1,function(typeind) PODsity[1,] * typeind))
PODsize1=(apply(PODsize1,1,sum))

###POD(size-type) for HI with USonic using PoD3 in (Chatterjee and Modarres, 2013)

PODsity[2,]=pUS
PODsize2=t(apply(typeind,1,function(typeind) PODsity[2,] * typeind))
PODsize2=(apply(PODsize2,1,sum))

##################PoD(distance,size-type)#####################################

PODT1=matrix(0.01,length(x),length(y))
PODT1=PODdist*matrix(rep(PODsize1,length(y)),ncol=length(y))


####negative log to be used in the GAMS (for Acoustic Emission)
LPODT1=-log(1-PODT1)

####negative log to be used in the GAMS (for HI with Ultra sonic)
PODT2=matrix(0.01,length(x),length(y))
PODT2=1*matrix(rep(PODsize2,length(y)),ncol=length(y))
LPODT2=-log(1-PODT2)


##############################################ME Intermediate Matrices####################################
MEinter=matrix(0,nsensor,ns)

##Measurement Error of AE sensors
measx=xLHS+MEAE
ME1=matrix(0,1,n)

for (i in 1:length(yquan)-1)
{
  xclass=measx[ xLHS > yquan[i] & xLHS < yquan[i+1]]
  MEinter[1,i]=1-sum(xclass > yquan[i] & xclass < yquan[i+1])/length(xclass)
}

ME1=t(apply(typeind,1,function(typeind) MEinter[1,] * typeind))
ME1=(apply(ME1,1,sum))


####>>>>>HI with Ultrasonic ME
####>>>>>HI with Ultrasonic ME
MEUS=-0.1*(xLHS)+0.01/1+rnorm(length(nsam*ns),0,0.028)
measx=xLHS+MEUS


for (i in 1:length(yquan)-1)
{
  xclass=measx[ xLHS > yquan[i] & xLHS < yquan[i+1]]
  MEinter[2,i]=1-sum(xclass > yquan[i] & xclass < yquan[i+1])/length(xclass)
}
ME2int=t(apply(typeind,1,function(typeind) MEinter[2,] * typeind))
ME2int=(apply(ME2int,1,sum))

####!!!!NOTE: ME2 for each node (HI) is the average ME value for all detected damages by the HI at that node

ME2=matrix(0,1,n)

for (i in 1:length(ME2int))
{ for (j in 1:length(ME2int))
  {
     if (dist[i,j]<Coverage) ME2[i]=ME2[i]+ME2int[j]
  }
  ME2[i]=ME2[i]/sum(dist[i,]<Coverage)
}

##########################
#MEinter[1,]
#MEinter[2,]


#########################Number of clusters#############################
#Using Minimum Spanning Tree-Kruskal Algorithm

edge=c()
for (i in 1:dim(dist)[1])
{
  for (j in 1:i)
    edge=rbind(edge,c(i,j,dist[i,j]))
  
}
clus0= getMinimumSpanningTree(nodes=1:nrow(dist),arcs=edge,algorithm = 'Kruskal')
#clus=msTreeKruskal(nodes=1:nrow(dist),arcs=edge)
clus=as.data.frame(clus0$tree.arcs)
clus=clus[order(clus$weight),]

#############OUTPUTS#####################OUTPUTS#####################OUTPUTS#####################OUTPUTS########
#############OUTPUTS#####################OUTPUTS#####################OUTPUTS#####################OUTPUTS########
#############OUTPUTS#####################OUTPUTS#####################OUTPUTS#####################OUTPUTS########

#LPOND
write.table(LPOND,file.path(getwd(),'output/LPOND.txt'),sep = "  ")

#ME
write.table(ME1,file.path(getwd(),'output/ME1.txt'),sep = "  ")
write.table(t(ME2),file.path(getwd(),'output/ME2.txt'),sep = "  ")

#Delta
write.table(deltaAE,file.path(getwd(),'output/DeltaAE.txt'),sep = "  ")
write.table(deltaHI,file.path(getwd(),'output/DeltaHI.txt'),sep = "  ")

#LPODT
write.table(LPODT1,file.path(getwd(),'output/LPODT1.txt'),sep = "  ")
write.table(LPODT2,file.path(getwd(),'output/LPODT2.txt'),sep = "  ")

View(clus)

####################################################################################################################

