#############################################
##########Damage Generator Module############ Copyright 2018, Amin Aria, All rights reserved.
#############################################



#############INPUTs#################INPUTs#################INPUTs#################INPUTs#################INPUTs####
#############INPUTs#################INPUTs#################INPUTs#################INPUTs#################INPUTs####

##damage denisty per meter
rate= 12/50;

##Length and radius of pipeline
l=50;
r=1;

##Number of samples to be produced 
#Refer to the journal paper for more info on number of samples required and Wilks method)
nwilks=93;

##Number of size-type catagories 
#Considering the size distriubtion, size-type classes are produced. Refer to journal paper for more info.
ns=4 

###Energy cosntraint
#sensor energy consumption is an exponential function of its distance to closest sensor.
#Refer to the journal paper for more info. (econ = 18 meters)
econ=18

###Inference constraint. 
#Considering the similarity matrix that should be defined for the structure, results of inference are not reliable if 
#the distance between monitored damage and infered damage exceedes a limit, icon. 
icon=20

###Cluster distance constraint
#Cluster distance limit is defined considering energy and inference constraints.
#Refer to the journal paper for more info. 
clcon=min(econ,icon)

############End of INPUTs#####################End of INPUTs#####################End of INPUTs######################
############End of INPUTs#####################End of INPUTs#####################End of INPUTs######################





#Determinatin of mesh size with probability of having more than 1 damage smaller than 0.01 
a=qpois(0.99,rate)
meshsize=1/a;
1-ppois(1,rate*meshsize)

#A damage is placed in a mesh if the result of a binomial test is success.
p=dpois(1,rate*meshsize)

#Number of meshes(strips) along a pipeline segment under study is determined here.
#Refer to the journal paper for more information about meshing. 
Meshnum=floor(l/meshsize);


dexist=matrix(0,Meshnum,nwilks)
  for (j in 1:nwilks)
    for (i in 1:Meshnum) dexist[i,j]=rbinom(1,1,p)

######Longitudinal location of damages############################################
#Poisson Distribution is used here for longitudinal placement of damages
x=seq(0,l,l/(Meshnum-1));
x=x+runif(1,0,meshsize);
xexist=apply(dexist,2,function(y) y*t(x))
xexist2=xexist[dexist>0]

#################################################################################

#####Circumferential location of damages
#It is assumed that longitudinal and circumferential distriubtion of damages 
#along a pipeline segmentare independent. Refer to the journal paper for more info. 

####################Circumferential density######################################

yyc=rnorm(Meshnum,0,pi*r/3)
yyc[yyc<0]=2*pi*r+yyc[yyc<0]
yyc[yyc>2*pi*r]=yyc[yyc>2*pi*r]-2*pi*r
yexist=apply(dexist,2,function(y) y*t(yyc))
yexist2=yexist[yexist>0]


###########################circumferential location of nodes####################
#Nodes are placed considering the final placement of damages in a pipeline sample (realization)
#Refer to the journal paper for more info. 
ynode=yexist2+runif(length(yexist2),-0.5*r,0.5*r)
ynode[ynode<0]=2*pi*r+ynode[ynode<0]
ynode[ynode>2*pi*r]=ynode[ynode>2*pi*r]-2*pi*r


##########################Defining damage type-size class#######################

detind=matrix(0,nrow=length(ynode),ncol=ns)

#NOTE: next line should be changes considering number of classes defined in the INPUTs
colnames(detind)=c('c#1','c#2','c#3','c#4')

indice=sample(1:ns,length(ynode),replace = TRUE)

for ( i in 1: length(ynode)) detind[i,indice[i]]=1


##########################Cannot-link damages (for clustering)##################
#Distance between consecutive damages is claculated to determin Cannot-link damages
#clsutering constraint(clcon). Constraint K-means will be used to cluster damages 
#of each sample considering cannot-link damages. Refer to the journal paper for more info.

##Distance between consecutive damages
dist2=matrix(0,length(xexist2),1)
clmat=matrix(0,length(xexist2),1)
colnames(clmat)=c('cannot-link')
for (i in 1:(length(xexist2)-1))
{

  dist2[i,1]=sqrt((xexist2[i]-xexist2[i+1])^2+min(abs((yexist2[i]-ynode[i+1])),6.28*r-abs((yexist2[i]-ynode[i+1])))^2)
   if (dist2[i,1] > clcon) {clmat[i,1]<-1 }
}

outloc=cbind(xexist2,yexist2,xexist2,ynode,detind,clmat)

####################################################################################################################
####################################################################################################################



##################OUTPUT##################OUTPUT##################OUTPUT##################OUTPUT##################
##################OUTPUT##################OUTPUT##################OUTPUT##################OUTPUT##################


#Change the directory considering your preferences.
write.table(outloc,file.path(getwd(),'wilks.txt'),sep = "  ")



#The output is a table including 'nwilks' samples with columns:
#"xexist2": Longitudinal location of damage
#"yexist2": Circumferential location of damage
#"xexist2": Longitudinal location of node which is equal to that of corresponding damage
#"ynode":   Circumferential location of Node
#"c#1", "c#2","c#3", "c#4": Binary variables that indicate the type-size class of a damage
#"cannot-link": a binary indicator that shows if the corresponding node cannot be connected to the next node

#NOTE: A new sample starts whenever "xexist2" values drop (damages of a sample are sorted based on their longitudinal location)


