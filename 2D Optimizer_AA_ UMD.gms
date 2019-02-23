$Title Sensors selection optimization

$Ontext
#############################################
##########2D Optimizer Module################ Copyright 2018, Amin Aria, All rights reserved.
#############################################


$Offtext

*********INPUTs that their change lead to changes in the final data****************
*********INPUTs that their change lead to changes in the final data****************
*********INPUTs that their change lead to changes in the final data****************
*********INPUTs that their change lead to changes in the final data****************
*********INPUTs that their change lead to changes in the final data****************
*********INPUTs that their change lead to changes in the final data****************


Set j  candidate nodes number  /1*11/;
Set i  defects number /1*11/;

Scalar L length of pipeline
       r radius of pipe;

****damagae spatial clustering********damagae spatial clustering
set node11(i) /1*5/;
set node22(i) /6*9/;
set node33(i) /10*11/;
*set node44(i) /18*25/;
*set node55(i) /26*29/;
*set node66(i) /30*38/;
*set node77(i) /39*44/;
*set node88(i) /45*51/;
*set node99(i) /52*60/;
*set node1010(i) /61*67/;


*****Constraints limits*****************************
scalar nHI   human inspection reduction factor
       nPig  Pigging reduction factor
       b     Cost limit
       nf    infinity value
       DLL   lower limit for detection
       DLU   Upper limit for detection
       NUL   Node Usage upper Limit
       PL    Prognosis Lower Limit (PoND exponent)
       PU    Prognosis Upper Limit (PoND exponent)

       Cset1 cluster 1 limit
       Cset2 cluster 2 limit
       Cset3 cluster 3 limit
       Cset4 cluster 1 limit
       Cset5 cluster 2 limit
       Cset6 cluster 3 limit
       Cset7 cluster 1 limit
       Cset8 cluster 2 limit
       Cset9 cluster 3 limit
       Cset10 cluster 1 limit

       Rulim Upper bound on Redundancy
       MdetU Max detection of a damage
       COSw  Cost     weight in utility function
       IVw   IV       weight in utility function
       FRw   Freq     weight in utility function
       COVw  Coverage weight in utility function
       MEw   ME       weight in utility function
       HImax maximum number of damages that can be in the coverage area of HI
       Uw   0-1 normalizer for total utility
;
****pipe dimensions****pipe dimensions****pipe dimensions****pipe dimensions****pipe dimensions
L=50;
r=1;
nf=L;

**HI prognosis constants****HI prognosis constants** **HI prognosis constants** **HI prognosis constants**
nHI=1.5;
nPig=1;


****cost limit****cost limit******cost limit******cost limit******cost limit******cost limit**
b=0.2*card(j);

****Detection limits****Detection limits****Detection limits****Detection limits****Detection limits
DLL=0.5;
DLU=1;
NUL=DLL+0.1;
MdetU=3;

****prognosis limits****prognosis limits****prognosis limits****prognosis limits****prognosis limits
PL=1.5;
PU=12;

****clusters****clusters ****clusters ****clusters ****clusters ****clusters ****clusters ****clusters
*Refer to SHM2Opt paper to know more about clustering bounds
Cset1=2;
Cset2=3;
Cset3=1;


****Redundancy****Redundancy ****Redundancy ****Redundancy ****Redundancy ****Redundancy ****Redundancy
Rulim=1.5;


****Utility Weights****Utility Weights ****Utility Weights ****Utility Weights ****Utility Weights
COSw=6;
IVw=3;
FRw=2;
COVw=6;
MEw=3;

Uw=COSw+IVw+FRw+COVw+MEw;

**********************************
Parameter cov11 coverage of AE 1
          cov12
          cov2
          cov3  pigging coverage;
cov11=0.4;
*cov12=0.8;
cov2=20;
*cov3=100;

HImax=(card(i)*cov2)/l;

Parameter cost11 cost of AE 1
          cost12
          cost2
          cost3  pigging direct cost;
cost11=1;
*cost12=15;
cost2=20;
*cost3=300;

Parameter IV11 Information value of AE 1
          IV12
          IV2
          IV3  pigging Information value;
IV11=1;
*IV12=1;
IV2=50;
*IV3=1;

Parameter FR11 data frequency of AE 1
          FR12
          FR2
          FR3  pigging frequnecy;
FR11=100;
*FR12=1;
FR2=1;
*FR3=1;

**********************************************************************************************
****End of Constraints,bounds and constants****End of Constraints,bounds and constants *******
**********************************************************************************************



********LPoND matrix**************************************************************************
$CALL GDXXRW.EXE Optimizer-Feeder_AA_UMD.xlsx trace=3 Squeeze=N par=LPOND rng=Sheet1!M25 maxDupeErrors=1000000
parameter LPOND(i,j) -logarithm of  missing a defect at i when a sensor is used at j for defect type 1

$GDXIN Optimizer-Feeder_AA_UMD.gdx
$LOAD LPOND
$GDXIN

*****Pexist(DELTA)matrices for AE and HI*****************************************************

$CALL GDXXRW.EXE Optimizer-Feeder_AA_UMD.xlsx trace=3 Squeeze=N par=DeltaAE rng=Sheet1!M52  maxDupeErrors=1000000
parameter DeltaAE(i,j) -log(1-PoD(distance size-type)) for AE type 1

$GDXIN Optimizer-Feeder_AA_UMD.gdx
$LOAD DeltaAE
$GDXIN

$CALL GDXXRW.EXE Optimizer-Feeder_AA_UMD.xlsx trace=3 Squeeze=N par=DeltaHI rng=Sheet1!AJ52  maxDupeErrors=1000000
parameter DeltaHI(i,j) -log(1-PoD(distance size-type)) for AE type 1

$GDXIN Optimizer-Feeder_AA_UMD.gdx
$LOAD DeltaHI
$GDXIN



*****PoD(distance,size-type) for AE *********************************************************

$CALL GDXXRW.EXE Optimizer-Feeder_AA_UMD.xlsx trace=3 Squeeze=N par=LPoDT1 rng=Sheet1!M76  maxDupeErrors=1000000
parameter LPoDT1(i,j) -log(1-PoD(distance size-type)) for AE type 1

$GDXIN Optimizer-Feeder_AA_UMD.gdx
$LOAD LPoDT1
$GDXIN

$CALL GDXXRW.EXE Optimizer-Feeder_AA_UMD.xlsx trace=3 Squeeze=N par=LPoDT2 rng=Sheet1!AJ76  maxDupeErrors=1000000
parameter LPoDT2(i,j) -log(1-PoD(distance size-type)) for AE type 1

$GDXIN Optimizer-Feeder_AA_UMD.gdx
$LOAD LPoDT2
$GDXIN


Parameter MEvalue11(j)   Measurement Error for AE sensor type 1
/
1        0.065
2        0.041
3        0.065
4        0.014
5        0.038
6        0.041
7        0.014
8        0.065
9        0.038
10        0.038
11        0.041
/

Parameter MEvalue2(j)   Measurement Error for HI with Ultra Sonic
/
1        0.341
2        0.286833333
3        0.286833333
4        0.286833333
5        0.286833333
6        0.30075
7        0.273333333
8        0.273333333
9        0.273333333
10       0.3248
11       0.3248
/

Parameter maxME ME normalizer;
maxME=0.34;


****END OF INPUTS****END OF INPUTS****END OF INPUTS****END OF INPUTS****END OF INPUTS****END OF INPUTS****END OF INPUTS
****END OF INPUTS****END OF INPUTS****END OF INPUTS****END OF INPUTS****END OF INPUTS****END OF INPUTS****END OF INPUTS
****END OF INPUTS****END OF INPUTS****END OF INPUTS****END OF INPUTS****END OF INPUTS****END OF INPUTS****END OF INPUTS



binary variables s0(j)   no sensor at node j indicator
                 s11(j)  AE 1 at node j indicator
*                 s12(j)  AE 2 at node j indicator
                 s2(j)   HI with USonic at node j indicator;
*                 s3(j)   Pigging at node j indicator;




variable cov(j)  coverage at node j
         cost(j) cost of tool at node j
         IV(j)   inoformation value at node j
         FR(j)   sensing frequnecy at node j
         ME(j)   measurement error for sensor at node j and the closest damage;


        variables obj            likelihood of detection and cost
                  totaldetection total number of detected defects
*                  totalcost      total cost of sensors
                  Redun          Redundancy of detection (hit-miss modeling)
                  avgdet         average detection
                  misslog(i)     logarithm of prognosis failure for defect i
                  nPOD           avelog(probability of missing defects)
*                  PODdist(i)     probability of detection (distribution model)
                  LPoDT(i)       -log(1- PoD(distance)*PoD(size)*PoD(hit-miss)) for damage i
                  utility(j)     a combination of cost and other data gathering tools attributes
                  totdety        total usage of nodes
                  totalut        total utiliy
                  aveut          average utility per used node;


binary variables delta(i,j)     detection indicator (damage i by tool at node j)
                 detyuse(j)     jth sensor usage
                 detect(i,j)    detection of defect i by sensor j
                 gendetect(i)   general detection of damage i;




Equations

           sensind(j)            only one sensor at each node
           coveq(j)              defining cov(j)
           costeq(j)             defining cost(j)
           IVeq(j)               defining IV(j)
           FReq(j)               defining Fr(j)
           MEeq(j)               defining ME(j)

           totcost               cost should be smaller than a predetermined value (exponential relation)
*           costlim               cost smaller than a limit

           utdef(j)              defenition of utility function (to be optimized)
           totutility            total utility equation
           aveutility            average utility per used node


           Gdetect(i)            definition of general detect(i)
           GdetectUp(i)          upperbound for gendetect(i)
           maxdet(i)             maximum detection of each defect

           usage(j)              each sensor much be used at least once

           senser(j)              A sensor shouldn't be used if it is not detecting

           nodeg1                node group 1 *set partitioning problem
           nodeg2                node group 2
           nodeg3                node group 3
*           nodeg4                node group 4
*           nodeg5                node group 5
*           nodeg6                node group 6
*           nodeg7                node group 4
*           nodeg8                node group 5
*           nodeg9                node group 6
*           nodeg10                node group 6

           usetot                rate of usage should be greate than some lower bound
           equality(i,j)         if damage i can be detected by sensor at node j
           Redundef              definition of redundancy of detection
           detavg                average detection of defects

           misslog1(i)           definition of misslog(i)
           missdef1(i)           lower bound for missing a defect
           missdef2(i)           upper bound for missing a defect
           missdefect            probability of missing existing defects should be smaller an upper bound
           PoDTeq(i)             PoD(hit-miss) * PoD (dist) equation

*          PODdisteq(i)          probability of detection (distribution model) definition
           redlim                limit on redundancy

           avgliml               lower limit for average detection
           avglimu               upper limit for average detection
           objective             likelihood of detection and cost as objective ;



sensind(j)..     s0(j)+s11(j)+s2(j) =e= 1;

 coveq(j)..      cov(j)  =e= s11(j)*cov11+s2(j)*cov2;
 costeq(j)..     cost(j) =e= s11(j)*cost11+s2(j)*cost2;
 IVeq(j)..       IV(j)   =e= s11(j)*iv11+iv2*(sum(i,deltaHI(i,j))/HImax)*s2(j);
 FReq(j)..       FR(j)   =e= s11(j)*fr11+s2(j)*fr2;
 MEeq(j)..       ME(j)   =e= s11(j)*MEvalue11(j)+s2(j)*MEvalue2(j);

totcost..              sum(j,cost(j)) =l= cost2*b;
*costlim..              totalcost =l= 1;

utdef(j)..      utility(j) =e= (-COSw*(cost(j)/cost2)+COVw*(cov(j)/cov2)+IVw*(IV(j)/IV2)+FRw*(FR(j)/FR11)-(MEw*ME(j)+(3-MEw)*MaxME)/MaxME+COSw+MEw)/Uw;
***************************************************************
totutility..    totalut    =e= sum(j, utility(j));
aveutility..    aveut * card(j) =e= totalut;
***************************************************************

*********************************************************************************************
equality(i,j)..        delta(i,j)      =e= s11(j)*deltaAE(i,j)+s2(j)*deltaHI(i,j);
*********************************************************************************************


maxdet(i)..            sum(j,delta(i,j))          =l= MdetU;
Gdetect(i)..           sum(j,delta(i,j))/1000     =l= gendetect(i);
GdetectUp(i)..         sum(j,delta(i,j))          =g= gendetect(i);

usage(j)..             sum(i,delta(i,j))=l= card(i)* detyuse(j) ;

senser(j)..            1-s0(j) =l= detyuse(j);

nodeg1..               sum(node11,gendetect(node11))=g= Cset1;
nodeg2..               sum(node22,gendetect(node22))=g= Cset2;
nodeg3..               sum(node33,gendetect(node33))=g= Cset3;
*nodeg4..               sum(node44,gendetect(node44))=l= Cset4;
*nodeg5..               sum(node55,gendetect(node55))=l= Cset5;
*nodeg6..               sum(node66,gendetect(node66))=l= Cset6;
*nodeg7..               sum(node77,gendetect(node77))=l= Cset7;
*nodeg8..               sum(node88,gendetect(node88))=l= Cset8;
*nodeg9..               sum(node99,gendetect(node99))=l= Cset9;
*nodeg10..              sum(node1010,gendetect(node1010))=l= Cset10;


usetot..               sum(j,detyuse(j)) =l= NUL*card(j);

detavg..               avgdet =e= sum(i,gendetect(i))/card(i);
avgliml..              avgdet =g= DLL;
avglimu..              avgdet =l= DLU;

**********************************************************************************************************
misslog1(i)..          misslog(i) =e= sum(j,detyuse(j)*LPOND(i,j))+ sum(j,DeltaHI(i,j)*s2(j)*nHI);
**********************************************************************************************************
missdef1(i)..          misslog(i) =g= PL;
missdef2(i)..          misslog(i) =l= PU;


PoDTeq(i)..            LPoDT(i)=e=  sum(j,(DeltaAE(i,j)*s11(j)*LPoDT1(i,j)+DeltaHI(i,j)*s2(j)*LPODT2(i,j)));
redlim..               redun =l= Rulim;

Redundef..             sum((i,j),delta(i,j))=e= Redun*card(i);
missdefect..           sum(i,misslog(i)+LPoDT(i)) =e= nPOD*card(i);

objective..            obj   =e= 0.5*totalut/(card(i)*NUL)+0.5*nPOD/(PU);


Model  sensorarr /all/;
option reslim = 5000;
*sensorarr.optfile = 1;
*Option MINLP=BARON;

Solve sensorarr using MIP maximizing obj;
*execute_unload 'sensors.gdx';
option profile=3;
Display avgdet.l,Redun.l,aveut.l,nPoD.l,s11.l,s2.l,ME.l,misslog.l,detyuse.l;
*sensorarr.MODELSTAT, sensorarr.SOLVESTAT,delta.l,nodelta.l,usetot.l,detyuse.l;




