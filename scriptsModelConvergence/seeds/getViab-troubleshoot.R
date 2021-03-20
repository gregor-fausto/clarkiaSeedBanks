# script to calculate viabilities

PVinit = if(a==1){
  1
} else if(a==2) {
  PV1[p,r]
} else {
  PV2[p,r]
}

Viab=0; Start=0; NB=0; PVmean=0; PVwts=0;

# …
I = which(Pop %in% PopNow & Round==r & Age==a & !is.na(JanIntact) & !is.na(OctIntact))

# pv=c()
# PVtot=c()

for(b in 1:length(I)){
  Ib = I[b];
  J = which(PopV %in% Pop[Ib] & BagV==Bag[Ib])
  if(!(length(J)==0)){
    NB = NB + 1;
    pv = sum(NumStain[J])/sum(ViabStart[J]);
    GS = sum(GermStart[J]);
    NG = sum(NoGerm[J]);
    printit = 0;
    if(GS>0){
      if(!is.na(pv)){
        PVtot = ( GS + NG*(pv - 1))/GS;
        Viab = Viab + OctIntact[Ib]*PVtot;
      } else {
        printit=0;
        PVtot = (GS - NG)/GS;
        Viab = Viab + OctIntact[Ib]*PVtot;
      }
      if( PVtot < PVinit){
        mult = (PVinit^(2/3))*(PVtot^(1/3));
      } else {
        mult = PVtot^(1/3);
      }
      Start = Start + JanIntact[Ib]*mult;
      PVmean = PVmean + PVtot*GS;
      PVwts = PVwts+GS;
    }
    
    bag = unique(BagV[which(PopV %in% Pop[Ib] & BagV==Bag[Ib])]) 
  } else {
    pv=NA;PVtot=NA;mult=NA;
    bag=NA;
  }
  
  
  mat.list[index,] = data.frame(PopNow,r,a,bag,pv,PVtot,mult,
                                      sum(NumStain[J]),sum(ViabStart[J]),
                                      sum(GermStart[J]),
                                      sum(NoGerm[J]), PVmean, PVwts) %>% dplyr::mutate(PopNow=as.character(PopNow)) 
  index = index+1
}     

if(Start>0){
  SS = Start;
  surv = Viab/Start
} else {
  SS = 0;
  surv=0;
}

pv1.list[index2,] = data.frame(PopNow,r,a,PVmean,PVwts,PVratio=PVmean/PVwts) %>% 
  dplyr::mutate(PopNow=as.character(PopNow)) 
index2 = index2+1

PVmean = PVmean/PVwts; 


pv.list[[p]] = pv

