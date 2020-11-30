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
    
    
  }
}     

if(Start>0){
  SS = Start;
  surv = Viab/Start
} else {
  SS = 0;
  surv=0;
}

PVmean = PVmean/PVwts; 