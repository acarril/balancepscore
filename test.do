clear all
cd "C:/Users/alvar/Dropbox (JPAL LAC)/Dina/procurement/ado/balancepscore"
use ../../data/procurement_test, clear
discard
balancepscore ///
high_direct p1  p2 p3 size_PRE1 size_PRE2 audited dTR1- dTR3 Year2 Year1 Zone1- Zone3 /// 
if (dis_cutoff2>-4 & dis_cutoff2<4), ///
  weight(weight4) pscore(ps_flexmodel41) comsup logit addnamtex(_4) ///
  namgroup(Low share/High share)
