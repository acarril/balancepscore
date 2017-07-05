program define balancepscore3, rclass byable(recall)
	version 11.1
syntax varlist [fweight iweight pweight] [if] [in] , [pscore(string) blockid(string) numblock(int 5) optblock comsup logit latex(string) level(real 5) dir(string) namgroup(string) bdec(int 3)]

tokenize `varlist'
local numvar `: word count `varlist''

/* retrieve the treatment indicator */
local T  `1' 

*** save covariates
macro shift 
local cov "`*'"
local numcov `: word count `cov''


tempvar touse
g `touse'=0

qui replace `touse'=1 `if' `in'

/* if weights are specified, create local variable */
if "`weight'" != "" { 
   tempvar wv
   qui gen double `wv' `exp'
   local w [`weight'=`wv']
   replace `touse'=0 if `wv'==0
}

**** Name group
if `"`namgroup'"' != `""'  {
local pos=strpos("`namgroup'","/")
local G0=substr("`namgroup'",1,`pos'-1)
local G1=substr("`namgroup'",`pos'+1,.)
}
else {
local G0="G0"
local G1="G1"
}

/*******/
/* NEW */
/*******/

***Define PSCORE
if `"`pscore'"' != `""'  { /* BEGINDETAIL */
   confirm new variable `pscore'
} 
else{
tempvar pscore
}

if `"`blockid'"' != `""'  { /* BEGINDETAIL */
   confirm new variable `blockid'
} /* ENDDETAIL */
else{
tempvar blockid
}



qui tab `T'  if `touse'==1

if `"`logit'"' != `""'  { 
   capture drop comsup
   qui logit `varlist'  if `touse'==1       
}
else {
   capture drop comsup
 qui  probit `varlist' [`weight'`exp'] if `touse'==1

}


tempvar epscore

qui predict double `epscore' if `touse'==1

ereturn clear

/* NEW */
/*******/

*capture drop `pscore'

qui gen double `pscore' = `epscore'
label var `pscore' "Estimated propensity score"


/* REGION OF COMMON SUPPORT */
if `"`comsup'"' != `""'  {
   qui sum `pscore' if `T'==1
   tempname mintreat maxtreat
   tempvar COMSUP
   scalar `mintreat'  = r(min)
   scalar `maxtreat'  = r(max)

   qui g `COMSUP'=(`pscore'>=`mintreat' & `pscore'<=`maxtreat')
   qui replace `COMSUP'=. if `touse'==0
   tempvar touse2
   qui gen `touse2'=`touse'
   qui replace `touse'=0 if `COMSUP'==0
}



*** OPTIMAL NUMBER OF BLOCKS:
if  `"`optblock'"' != `""' {
local optblock=0
local numblock=1
while `optblock'!=`numblock' {

local optblock=`optblock'+1
local numblock=`numblock'+1

/*******/

** Generate percentile
local k=0
local j=`optblock'-1
forval i=1/`j' {
local k=`k'+1
local perc`i'=`k'*100/`optblock'
qui centile `pscore' if `touse'==1, c(`perc`i'')
local p`i'=r(c_1)
}


***Generate blockid
local p0=0
local p`optblock'=`optblock'
qui cap gen `blockid'=.
local k=-1
forval i=1/`optblock' {
local k=`k'+1
qui replace `blockid'=`i' if `p`k''<=`pscore' & `pscore'<=`p`i'' & `touse'==1
*** In case of random selection of Stata
local p`k'=`p`k''+0.0000001
}

***  P-VALUE balance global of the covariates
local int 
local int_test
forvalue i=1/`optblock' {
qui tempvar block`i'
qui gen `block`i''=(`blockid'==`i')
local int_test`i'
 foreach x in `cov' {
 qui tempvar block`i'_`x'
 qui gen `block`i'_`x''=`block`i''*`x'
 local int_test`i' `int_test`i'' `block`i'_`x''
 }
 local int_test `int_test' `int_test`i''
}

qui reg `T'  i.`blockid' `int_test' if `touse'
qui test `int_test'
local pvalue: di %12.`bdec'fc  1-F(r(df),r(df_r),r(F))
local sign=`level'/100
if `pvalue'>`sign' {
local numblock= `optblock'
}
cap drop `block`i'_`x''  `block`i''  cap 
}
}
else {
/*******/

** Generate percentile
local k=0
local j=`numblock'-1
forval i=1/`j' {
local k=`k'+1
local perc`i'=`k'*100/`numblock'
qui centile `pscore' if `touse'==1, c(`perc`i'')
local p`i'=r(c_1)
}

***Generate blockid
local p0=0
local p`numblock'=`numblock'
qui cap gen `blockid'=.
local k=-1
forval i=1/`numblock' {
local k=`k'+1
qui replace `blockid'=`i' if `p`k''<=`pscore' & `pscore'<=`p`i'' & `touse'==1
*** In case of random selection of Stata
local p`k'=`p`k''+0.0000001
}
}


*** Table for distribution by block and treatment/group:

di in ye _newline(3) "**************************************************** "
di in ye	     "BLOCK DISTRIBUTION "
di in ye	     "**************************************************** "


forval j=0/1 {
forval i=1/`numblock' {
qui summarize `T' if `T'==`j' & `blockid'==`i'
local m`i'`j'=r(N)
}
}

*** Total by block
forval i=1/`numblock' {
local Total`i'=`m`i'0'+`m`i'1'
}

*** Total by group

local TotalSample0=0
local TotalSample1=0
local TotalSample=0
forval i=1/`numblock' {
local TotalSample0=`TotalSample0'+`m`i'0'
local TotalSample1=`TotalSample1'+`m`i'1'
local TotalSample=`TotalSample'+`Total`i''
}


	
		tempname dblock
		matrix `dblock' = J(`numblock'+1,3,.)
		
	
		forvalue i= 1/`numblock' {
			matrix `dblock'[`i',1] = `m`i'0'	
			matrix `dblock'[`i',2] = `m`i'1'
			matrix `dblock'[`i',3] = `Total`i''				
			local rown "`rown' "Block `i'""
		}
		
		** Add Total to the matrix:
		local k=`numblock'+1
		matrix `dblock'[`k',1] = `TotalSample0'
		matrix `dblock'[`k',2] = `TotalSample1'
		matrix `dblock'[`k',3] = `TotalSample'		

		matrix colnames `dblock' = "`G0'" "`G1'" Total
		matrix rownames `dblock' = `rown' "Total Sample"
		
			local form ", noheader"
			if "`format'" != "" {
				local form "`form' `format'"
			}
			 matrix list `dblock' `form'
		
		if "`matrix'" != "" {
			 matrix `matrix' = `dblock'
		}
	

	 return matrix dblock = `dblock'
		
		
*** LATEX

if "`latex'"=="blockdist"	| "`latex'"=="all" {
*** Row titles

forval i=1/`numblock' {
local lbl_`i' "Block `i'"	
}
local lbl_`k' "Total Sample"
if `"`dir'"' != `""'  {
texdoc init `dir'/blockdist_`T'.tex, replace
}
else {
texdoc init blockdist_`T'.tex, replace
}
tex \hline\hline \\ [-1.5ex]
tex {} & (1) & (2) & (3)  \\
tex [1ex] \\ [-1.5ex]
tex  & `G0' & `G1' & Total \\\\
forvalue i=1/`numblock' {
tex `lbl_`i'' & `m`i'0' & `m`i'0' & `Total`i''  \\
}
tex `lbl_`k'' & `TotalSample0' & `TotalSample1' & `TotalSample' \\
tex [1ex] \hline\hline \\ [-1.5ex]
texdoc close
}	

** Balance within blocks

di in ye _newline(3) "**************************************************** "
di in ye	     "BALANCE WITHIN BLOCKS "
di in ye	     "**************************************************** "


qui ssc install orth_out
foreach i of numlist 1(1)`numblock'{
di in ye "**************************************************** "
di in ye	     "BLOCK `i'"
di in ye	     "**************************************************** "


qui orth_out `cov' if `touse' & `blockid'==`i' , replace by(`T') test count
qui mat m=r(matrix)

qui tabstat `cov' if `touse', stat(sd) save 
qui matrix overall= r(StatTotal)'

**Number of observations
local Block`i'N0=m[`numcov'+1,1] 
local Block`i'N1=m[`numcov'+1,2]

*** mean, std mean difference, pvalues
forvalues j=1/ `numcov' {
local Block`i'_`j'_1: di m[`j',1]
local Block`i'_`j'_2: di m[`j',2]
local Block`i'_`j'_3: di (m[`j',1]-m[`j',2])/overall[`j',1]
local Block`i'_`j'_4: di m[`j',3]
}

*** abs sd mean difference 
local k=0 
forvalue j=1/`numcov'{ 
foreach x  of numlist `Block`i'_`j'_3' {
 local k=abs(`x')+`k' 
 } 
 }
 
 local m=`numcov'+1
 local Block`i'_`m'_3: di `k'/`numcov'
 
  *** F statistics 
  qui reg `varlist' if `blockid'==`i' & `touse' 
  local m=`m'+1
  local Block`i'_`m'_4: di e(F)

  *** p-value
  local m=`m'+1
  local Block`i'_`m'_4: di 1-F(e(df_m),e(df_r),e(F))
  
		tempname bblock`i'
		matrix `bblock`i'' = J(`numcov'+4,4,.)
		
	local j=0                              
	foreach var of varlist `cov' {
	local j=`j'+1  
			matrix `bblock`i''[`j',1] = round(`Block`i'_`j'_1',10^(-`bdec'))	
			matrix `bblock`i''[`j',2] = round(`Block`i'_`j'_2',10^(-`bdec'))	
			matrix `bblock`i''[`j',3] = round(`Block`i'_`j'_3',10^(-`bdec'))
		
			matrix `bblock`i''[`j',4] = round(`Block`i'_`j'_4',10^(-`bdec'))
			local rown2`i' "`rown2`i'' `var'"
	}
	
			matrix `bblock`i''[`numcov'+1,1] = `Block`i'N0'
			matrix `bblock`i''[`numcov'+1,2] = `Block`i'N1'
			local m=`numcov'+1
			matrix `bblock`i''[`numcov'+2,3] = round(`Block`i'_`m'_3',10^(-`bdec'))
			local m=`m'+1		
			matrix `bblock`i''[`numcov'+3,4] = round(`Block`i'_`m'_4',10^(-`bdec'))
			local m=`m'+1			
			matrix `bblock`i''[`numcov'+4,4] = round(`Block`i'_`m'_4',10^(-`bdec'))
	

		matrix colnames `bblock`i'' = "Mean `G0'" "Mean `G1'" "StMeanDiff" p-value	
		matrix rownames `bblock`i'' = `rown2`i'' Observations Abs(StMeanDiff) F-statistic p-value
	
			local form ", noheader"
			if "`format'" != "" {
				local form "`form' `format'"
			}
			matrix list `bblock`i'' `form'
		
		if "`matrix'" != "" {
			matrix `matrix' = `bblock`i''
		}

		return matrix bblock`i' = `bblock`i''
		
*** LATEX if numblock is 1 or >3		
if ("`latex'"=="balblock" | "`latex'"=="all") & (`numblock'==1 | `numblock'>3) {
*** Row titles
local j=0                              
foreach var of varlist `cov' {
local j=`j'+1
local label`j': variable label `var'
local lbl_`j' "`label`j''"	
}

local k=`numcov'+1
local lbl_`k'"Abs(St. mean diff.)"
local m=`k'+1
local lbl_`m' "F-statistic"
local m=`m'+1
local lbl_`m' "P-value"

*** Rounding values:
forval j=1/`m' {
forval l=1/4 {
local Block`i'_`j'_`l': di %12.`bdec'f `Block`i'_`j'_`l''
}
}
if `"`dir'"' != `""'  {
texdoc init `dir'\BalanceBlock`i'_`T'.tex, replace
}
else {
texdoc init BalanceBlock`i'_`T'.tex, replace
}
tex \hline\hline \\ [-1.5ex]
tex {} & (1) & (2) & (3) & (4)  \\
tex [1ex] \\ [-1.5ex]
tex  & `G0' & `G1' & &  \\
tex & (n=`Block`i'N0') & (n=`Block`i'N1') & &   \\
tex [1ex] \\ [-1.5ex]
tex & Mean & Mean & St.Mean Diff. & P-value  \\\\
forvalue j=1/`numcov' {
tex `lbl_`j'' & `Block`i'_`j'_1' & `Block`i'_`j'_2' & `Block`i'_`j'_3' & `Block`i'_`j'_4' \\
}
tex \hline \\
forvalue t=`k'/`m' {
tex `lbl_`j'' &  &  & `Block`i'_`t'_3'  & `Block`i'_`t'_4'   \\
}
tex [1ex] \hline\hline \\ [-1.5ex]
texdoc close
}	
}

*** LATEX if numblock is 2 or 3
if ("`latex'"=="balblock"	| "`latex'"=="all") & (`numblock'==2 | `numblock'==3) {
*** Row titles
local j=0                              
foreach var of varlist `cov' {
local j=`j'+1
local label`j': variable label `var'
local lbl_`j' "`label`j''"	
}

local k=`numcov'+1
local lbl_`k'"Abs(St. mean diff.)"
local m=`k'+1
local lbl_`m' "F-statistic"
local m=`m'+1
local lbl_`m' "P-value"

*** Rounding values:
forval j=1/`m' {
forval l=1/4 {
local Block1_`j'_`l': di %12.`bdec'f `Block1_`j'_`l''
local Block2_`j'_`l': di %12.`bdec'f `Block2_`j'_`l''
local Block3_`j'_`l': di %12.`bdec'f `Block3_`j'_`l''
}
}

if `"`dir'"' != `""'  {
texdoc init `dir'\BalanceBlock_`T'.tex, replace
}
else {
texdoc init BalanceBlock_`T'.tex, replace
}
tex \hline\hline \\ [-1.5ex]
if `numblock'==2 {
tex {} & (1) & (2) & (3) & (4) & (5) & (6) & (7) & (8)  \\
}
else {
tex {} & (1) & (2) & (3) & (4) & (5) & (6) & (7) & (8) & (9) & (10) & (11) & (12) \\
}
tex [1ex] \\ [-1.5ex]
if `numblock'==2 {
tex &\multicolumn{4}{c}{Block 1} &\multicolumn{4}{c}{Block 2} \\\\
}
else {
tex &\multicolumn{4}{c}{Block 1} &\multicolumn{4}{c}{Block 2} &\multicolumn{4}{c}{Block 3} \\\\

}
tex [1ex] \\ [-1.5ex]
if `numblock'==2 {
tex  & `G0' & `G1' & & & `G0' & `G1' & & \\
}
else{
tex & `G0' & `G1' & & & `G0' & `G1' & && `G0' & `G1' \\
}
if `numblock'==2 {
tex & (n=`Block1N0') & (n=`Block1N1') & & & (n=`Block2N0') & (n=`Block2N') & &  \\ 
}
else{
tex & (n=`Block1N0') & (n=`Block1N1') & & & (n=`Block2N0') & (n=`Block2N1') & & & (n=`Block3N0') & (n=`Block3N1') & &  \\ 
}
tex [1ex] \\ [-1.5ex]
if `numblock'==2 {
tex & Mean & Mean & St.Mean Diff. & P-value & Mean & Mean & St.Mean Diff. & P-value  \\\\
}
else {
tex & Mean & Mean & St.Mean Diff. & P-value & Mean & Mean & St.Mean Diff. & P-value & P-value & Mean & Mean & St.Mean Diff. & P-value \\\\
}
forvalue j=1/`numcov' {
if `numblock'==2 {
tex `lbl_`j'' & `Block1_`j'_1' & `Block1_`j'_2' & `Block1`j'3' & `Block1_`j'_4' & `Block2_`j'_1' & `Block2_`j'_2' & `Block2_`j'_3' & `Block2_`j'_4' \\
}
else{
tex `lbl_`j'' & `Block1_`j'_1' & `Block1_`j'_2' & `Block1`j'3' & `Block1_`j'_4' & `Block2_`j'_1' & `Block2_`j'_2' & `Block2_`j'_3' & `Block2_`j'_4' & `Block3_`j'_1' & `Block3_`j'_2' & `Block3_`j'_3' & `Block3_`j'_4'\\
}
}
tex \hline \\
forvalue t=`k'/`m' {
if `numblock'==2 {
tex `lbl_`t'' &  &  & `Block1_`t'_3'  & `Block1_`t'_4' &  &  & `Block2_`t'_3'  & `Block2_`t'_4'  \\
}
else{
tex `lbl_`t'' &  &  & `Block1_`t'_3'  & `Block1_`t'_4' &  &  & `Block2_`t'_3'  & `Block2_`t'_4' &  &  & `Block3_`t'_3'  & `Block3_`t'_4' \\
}
}
tex [1ex] \hline\hline \\ [-1.5ex]
texdoc close
}


di in ye _newline(3) "**************************************************** "
di in ye	     "BALANCE IMPROVE "
di in ye	     "**************************************************** "

*** Original balance

qui orth_out `cov' if `touse2', replace by(`T') test count
qui mat m=r(matrix)


**Observations:
local N0=m[`numcov'+1,1]
local N1=m[`numcov'+1,2]

***** Overall SD *******
qui tabstat `cov' if `touse2', stat(sd) save
matrix overall_org= r(StatTotal)'

forvalues j=1/`numcov'  {
* Mean G0
local m`j'1:  di m[`j',1]
* Mean G1
local m`j'2:  di m[`j',2]
* Std mean diff.
local m`j'3: di (m[`j',1]-m[`j',2])/overall_org[`j',1]
* P-value
local m`j'4:  di m[`j',3]
}

*** abs sd mean difference 
local k=0
forvalue j=1/`numcov' {
foreach x  of numlist `m`j'3' {
local k=abs(`x')+`k'
}
}

** Compute
local l=`numcov'+1
local m`l'3:di `k'/`numcov'

*** F statistics

qui reg `varlist' if `touse2'
local l=`l'+1
local m`l'4: di e(F)

** p-valor
local l=`l'+1
local m`l'4: di 1-F(e(df_m),e(df_r),e(F))

**************************
*** After propensity score
**************************

forvalues j=1/`numcov' {
* Mean G0
local Block`j'_1=0
* Mean G1
local Block`j'_2=0

foreach i of numlist 1(1)`numblock' {
local Block`j'_1: di `Block`j'_1' + (`Block`i'_`j'_1'*`Total`i'')/`TotalSample'	
local Block`j'_2: di `Block`j'_2' + (`Block`i'_`j'_2'*`Total`i'')/`TotalSample'
}

*Std mean difference
local Block`j'_3: di (`Block`j'_1'-`Block`j'_2')/overall[`j',1]
}

*p-value
tempvar Block_T
qui gen `Block_T'=`blockid'*`T'
local j=0
 foreach var in `cov' {
local j=`j'+1
 qui reg `var' i.`blockid' i.`Block_T' if `touse'
 qui testparm i(1/`numblock').`Block_T'
 local Block`j'_4: di  1-F(r(df),r(df_r),r(F))
}

*** abs sd mean difference 
local k=0
forvalue j=1/`numcov' {
foreach x  of numlist `Block`j'_3' {
local k=abs(`x')+`k'
}
}

** Compute
local l=`numcov'+1
local Block`l'_3:di `k'/`numcov'

*** F-STATISTIC and P-VALUE global
local int 
local int_test
forvalue i=1/`numblock' {
 qui tempvar block`i'
 qui gen `block`i''=(`blockid'==`i')
local int_test`i'
 foreach x in `cov' {
 qui tempvar block`i'_`x'
 qui gen `block`i'_`x''=`block`i''*`x'
 local int_test`i' `int_test`i'' `block`i'_`x''
 }
 local int_test `int_test' `int_test`i''
}

 qui reg `T'  i.`blockid' `int_test' if `touse'
 qui test `int_test'
 
local l=`l'+1
local Block`l'_4: di  r(F)
local l=`l'+1
local Block`l'_4: di  1-F(r(df),r(df_r),r(F))


di in ye             "**************************************************** "
di in ye	     "ORIGINAL BALANCE "
di in ye	     "**************************************************** "


		tempname orbal
		matrix `orbal' = J(`numcov'+4,4,.)
		
	local j=0                              
	foreach var of varlist `cov' {
	local j=`j'+1  
			matrix `orbal'[`j',1] = round(`m`j'1',10^(-`bdec'))	
			matrix `orbal'[`j',2] = round(`m`j'2',10^(-`bdec'))
			matrix `orbal'[`j',3] = round(`m`j'3',10^(-`bdec'))
			matrix `orbal'[`j',4] = round(`m`j'4',10^(-`bdec'))
			local rown3 "`rown3' `var'"
	}
	
			matrix `orbal'[`numcov'+1,1] = `N0'
			matrix `orbal'[`numcov'+1,2] = `N1'			
			local l=`numcov'+1
			matrix `orbal'[`numcov'+2,3] = round(`m`l'3',10^(-`bdec'))
			local l=`l'+1		
			matrix `orbal'[`numcov'+3,4] = round(`m`l'4',10^(-`bdec'))
			local l=`l'+1			
			matrix `orbal'[`numcov'+4,4] = round(`m`l'4',10^(-`bdec'))

	
		matrix colnames `orbal' = "Mean `G0'" "Mean `G1'" "StMeanDiff" p-value 
		matrix rownames `orbal' = `rown3' Observations Abs(StMeanDiff) F-statistic p-value
	
			local form ", noheader"
			if "`format'" != "" {
				local form "`form' `format'"
			}
			matrix list `orbal' `form'
		
		if "`matrix'" != "" {
			matrix `matrix' = `orbal'
		}

		return matrix orbal = `orbal'

di in ye             "**************************************************** "
di in ye	     "BALANCE AFTER PROPENSITY SCORE "
di in ye	     "**************************************************** "


		tempname balimp
		matrix `balimp' = J(`numcov'+4,4,.)
		
	local j=0                              
	foreach var of varlist `cov' {
	local j=`j'+1  
			matrix `balimp'[`j',1] = round(`Block`j'_1', 10^(-`bdec'))	
			matrix `balimp'[`j',2] = round(`Block`j'_2', 10^(-`bdec'))	
			matrix `balimp'[`j',3] = round(`Block`j'_3', 10^(-`bdec'))
			matrix `balimp'[`j',4] = round(`Block`j'_4', 10^(-`bdec'))	
			
			local rown4 "`rown4' `var'"
	}
	
			matrix `balimp'[`numcov'+1,1] = `TotalSample0'
			matrix `balimp'[`numcov'+1,2] = `TotalSample1'			
			local l=`numcov'+1
			matrix `balimp'[`numcov'+2,3] = round(`Block`l'_3',10^(-`bdec'))
			local l=`l'+1		
			matrix `balimp'[`numcov'+3,4] = round(`Block`l'_4',10^(-`bdec'))			
			local l=`l'+1			
			matrix `balimp'[`numcov'+4,4] = round(`Block`l'_4',10^(-`bdec'))			


		matrix colnames `balimp' = "Mean `G0'" "Mean `G1'" "StMeanDiff" p-value
		matrix rownames `balimp' = `rown4' Observations Abs(StMeanDiff) F-statistic p-value
	
			local form ", noheader"
			if "`format'" != "" {
				local form "`form' `format'"
			}
			matrix list `balimp' `form'
		
		if "`matrix'" != "" {
			matrix `matrix' = `balimp'
		}

		return matrix balimp = `balimp'
	
*** TEX
if ("`latex'"=="balimp"	| "`latex'"=="all")  {


*** Row titles
local j=0                              
foreach var of varlist `cov' {
local j=`j'+1
local label`j': variable label `var'
local lbl_`j' "`label`j''"	
}

local k=`numcov'+1
local lbl_`k'"Abs(St. mean diff.)"
local m=`k'+1
local lbl_`m' "F-statistic"
local m=`m'+1
local lbl_`m' "P-value"

** Rounding values

forval j=1/`m' {
forval l=1/4 {
local m`j'`l': di %12.`bdec'f  `m`j'`l''
local Block`j'_`l': di %12.`bdec'f `Block`j'_`l''
}
}
if `"`dir'"' != `""'  {
texdoc init `dir'\BalanceImprove_`T'.tex, replace
}
else {
texdoc init BalanceImprove_`T'.tex, replace
}

tex \hline\hline \\ [-1.5ex]
tex {} & (1) & (2) & (3) & (4) & (5) & (6) & (7) & (8)  \\
tex [1ex] \\ [-1.5ex]
tex &\multicolumn{4}{c}{Original balance} &\multicolumn{4}{c}{Balance after propensity score stratification} \\\\
tex [1ex] \\ [-1.5ex]
tex  & `G0' & `G1' & & & `G0' & `G1' & & \\
tex & (n=`N0') & (n=`N1') & & & (n=`TotalSample0') & (n=`TotalSample1') & &  \\ 
tex [1ex] \\ [-1.5ex]
tex & Mean & Mean & St.Mean Diff. & P-value & Mean & Mean & St.Mean Diff. & P-value  \\\\
forvalue j=1/`numcov' {
tex `lbl_`j'' & `m`j'1' & `m`j'2' & `m`j'3' & `m`j'4' & `Block`j'_1' & `Block`j'_2' & `Block`j'_3' & `Block`j'_4' \\
}
tex \hline \\
forvalue t=`k'/`m' {
tex `lbl_`t'' &  &  & `m`t'3'  & `m`t'4' &  &  & `Block`t'_3'  & `Block`t'_4'  \\
}
tex [1ex] \hline\hline \\ [-1.5ex]
texdoc close
}
	
*** RETURN

eret clear 
if  `"`optblock'"' != `""' {

return scalar optblock=`numblock'

}

if `"`comsup'"' != "" {
   qui g comsup = `COMSUP'
   label var comsup "Dummy for obs. in common support"
}
	

end
