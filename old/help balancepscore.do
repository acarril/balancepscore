********************************************************************************
************************** BALANCEPSCORE COMMAND *******************************
********************************************************************************


*************************
****** Help of the commands:

/*

Syntax

			balancepscore2 subgroup [covariates] [if] [in] [, pscore(string)
			blockid(string) numblock(int) optblock comsup logit latex(string)
			dir(string) namgroup(group1/group2) bdec(int)]
		


Options 	Description

pscore      generates a new variable with propensity score estimation
    
blockid		generates a new variable with the block for each observation 
			
numblock    specifies the number of blocks (manually). The default value is 5.

optblock	the number of blocks is choosen optimally

comsup 		sample is restricted to the region of common support. 

logit	    pscore is calculated using logit regression. The defaults is a
			probit.
		
latex 		exports a latex table. "all" if you need export each table in latex; 
			"blockdist" if you only need export the table on the block distribution;
			"balblock" if you only need export the table of the balance by block;
			"balimp" if you only need export the table showing the balance improve. 
		
dir 		specifies the directory where each latex is placed
		

Options 	Description

namgroup	specifies the name of each group. Example: control/treatment. This is
			to be shown in the header of each table. The default is g0/g1
		
bdec		specifies the number of decimal places reported for all estimates. 
			The default value for bdec is 3.
			
		
Stored results

    balancepscore stores the following in r():

    Scalars   
      r(optblock)    optimal number of block

    Matrix   
      r(balimp)      balance after propensity score matching
      r(orbal)       original balance
      r(bblock"j")   balance in block j
      r(dblock)      distribution of observations by block.		
		
*/


*****************
** EXAMPLES:

clear all
set maxvar 30000

cd "C:\Users\jpaluser\Dropbox (JPAL LAC)\Chile Compra Project\Stata resources"

/* For these examples I use the "Comparing Inference Approaches for RD Designs:
A Reexamination of the effect of Head Start on Child Mortality" */

use "balancepscore\Cattaneo2016\headstart.dta", clear

*** Generate running variable

gl R povrate60
gl c = 59.1984
gen double R = $R - $c

*** Use restricted sample from the paper:
local BW 13.544

*** generate subgroup variable
centile mort_age59_related_preHS, c(50)
gen high_mortalityPRE=(mort_age59_related_preHS>=r(c_1))

*** Keep subgroup variable and covariates:
keep high_mortalityPRE census1960_pctsch534 census1960_pop534 census1960_pcturban census1960_pctblack R

** Label for formatting in latex:

** Subgroup dummy
label variable high_mortalityPRE  		`"Mortality rate in t-1 > p(50)"'

rename census1960_pop534    pop534
rename census1960_pctsch534 pctsch534
rename census1960_pctblack  pctblack
rename census1960_pcturban  pcturban

** log of population
replace pop534=log(pop534)

** Covariates:
label variable pop534	 	`"Log of Population between 5-34 year"'
label variable pctsch534  	`"% attending school between 5-34 year"'
label variable pctblack		`"% black population"'
label variable pcturban  	`"% urban population"' 

*****************************
**** CASE 1: Manual choise of the number of blocks:

** log file

log using "C:\Users\jpaluser\Dropbox (JPAL LAC)\Chile Compra Project\Stata resources\balancepscore\Example Head Start_log.log", replace

balancepscore high_mortalityPRE pop534 pctsch534 pctblack pcturban if R>=-`BW' & R<=`BW', comsup numblock(5) namgroup(Low Mortality/High Mortality)

log close 

** If you want export tables to latex: 

cd "C:\Users\jpaluser\Dropbox (JPAL LAC)\Chile Compra Project\Stata resources\balancepscore\tex"
balancepscore high_mortalityPRE pop534 pctsch534 pctblack pcturban if R>=-`BW' & R<=`BW', comsup numblock(5) namgroup(Low Mortality/High Mortality) latex(all)
syop
*****************************
**** CASE 2: Optimal number of blocks:

balancepscore high_mortalityPRE pop534 pctsch534 pctblack pcturban if R>=-`BW' & R<=`BW', comsup optblock
di r(optblock) // 12 blocks

