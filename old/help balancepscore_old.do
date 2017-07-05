********************************************************************************
************************** BALANCEPSCORE COMMAND *******************************
********************************************************************************



/*

We have 3 commands. The only thing that differs is in the choice of the optimal 
number of blocks.:


1) balancepscore1: The optimal number of blocks is choosen according to the global
balance of the covariates. This is, the number of blocks with the lowest F-statistic 

2) balancepscore2: The optimal number of blocks is choosen according to the overall 
standardized mean difference. This is, the number of blocks with the lowest std.
mean difference. 

3) balancepscore3: The optimal number of blocks is choosed according to the 
significance level. This is, the first number of blocks in which its p-value is 
higher than a certain significance level.

*/


*************************
****** Help of the commands:

/*

Syntax

		balancepscore1 subgroup [covariates] [if] [in] [, pscore(string)
		blockid(string) numblock(int) optblock comsup logit latex(string)
                dir(string) namgroup(group1/group2) bdec(int)]
		
		balancepscore2 subgroup [covariates] [if] [in] [, pscore(string)
		blockid(string) numblock(int) optblock comsup logit latex(string)
                dir(string) namgroup(group1/group2) bdec(int)]
		
		balancepscore3 subgroup [covariates] [if] [in] [, pscore(string)
		blockid(string) numblock(int) optblock comsup logit latex(string)
                dir(string) namgroup(group1/group2) bdec(int) level(real)


Options 	Description

pscore          generates a new variable with propensity score estimation
    
blockid		generates a new variable with the block for each observation 
			
numblock        specifies the number of blocks (manually). The default value is 5.

optblock	the number of blocks is choosen optimally

comsup 		sample is restricted to the region of common support. 

logit	        pscore is calculated using logit regression. The defaults is a
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


Additionally for the balancepscore3 is added as option:

level 		specifies the significance level to be used in the optimal number of 
		blocks. The default value is 5%.
		
		
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

/* For these examples I use the "Lalonde Experimental Data (Dehejia-Wahba
    Sample)" */

use nsw, clear

*****************************
**** CASE 1: Manual choise of the number of blocks:

*** The 3 commands generate the same result:

balancepscore1 treat age education nodegree, comsup numblock(3) 
dyop
** If you want export tables to latex: 

** Label for formatting in latex:

label variable age `"Age"'
label variable education `"Years of schooling"'
label variable nodegree `"Dummy(no degree)"'

balancepscore1 treat age education nodegree, comsup numblock(3) latex(all)


*****************************
**** CASE 2: Optimal number of blocks:


***  OPTION 1:
balancepscore1 treat age education nodegree, comsup optblock
di r(optblock) // 3 blocks

***  OPTION 2:
balancepscore2 treat age education nodegree, comsup optblock
di r(optblock) // 13 blocks

***  OPTION 3:
balancepscore3 treat age education nodegree, comsup optblock level(10)
di r(optblock) // 2 blocks

