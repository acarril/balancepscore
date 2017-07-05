********************************************************************************
************************** BALANCEPSCORE COMMAND *******************************
********************************************************************************



*************************
****** Help of the commands:


Syntax
		
		balancepscore subgroup [covariates] [if] [in] [, pscore(string)
		blockid(string) numblock(int) optblock comsup logit latex(string)
                dir(string) namgroup(group1/group2) bdec(int)]
		

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
		
