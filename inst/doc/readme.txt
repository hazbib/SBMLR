
After installing XML from http://www.omegahat.org/RSXML, 
an error message from library(XML)can be resolved 
by copying the *.dll files of the XML package 
(e.g. in C:\Program Files\R\rw2000\library\XML\libs) 
to the C:\windows directory (where Windows can find them). 

version changes
1.00 	initial release (11/05/2004)

1.01 	initial release trivially fixed to make SBML code pass through 
	the SBML.org validator by switching the order of prods and mods 

1.10 	substantive upgrade of fderiv to handle steady state and transient 
	microarray data perturbations (11/27/2004)

1.12 	Added the BMC Cancer 2004 example directory and removed similar 
	preliminary scripts from the demo directory of 1.10 (1/04/2005). 
     	Manual.doc was shrunken down to avoid redundancy with publication-based documentation

1.15 	Upgraded manual.doc to include a new fderiv description since the initial 
	release publication documentation was outdated.(1/17/2005) 

1.16 	Minor further upgrade of manual.doc to include figure of functions and objects (2/1/2005). 
