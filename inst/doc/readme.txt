
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

1.20    Major upgrade. Model object of class SBML defined and summary and == methods defined for it. 
	Read and writes to files are now to and from this SBML model object. 
        The R model definition file has been maintained for model editing, but it has been streamlined.
        For example, functions definitions now only need the bottom line string expression. Conversions 
        to and from SBML and SBMLR are no longer lossy in appeal since expressions are used instead of 
        string with too many parentheses. Many new functions defined to make script writing clear. (5/4/2005)
        
1.21    Cleaning. The model field names rxns, prods, mods, reacts => reactions, products, modifiers and reactants.
	Better names for functions by starting them with action verb, e.g. readSBML, readSBMLR, saveSBML, saveSBMLR.  
	getIncidenceMatrix was incoporated into getModelInfo. (5/5/2005)

1.22    Simulate fixed to have vector rather than matrix passed for mod=1 with "Control" column now automatic for t<0. Simulate
	was also modified to take an initial state vector override of the one in the model. (5/28/2005)
 