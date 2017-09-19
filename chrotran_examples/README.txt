To run the example problems, execute the shell script ./multirun.sh
Python scripts are included, which aid in visualizing simulation results

EXAMPLE_1:

This considers a 2x2 matrix of cases, all considering remediation of Cr(VI) in a heterogeneous aquifer:

CT1di.in -- YES direct reaction YES inhibitor (i.e., alcohol)
CT1dx.in -- YES direct reaction NO  inhibitor
CT1xi.in -- NO  direct reaction YES inhibitor
CT1xx.in -- NO  direct reaction NO  inhibitor

EXAMPLE_2:

This illustrates permeability reduction due to bio-fouling and the use of biocide for unclogging purposes. 
A simulation of constant-head injection into a homogeneous aquifer is shown. 
Injection fluid is amended initially with the biostimulant acetate 400 days.
The acetate amendment is then replaced with the biocide dithionite for the remainder of the simulation.
