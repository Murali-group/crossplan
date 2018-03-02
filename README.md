### CrossPlan
The CrossPlan source code for systematically planning genetic cross experiments
#### Overview 
 CrossPlan is a novel methodology for systematically planning genetic crosses to make a set of target mutants from a set of source mutants. We base our approach on a generic experimental workflow used in performing genetic crosses in budding yeast. CrossPlan uses an integer-linear-program (ILP) formulation to maximize the number of target mutants that we can make under certain experimental constraints. Specifically, CrossPlan computes a sequence of genetic cross experiments organized into batches such that we can perform the crosses in each batch in parallel. CrossPlan takes as input a source set *S* of mutants that are available in the lab, a set *T* of target mutants whose phenotypes we are interested in characterizing experimentally, and the number *k* of batches (which reflects the experimental budget) with at most *s* crosses per batch. The plan computed by CrossPlan maximizes the number of target mutants that can be made from the source set in *k* batches.
 
### Dependencies 
Python packages:
  * <a href= https://www.ibm.com/support/knowledgecenter/SSSA5P_12.7.0/ilog.odms.cplex.help/CPLEX/GettingStarted/topics/set_up/Python_setup.html> CPLEX </a>
#### Requirements
Mutant Information file: A tab-delimited file with one mutant per line. Each line must contain a mutant id, mutant name, genes mutated (ids) in that mutant, mutant viability, whether the mutant is source or target. An example file can be found under inputs/example-infoFile.txt

### Usage: 
```
python master-script.py --mutantInfoFile=inputs/example-infoFile.txt
```
For help on other options use:
```
python master-script.py -h
```
### License

GNU GPLv3.0

### How to cite CrossPlan

We will be very glad to hear from you if you use CrossPlan in your work. If you publish a paper that uses CrossPlan, please cite:
1. <a href="https://doi.org/10.1093/bioinformatics/bty072"> CrossPlan: Systematic Planning of Genetic Crosses to Validate Mathematical Models</a>. Aditya Pratapa, Neil Adames, Pavel Kraikivski , Nicholas Franzese, John J Tyson, Jean Peccoud, TM Murali. *Bioinformatics*, 2018. 
