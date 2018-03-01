# Main script to run CrossPlan for systematically planning genetic cross experiments
# Author: Aditya Pratapa
# Date Created: 28 June 2017
# Example run: python master-script.py --mutantInfoFile=inputs/example-infoFile.txt

import sys
sys.path.insert(0, 'src/')
from optparse import OptionParser
from CrossPlan import *

def main(args):
    usage = '''
    master-script.py [options]
    Required Argument:
        mutantInfoFile - 
            A tab-delimited file with one mutant per line. Each line must contain a mutant id, 
            mutant name, genes mutated (ids) in that mutant, mutant viability, whether the mutant
            is source or target. An example file can be found under inputs/example-infoFile.txt
    '''
    parser = OptionParser(usage=usage)


    parser.add_option('', '--mutantInfoFile', dest="mutantInfoFile", type='string', metavar='STR',
                      help='File containing all the mutant information such as viability, sources/targets and mutated genes . Required.')

    parser.add_option('', '--numBatches', dest="numBatches", type='int',
                      help='Number of batches of experiments to be planned. Default=2', default=2)

    parser.add_option('', '--maxCrossesPB', dest="maxCrossesPB", type='int',
                      help='Max number of crosses/mutants made per batch. Default=2', default=2)

    parser.add_option('', '--maxNumMutations', dest="maxNumMutations", type='int',
                      help='Max number of mutations in target genes. Default=4', default=4)

    parser.add_option('', '--outPrefix', dest="outPrefix", type='string', metavar='STR',
                      help='Name to prepend to output files. Default=outputs/out.',default='outputs/out')

    parser.add_option('', '--timeLimit', dest="timeLimit", type='int',
                      help='Optimizer time limit value in seconds. Default=False', default=3600)

    # If no/wrong options are provided, print help
    if len(sys.argv) == 1:
        parser.print_help()
        exit(-1)

    (opts, args) = parser.parse_args()
    print("Running CrossPlan...\n")
    CrossPlan(opts)


if __name__ == '__main__':
    main(sys.argv)
