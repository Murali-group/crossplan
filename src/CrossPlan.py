

import os.path
from itertools import chain, combinations
import cplex
import numpy as np
import re

def all_subsets(S):
    '''
    Returns all subsets of a given set S
    '''
    return chain(*map(lambda x: combinations(S, x), range(0, len(S)+1)))       


def CrossPlan(opts):
    """
    :param opts: Options from parsed options.
    :Output: List of crosses to be made in outPrefix-plan.txt
    """
    
    numBatches = int(opts.numBatches)
    maxCrossesPB = int(opts.maxCrossesPB)
    ilpFile = str(opts.outPrefix + '-' + str(numBatches) + '-' + str(maxCrossesPB) + '.lp')
    CrossPlanFile = str(opts.outPrefix + '-' + str(numBatches) + '-' + str(maxCrossesPB) + '-plan.txt')
    #Initialize mutant id lists
    allMutants = []
    inviableMutants = []
    viableMutants = []
    sourceMutants = []
    targetMutants = []

    #Map mutant ids to their names from mutant information file
    mutantNameDict = {}
    
    #Mutant composition (mutated genes) dictionary
    mutatedGenes = {}


    #Read mutant information file
    infile = open(opts.mutantInfoFile, 'r')
    
    for line in infile:
        items = [x.strip() for x in line.rstrip().split('\t')]
        # Skip empty lines or those beginning with '#' comments
        if line == '':
            continue
        if line[0] == '#':
            continue
        allMutants.append(items[0])
        mutantNameDict[items[0]] = items[1]
        # For each mutant m, the mutatedGenes dictionary is the set of genes mutated in m1, i.e., G(m)
        mutatedGenes[items[0]] = set(items[2].rstrip().split(','))
        if items[3] == 'viable':
            viableMutants.append(items[0])
        else:
            inviableMutants.append(items[0])
        if items[4] == 'source':
            sourceMutants.append(items[0])
        elif items[4] == 'target':
            targetMutants.append(items[0])
    print("Identified %d mutants, of which %d are viable and %d are inviable.\nThere are %d source mutants and %d target mutants."
             %(len(allMutants), len(viableMutants), len(inviableMutants), len(sourceMutants), len(targetMutants)))

    # geneticCrossDictionary is a dictionary for genetic cross graph.
    # The keys in this dictionary correspond to a genetic cross.
    # The values in this dictionary are correspond to all mutants that are a result of this genetic cross.
    # A genetic cross is represented by mutantID1_mutantID2 
    geneticCrossDictionary = {}   
    # The two loops run for all pairs of viable mutants
    for mutantID1 in viableMutants:
        for mutantID2 in viableMutants: 
            # However we only need nC2 combinations of mutants, hence mutantID1<mutantID2
            # Two mutants m1 and m2 are crossed only if |G(m1) intersection G(m2)| == 0
            # Two mutants m1 and m2 are crossed only if |G(m1) union G(m2)| has at most maxNumMutations
            if mutantID1<mutantID2 and len(mutatedGenes[mutantID1].intersection(mutatedGenes[mutantID2])) == 0 \
                        and len(mutatedGenes[mutantID1].union(mutatedGenes[mutantID2])) <= opts.maxNumMutations :
                # mutantPowerset: Power set of G(m1) union G(m2), i.e., mutants that are produced as a result of the genetic cross between m1 and m2
                # mutantPowerset_ids: List of all mutant ids corresponding to the genes that are produced as a result of the genetic cross
                mutantPowerset = all_subsets(mutatedGenes[mutantID1].union(mutatedGenes[mutantID2]))
                mutantPowerset_ids = []
                for mutant in mutantPowerset:
                    mutant_id = None
                    for key, val in mutatedGenes.items(): #Map from set of mutated genes to mutant ids
                        if val == set(mutant):
                            mutant_id = key
                            break
                    if mutant_id != None:
                        mutantPowerset_ids.append(mutant_id) 
                    geneticCrossDictionary[str(mutantID1)+'_'+str(mutantID2)] = mutantPowerset_ids


    print '\nWriting ILP...'
    out = open(ilpFile, 'w')
    # Objective funtion goes here.
    # Our goal is maximize the number of target mutants
    out.write('Maximize\n')
    out.write(' obj: ')
    for target in targetMutants:
        out.write(' + a' + target)
    out.write('\n')

    # Constraints begin here. A counter "count" keeps track of the constraint number.
    out.write('Subject To \n')
    count = 1

    # Source mutants constraints
    # Source mutants are marked as available in batch 0.
    # This is achieved by setting the binary variable for source mutants in batch 0 to 1.
    # Non-source mutants are marked as un-available in batch 0.
    # This is achieved by setting the binary variable for non-source mutants in batch 0 to 0.
    for mutant in allMutants:
        if mutant in sourceMutants:
            out.write(' c' + str(count) + ': ' + 'm' + mutant + 'b' + str(0) + '= 1\n')
            count += 1
        else:
            out.write(' c' + str(count) + ': ' + 'm' + mutant + 'b' + str(0) + '= 0\n')
            count += 1

    # Batch size constraints
    # Max crosses in a batch are set to maxCrossesPB   
    for j in range(1, numBatches + 1):
        out.write(' c' + str(count) + ': ')
        for cross in geneticCrossDictionary.keys():
            out.write('+ x' + cross + 'b' + str(j))
        out.write(' <= ' + str(maxCrossesPB) + '\n')
        count += 1

    # Cross input constraints
    for cross in geneticCrossDictionary.keys():
        mutantsCrossed = cross.rstrip().split('_')
        for j in range(1, numBatches + 1):
            out.write(' c' + str(count) + ': + x' + cross + 'b' + str(j))
            for i in range(j):
                out.write('- m' + mutantsCrossed[0] + 'b' + str(i))
            out.write(' <= 0 \n')
            count += 1
    
            out.write(' c' + str(count) + ': + x' + cross + 'b' + str(j))
            for i in range(j):
                out.write('- m' + mutantsCrossed[1] + 'b' + str(i))
            out.write(' <= 0 \n')
            count += 1

    # Mutant input constraints
    for mutant in allMutants:
        for j in range(1, numBatches + 1):
            out.write(' c' + str(count) + ': m'+ mutant + 'b' + str(j))
            for cross, mutantList in geneticCrossDictionary.items():
                if mutant in mutantList:                   
                    out.write(' - x'+ cross + 'b' + str(j))
            out.write(' <= 0\n')
            count += 1

    # Making mutant constraints
    for target in targetMutants:
        out.write(' c' + str(count) + ': -a' + target)
        for k in range(1, numBatches + 1):
            out.write(' + ' + 'm'+ target + 'b' + str(k))
        out.write(' >= 0 \n')
        count += 1
   
    # Specify that all variables in this .lp file are Binary variables
    out.write('Binary\n')
    
    for target in targetMutants:
        out.write(' a' + target + '\n')

    for mutant in allMutants:
        for j in range(numBatches + 1):
            out.write(' m' + mutant + 'b' + str(j) + '\n')

    for cross in geneticCrossDictionary.keys():
        for j in range(1, numBatches + 1):
            out.write(' x' + cross + 'b' + str(j) + '\n')

    out.write('End')
    # Ensure you close the file to avoid "missing max or min" errors.
    out.close()
    print("ILP written to %s\n" %(ilpFile))
    # Create cplex object
    cp = cplex.Cplex()

    # Read the ilpFile
    cp.read(ilpFile)
    start = cp.get_time()
    cp.parameters.timelimit.set(opts.timeLimit)
    # Solves the CrossPlan ILP problem
    cp.solve()
    end = cp.get_time()
    # Total time taken to solve the ILP
    elapsedtime = end - start

    # Obtain variable names, the order is not necessarily the same as it is in the ilp file.
    ilp_vars = cp.variables.get_names()
    num_vars = cp.variables.get_num()

    # Obtain values for each variable. Batch 0 mutants will be set to 1 by contraints.
    ilp_var_vals = cp.solution.get_values()
    objVal = cp.solution.get_objective_value()

    print("%d out of %d target mutants can be made in %d batches with %d crosses per batch" %(objVal, len(targetMutants), numBatches, maxCrossesPB))
    # List of ids of mutants made. Batch 0 mutants will be set to 1 by contraints.
    mutants_made = []
    # List of crosses made.
    crosses_made = []
    # List of target mutants made.
    target_mutants_made = []

    for i in range(num_vars):
        if int(ilp_var_vals[i]) == 1:
            if ilp_vars[i].find(str('m')) == 0:
                ids = re.findall(r'[0-9]+', ilp_vars[i])
                mutants_made.append(int(ids[0]))
            elif ilp_vars[i].find(str('x')) == 0:
                crosses_made.append(ilp_vars[i])
            else:
                target_mutants_made.append(ilp_vars[i])
    mutants_made=list(set(mutants_made))
    # Print final solution
    outFile = open(CrossPlanFile, 'w')
    for j in range(1,numBatches+1):
        outFile.write("Batch" + str(j) +" :\n")
        for cross in crosses_made:
            cross = re.findall(r'[0-9]+', cross) #extract numbers from cross string
            if int(cross[2]) == j:
                outFile.write("Cross mutants " + mutantNameDict[cross[0]]+ " and " + mutantNameDict[cross[1]] + " to obtain mutants: ")
                outMutantNames = []
                for mutant in geneticCrossDictionary[cross[0]+"_"+cross[1]]:
                    outMutantNames.append(mutantNameDict[mutant])
                outFile.write(', '.join(outMutantNames) +'\n')
    print("CrossPlan output written to %s\n" %(CrossPlanFile))
    return

