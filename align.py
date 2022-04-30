# Calculates an initial multiple alignment for sequences by building a profile using needleman-wunsch, then
# applies simulated annealing combined with iterative refinement to find an improved solution

import sys
import random
import math

from os.path import exists
from Bio.Align import substitution_matrices
from tqdm import tqdm

blosum50 = substitution_matrices.load("BLOSUM50")

def scoreAlignment(alignment, gapPenalty=8):
    """Scores a multiple alignment using sum of pairs scoring of all columns"""

    if gapPenalty < 0:
        gapPenalty *= -1  #  Gap penalty needs to be a positive integer

    total = 0
    
    for i in range(len(alignment[0])):
        for j in range(len(alignment)-1):
            for k in range(j+1, len(alignment)):
                if (alignment[j][i] == '-') ^ (alignment[k][i] == '-'):
                    total -= gapPenalty
                elif alignment[j][i] != '-' and alignment[k][i] != '-':
                    total += blosum50[alignment[j][i], alignment[k][i]]
    
    return total


def scoreColumns(x, y, i, j, gapPenalty=8):
    """Scores an alignment between column i in x and column j in y using the blosum50 matrix and sum of pairs scoring"""

    if gapPenalty < 0:
        gapPenalty *= -1  #  Gap penalty needs to be a positive integer

    total = 0

    for k in range(len(x)):
        for l in range(len(y)):
            a, b = x[k][i], y[l][j]
            if (a == '-') ^ (b == '-'):
                total -= gapPenalty
            elif a != '-' and b != '-':
                total += blosum50[x[k][i], y[l][j]]
    
    return total


def alignProfiles(x, y, d=8, printMatrices=False):
    """Computes one optimal alignment of profiles x and y"""
    lenX, lenY = len(x[0]), len(y[0])  # Length of a sequence in profiles x and y

    f = [[0 for i in range(lenY+1)] for i in range(lenX+1)]
    t = [['' for i in range(lenY+1)] for i in range(lenX+1)]

    # Compute the weights matrix
    #for i in tqdm (range (lenX+1), desc="Calculating Weights Matrix"):
    for i in range(lenX + 1):
        for j in range(lenY+1):
            if i == 0:
                f[i][j] = -j * d * len(x)
                t[i][j] = 'L'
            elif j == 0:
                f[i][j] = -i * d * len(y)
                t[i][j] = 'T'
            else:
                diag = f[i-1][j-1] + scoreColumns(x, y, i-1, j-1)
                top = f[i-1][j] - d * len(y)
                left = f[i][j-1] - d * len(x)

                if diag >= top and diag >= left:
                    t[i][j] = 'D'
                    f[i][j] = diag
                elif top >= diag and top >= left:
                    t[i][j] = 'T'
                    f[i][j] = top
                else:
                    t[i][j] = 'L'
                    f[i][j] = left
    
    if printMatrices:
        # print the scores matrix
        for i in range(len(f)):
            for j in range(len(f[i])):
                print('% 4d,' % f[i][j], end='')
            print("")
        print("\n")

        # print the traceback matrix
        for i in range(len(t)):
            for j in range(len(t[i])):
                print('% 4s,' % t[i][j], end='')
            print("")
    
    # Traceback

    xp = [[] for i in range(len(x))]
    yp = [[] for i in range(len(y))]

    i, j = lenX, lenY
    while i + j > 0:
        if t[i][j] == 'D':
            for k in range(len(x)):
                xp[k].insert(0, x[k][i-1])
            for k in range(len(y)):
                yp[k].insert(0, y[k][j-1])
            i-=1
            j-=1
        elif t[i][j] == 'T':
            for k in range(len(x)):
                xp[k].insert(0, x[k][i-1])
            for k in range(len(y)):
                yp[k].insert(0, '-')
            i -= 1
        else:
            for k in range(len(x)):
                xp[k].insert(0, '-')
            for k in range(len(y)):
                yp[k].insert(0, y[k][j-1])
            j -= 1

    s = f[lenX][lenY]
    return (s, xp, yp)


def iterate(alignment, t=10, m=0.95, y=0, d=8):
    """Improves a multiple alignments using simulated annealing and iterative refinement. t is the initial temperature, m is the amount to
    multiple t by to reduce it (must be < 1, eg 0.95), y is the proportion of sequences that will be re-aligned each cycle (if y=0, only 1 
    sequence will be realigned), and d is the gap penalty. Some 'good' values are t=10 and m=0.95, y=0 and d=8"""

    if m >= 1:
        m = 0.95
    
    # get the number of sequences to be realigned each cycle
    if y == 0:
        numRealignedSequences = 1
    else:
        numRealignedSequences = math.ceil(len(alignment) * y)
        if numRealignedSequences == len(alignment):
            numRealignedSequences -= 1

    #print("Number of sequences to be realigned: " + str(numRealignedSequences))

    counter = 0  # How many cycles where the solution hasn't changed

    # Calculations for the loading bar
    n = 0  # How many times the cycle would run uninterrupted
    while t > 1:
        n += 1
        t = t * m
    
    for i in tqdm (range (n), desc="Iterating"):
    #while t > 1 and counter < 7:
        if counter >= 7:  # If the solution hasn't changed
            break

        #print("t: " + str(t))

        prevProfile = list(alignment)  # retain the original profile and score in case the score doesn't improve
        prevScore = scoreAlignment(prevProfile, d)

        sequences = []
        for i in range(numRealignedSequences):
            index = random.randrange(len(alignment))
            sequence = [e for e in alignment.pop(index) if e != '-']  # remove a random sequence without the gaps
            sequences.append(sequence)
        
        # Remove any columns that only have gaps in them, created by removing sequences
        nullColumnIndexes = []
        for i in range(len(alignment[0])):
            columnHasSymbol = False
            for j in range(len(alignment)):
                if alignment[j][i] != '-':
                    columnHasSymbol = True
                    break
            if not columnHasSymbol:
                nullColumnIndexes.insert(0, i)
        for index in nullColumnIndexes:
            for i in range(len(alignment)):
                alignment[i].pop(index)   
        
        for sequence in sequences:
            s, alignment, y = alignProfiles(alignment, [sequence], d, False) # change align
            alignment.append(y[0])

        newScore = scoreAlignment(alignment, d)

        if newScore == prevScore:
            counter += 1
        else:
            counter = 0

        if newScore <= prevScore:
            p = math.exp((newScore - prevScore) / t)
            if random.random() <= p:
                prevProfile = alignment
                prevScore = newScore
        else:
            prevProfile = alignment
            prevScore = newScore
        t = t * m
    
    s = scoreAlignment(alignment, d)
    return (s, alignment)


def printAlignment(alignment):
    """Prints an alignment in a pretty way"""

    for i in range(len(alignment)):
        for j in range(len(alignment[i])):
            print('%s' % alignment[i][j], end='')
        print("")


def readFile(path, a):
    """Reads the text file given by path containing sequences"""

    with open(path) as f:
        numLines = int(f.readline())
        for x in range(numLines):
            a.append(list(f.readline().replace('\n', '')))


if __name__ == '__main__':

    hasAlignment = False
    t, m, y = 10, 0.95, 0

    try:
        i = 0
        if sys.argv[1] == '-y':
            hasAlignment = True
            i = 1
        path = sys.argv[1 + i]

        if len(sys.argv) > 3:
            t = int(sys.argv[2 + i])
            m = float(sys.argv[3 + i])
            y = float(sys.argv[4 + i])
    except:
        raise SystemExit(f'Usage: {sys.argv[0]} [-y] <file> <t> <m> <y> \n-y - Set this flag if the text file already contains an alignment. If the text file contains un-aligned sequences, leave this out\n\n<file> - The path to the file containing alignments\n\n<t> - The initial temperature for simulated annealing, default is 10\n\n<m> - The multiplier to reduce t by, default is 0.95\n\n<y> - The proportion of sequences that will be re-aligned every cycle, defaults to one sequence (Anything over 32 sequences and it is better to set this to 0)')
    
    if not exists(path):
        raise SystemExit(f'Error: {path} not found.')
    
    if path[-4:] != '.txt':
        raise SystemExit(f'Error: {sys.argv[1]} needs to be a text file.')
    
    sequences = []
    readFile(path, sequences)
    alignment = []

    if hasAlignment:
        alignment = sequences
    else:
        # Create an initial alignment using Needleman-Wunsch
        alignment = [sequences.pop(0)]
        for i in tqdm (range (len(sequences)), desc="Aligning"):
            s, alignment, x = alignProfiles(alignment, [sequences[i]], 8, False)
            alignment.append(x[0])

    # Try and improve the alignment using simulated annealing and iterative refinement
    score, alignment = iterate(alignment, t, m, y, 8)

    printAlignment(alignment)
    print("Score: " + str(score))

    
