import sys
from Bio.Align import substitution_matrices
from os.path import exists
from tqdm import tqdm

# Compute a highest scoring alignment of 2 profiles using blosum50, with linear gap penalty d and sum of pairs scoring
# Modified Needleman Wunsch

blosum50 = substitution_matrices.load("BLOSUM50")


def score(x, y, i, j, d=8):
    """Scores an alignment between column i in x and column j in y using the blosum50 matrix and sum of pairs scoring"""

    if gapPenalty < 0:
        gapPenalty *= -1  #  Gap penalty needs to be a positive integer

    total = 0

    for k in range(len(x)):
        for l in range(len(y)):
            a, b = x[k][i], y[l][j]
            if (a == '-') ^ (b == '-'):
                total -= d
            elif a != '-' and b != '-':
                total += blosum50[x[k][i], y[l][j]]
    
    return total


def align(x, y, d=8, printMatrices=False):
    """Computes one optimal alignment of profiles x and y"""
    lenX, lenY = len(x[0]), len(y[0])  # Length of a sequence in profiles x and y

    f = [[0 for i in range(lenY+1)] for i in range(lenX+1)]
    t = [['' for i in range(lenY+1)] for i in range(lenX+1)]

    # Compute the weights matrix
    for i in tqdm (range (lenX+1), desc="Calculating Weights Matrix"):
        for j in range(lenY+1):
            if i == 0:
                f[i][j] = -j * d * len(x)
                t[i][j] = 'L'
            elif j == 0:
                f[i][j] = -i * d * len(y)
                t[i][j] = 'T'
            else:
                diag = f[i-1][j-1] + score(x, y, i-1, j-1)
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

    # change to multiple sequences --------------------------
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


def readFile(path, a, b):
    """Reads the text file f containing profiles and puts them in arrays a and b"""

    with open(path) as f:
        numLines = int(f.readline())
        for x in range(numLines):
            a.append(list(f.readline().replace('\n', '')))
        
        numLines = int(f.readline())
        for x in range(numLines):
            b.append(list(f.readline().replace('\n', '')))


if __name__ == '__main__':

    try:
        x = sys.argv[1]
    except:
        raise SystemExit(f'Usage: {sys.argv[0]} <file> \n\t<file> - The path to the file containing alignments')
    
    if not exists(x):
        raise SystemExit(f'Error: {sys.argv[1]} not found.')
    
    if x[-4:] != '.txt':
        raise SystemExit(f'Error: {sys.argv[1]} needs to be a text file.')
    
    a = []
    b = []
    readFile(x, a, b)

    s, x, y = align(a, b, 8, False)

    print("Score: " + str(s) + "\n")
    
    for i in range(len(x)):
        for j in range(len(x[i])):
            print('%s' % x[i][j], end='')
        print("")
    print("")
    for i in range(len(y)):
        for j in range(len(y[i])):
            print('%s' % y[i][j], end='')
        print("")

    