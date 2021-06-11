import numpy as np

# Initialize matrix with zeros which don't have pairings
def InitializeMatrix(N):
    matrix = np.empty((N,N))
    matrix[:] = np.NAN
    np.fill_diagonal(matrix, 0)
    for i in range(1, N):
        matrix[i][i-1] = 0
    return matrix

def InitializeTracebackMatrix(N):
    matrix = np.empty((N,N), object)
    matrix[:] = np.NAN
    np.fill_diagonal(matrix, 0)
    for i in range(1, N):
        matrix[i][i-1] = 0
    return matrix

# Get maximum value of K from Bifurcation (merging two substructures)
def Get_Max_K_Value(i, j, matrix):
    maxK = 0
    maxKIndex = 0
    kList = [k for k in range(i,j) if k > i and k < j]
    for k in kList:
        maxK = max(maxK, matrix[i][k] + matrix[k+1][j])
        if maxK == matrix[i][k] + matrix[k+1][j]:
            maxKIndex = k
    return maxK, maxKIndex

# Check if matrix[i] and matrix[j] are paired or not
def Pair_check(tup):
    if tup in [('A', 'U'), ('U', 'A'), ('C', 'G'), ('G', 'C')]:
        return 1
    return 0

# Mapping each character in the input sequence to a unique identifier
# to access them when filling the matrix
def MapInput(sequence):
    rnaDict = {}
    for i in range(len(sequence)):
        rnaDict[i] = sequence[i]
    return rnaDict

def FillMatrix(matrix, mappedRnaDict, tracebackMatrix, N):
    j = 1
    tup = ()
    for k in range (j, N):
        for i in range(0, N - k):
            tup = (mappedRnaDict[i], mappedRnaDict[i+j])
            down = matrix[i+1][i+j]
            left = matrix[i][(i+j)-1]
            diagonal = matrix[i+1][(i+j)-1] + Pair_check(tup)
            bifurcation, maxKIndex = Get_Max_K_Value(i, i+j, matrix)
            
            matrix[i][i+j] = max(down, left, 
                  diagonal, bifurcation)
            
            traceList = []
            if matrix[i][i+j] == down:
                traceList = [i+1, i+j, 'down', -1]
                tracebackMatrix[i][i+j] = traceList
                
            elif matrix[i][i+j] == left:
                traceList = [i, (i+j)-1, 'left', -1]
                tracebackMatrix[i][i+j] = traceList
            
            elif matrix[i][i+j] == diagonal:
                traceList = [i+1, (i+j)-1, 'diagonal', -1]
                tracebackMatrix[i][i+j] = traceList
            
            else:
                traceList = [i, i+j, 'bifurcation', maxKIndex]
                tracebackMatrix[i][i+j] = traceList
                
        j += 1
 
def TraceBack(tracebackMatrix, dotList, i ,j):
    if i == j or i > j: 
        return dotList
    
    if tracebackMatrix[i][j][2] == 'bifurcation':
        TraceBack(tracebackMatrix,  dotList, tracebackMatrix[i][j][0], tracebackMatrix[i][j][3])
        TraceBack(tracebackMatrix, dotList, tracebackMatrix[i][j][3]+1, tracebackMatrix[i][j][1])
        
    elif tracebackMatrix[i][j][2] == 'diagonal':
        dotList[min(i,j)] = '('
        dotList[max(i,j)] = ')'
        TraceBack(tracebackMatrix, dotList ,tracebackMatrix[i][j][0], tracebackMatrix[i][j][1])
    
    else:
        TraceBack(tracebackMatrix , dotList, tracebackMatrix[i][j][0], tracebackMatrix[i][j][1])
        
    return dotList

# Input:
# GGGAAAUCC No Bifurcation
# CGGACCCAGACUUUC Bifurcation
# UAACGUACUGGAGUA Bifurcation
# GGAAUUAGUUAACC Bifurcation
if __name__ == "__main__":
    rna = input("please input the RNA Sequence: ")
    N = len(rna)
    matrix = InitializeMatrix(N)
    tracebackMatrix = InitializeTracebackMatrix(N)
    mappedRnaDict = MapInput(rna)
    FillMatrix(matrix, mappedRnaDict, tracebackMatrix, N)
    dotList = ['.'] * N
    finalList = TraceBack(tracebackMatrix, dotList, 0, N - 1)
    print(''.join(finalList))
