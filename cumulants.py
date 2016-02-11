# Library

from sympy import symbols, zeros, simplify, cancel, ratsimp,\
                    exp, log, diff, sqrt
from numpy import array, copy
from copy import deepcopy
import time

#######################################
# Helper functions
#######################################

### Get the state space (as a set) from a given model.
def getStateSpace(model):
    return(set(map(int,array(list(model.keys())).flatten())))

### Check a model for dynamical reversibility, i.e. whether
### for every forward transition, there is a backward transition
def isReversible(model):
    ### First transform the edge list into a set of sets,
    ### collapsing edges and their reverse
    uniquetransitions=set(map(frozenset,list(model.keys())))
    ### The above is a random choice for forward
    forwtransitions=list(map(list,list(uniquetransitions)))
    ### Initiate the backward transitions as a deep copy of the forward
    backtransitions=deepcopy(forwtransitions)
    ### revert backward edges
    for edge in backtransitions:
        edge.reverse()
    ### transform back to list of tuples and concatenate
    transitions = list(map(tuple,forwtransitions))
    transitions.extend(list(map(tuple,backtransitions)))
    ### check whether reconstructed transitions and
    ### the model keys agree (as sets)
    return(set(transitions) == set(model.keys()))

def isIndexed(model):
    space = getStateSpace(model)
    if( space == set(range(len(space))) ):
        return(True)
    return(False)

### A pair model and chords can be inconsistent.
### Here we check for Integrity and Consistency.
def isConsistent(model,chords):

    ### number of states for later use
    N = len(getStateSpace(model))

    ### Check dynamical reversibility
    if( not isReversible(model) ):
        print(" ERROR:  Given model does not pass the reversibility test.")
        return(False)

    ### Check whether the model has an even number of transitions.
    if( len(list(model.keys())) % 2 != 0 ):
        print(" ERROR:  Given models has an odd number of transition rates.")
        print(" ERROR:  Is it dynamically reversible?")
        return(False)

    ### Assert that len(chords) == len(model.keys)/2 - N + 1
    if( len(chords) > len(list(model.keys()))/2 - N + 1 ):
        print(" ERROR:  Given number of chords is too large.")
        return(False)
    if( len(chords) < len(list(model.keys()))/2 - N + 1 ):
        print(" ERROR:  Given number of chords is too small.")
        return(False)

    ### Check whether the state space is an integer range starting at 0
    if( not isIndexed(model) ):
        print(" ERROR:  Given state space is not an "
                    "integer range starting at 0.")
        return(False)

    ### Assert that chords are a subset of model.keys
    if( not set(chords).issubset(list(model.keys())) ):
        print(" ERROR:  Given set of chords is not contained in given model.  ")
        return(False)

    ### Verify that chords are consistent in themselves, i.e. no doubles when disregarding orientation
    if( len(chords) != len(set(map(frozenset,chords))) ):
        print(" ERROR:  List of chords is not consistent as a set of undirected edges.  ")
        return(False)

    ### Otherwise: Integrity probably OK
    ### Note: These tests do not check whether graph w/o chords actually is a tree.
    return(True)

    
#######################################
# Cumulants possibly with simplifications
#######################################


# 'model' describes the topology
# 'chords' are a list of chords [(0,1),(0,2)]
# 'param' is an optional substitution list for the parametrization of the model
# 'simp' is an optional substitution list for simplifications
# 'unsimp' should revert 'simp': expression.sub(simp).sub(unsimp) == expression

def getCumulants(model, chords, param=[], simp=[], unsimp=[]):

    doSimplify=not(simp==[] and unsimp == [])
    start_time = time.time()

    if( not isConsistent(model,chords) ):
        print(" ERROR:  Model and chords are not correct or inconsistent.  ")
        return( False )

    ### number of states
    N = len(getStateSpace(model))

    ### The first Betti number (cyclomatic number) of the model, equals the number of chords
    B = len(chords)
    
    ### Generate transition matrix W from model
    W = zeros(N) 
    for edge in model:
        W[edge] = (model[edge])
    for j in range(N):
        W[j,j] = -sum(W[j,i] for i in range(N)) 
        
    #display(W)
        
    ### Generate tilted matrix, from copy of W
    Wq = W[:,:]
    q=zeros(B,1)
    for i in range(B):
        name = 'q_{{{0}}}'.format(i)
        q[i] = symbols(name)
        Wq[chords[i]] =  W[chords[i]]*exp(q[i])
        Wq[chords[i][::-1]] =  W[chords[i][::-1]]*exp(-q[i])
        
    #display(Wq)
    
    ### Find coefficients of characteristic polynomial
    print(("--- %s seconds ---" % (time.time() - start_time)))
    print("Start calculating characteristic polynomial")
    a = simplify(Wq.berkowitz()[-1][::-1])  # symbolic simplification is *crucial* here!
        
    ### Initialize current vector and covariance matrix
    c = zeros(B,1)
    C = zeros(B,B)
    
    ### Calculate current vector
    print(("--- %s seconds ---" % (time.time() - start_time)))
    print("Start calculating current vector")
    for i in range(B):
        c[i] = -(diff(a[0],q[i])/a[1]).ratsimp() ## populate current vector
        
    for i in range(B):
        c=c.subs(q[i],0) ## subsitute q=0
    
    if(doSimplify):
        c = c.subs(param)
        
        print(("--- %s seconds ---" % (time.time() - start_time)))
        print("Start simplifying current vector")
        c = simplify(c.subs(simp)) #simplify cancel
    
    ### Calculate co-variance matrix
    print(("--- %s seconds ---" % (time.time() - start_time)))
    print("Start calculating covariance matrix")
    
    ### Do in-place parametrization, before simplification, if latter is demanded
    if(doSimplify):
        for i in range(B):
            for j in range(i+1):
                t1 = ratsimp(  diff(a[0],q[i],q[j]).subs(param).subs(simp) )
                t2 = ratsimp( (diff(a[1],q[i])*c[j]).subs(param).subs(simp) )
                t3 = ratsimp( (diff(a[1],q[j])*c[i]).subs(param).subs(simp) )
                t4 = ratsimp( (2*a[2]*c[i]*c[j]).subs(param).subs(simp) )
                t5 = ratsimp( a[1].subs(param).subs(simp) )
                C[i,j] = -(t1 + t2 + t3 + t4)/t5
    else:
        for i in range(B):
            for j in range(i+1):
                t1 = (  diff(a[0],q[i],q[j]) )#.ratsimp()
                t2 = ( (diff(a[1],q[i])*c[j]) )#.ratsimp()
                t3 = ( (diff(a[1],q[j])*c[i]) )#.ratsimp()
                t4 = ( (2*a[2]*c[i]*c[j]) )#.ratsimp()
                t5 = ( a[1] )#.ratsimp()
                C[i,j] = -(t1 + t2 + t3 + t4)/t5

    ### Populate Covariance Matrix
    
    ## If no simplification, perform parametrization now
    if(not doSimplify):
        c = c.subs(param)

    ### simplification of the expectation should be safe to do, in any case
    c = simplify(c)

    ### simplification of covariance is possibly very time consuming
    C = C.subs(param)
    if(doSimplify):
        print(("--- %s seconds ---" % (time.time() - start_time)))
        print("Start simplifying covariance matrix")
        C = simplify(C)
    
    ### subsitute q=0 into covariance
    for i in range(B):
        C = C.subs(q[i],0)    

    
    print(("--- %s seconds ---" % (time.time() - start_time)))
    print("All Done")
    
    ### Symmetrize covariance matrix:
    for i in range(B):
        for j in range(i):
            C[j,i] = C[i,j]
    
    ### return unsimplified expecation and covariance
    return( [c.subs(unsimp),C.subs(unsimp)]) 


##########################################################################
# Explicit calculation of the cumulants via the SCGF for a two-state model
# Does not work like this. Skrews up parameters. Calculate it directly
#########################################################################

def getSCGF(tiltedGenerator):
    trace = tiltedGenerator.trace()
    det = tiltedGenerator.det()
    scgf = trace/2 + sqrt(trace**2/4-det) #positive sign give largest EV
    return( scgf.simplify() )


