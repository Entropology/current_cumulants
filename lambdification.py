from sympy import symbols, diff, lambdify, Rational, simplify, N

import cumulants

### Import the kinesin model and its parameters
from models import *

# dummy variables -- lambdify does not like the latex expressions in symbol identifiers
x, y = symbols("x, y", real=True)
dummy = [(f,x),(mu,y)]

##########################################
#  Lau et. al., PRL 99 (2007)            #
##########################################

def lau():

    half = Rational(1,2)

    # calculate SCGF
    scgfLa = cumulants.getSCGF(TiltedWLa)

    # and get the cumulants
    cumDisLa = diff(scgfLa, l)
    cumDisLa = half*cumDisLa.subs([(l,0),(g,0)]) # because velo is measured in d=L/2

    cumDisDisLa = diff(scgfLa, l, l) 
    cumDisDisLa = (half**2)*cumDisDisLa.subs([(l,0),(g,0)]) #  because diffu is measured in d^2=(L/2)^2

    respDisDisLa = 2*diff(cumDisLa,f)/cumDisDisLa  #  non-dimensionalized quantity, obtained by deriving the right nondim velocity with the correct non-dimensionalized force

    # and lambdify them

    return [ lambdify( (x,y), N(thing).subs(dummy).subs(valsubs),"numpy")\
            for thing in (cumDisLa, 0.5*cumDisDisLa, respDisDisLa) ]


##########################################
#  Liepelt, Lipowsky PRL 98 (2007)       #
##########################################

def ll(quick=True):

    cums6_exact = cumulants.getCumulants(model6State,[(1,4),(3,4)], kinesin6_exact)

    vel6_exact = cums6_exact[0][0]
    hyd6_exact = vel6_exact + 2 * cums6_exact[0][1]
    dif6_exact = .5 * cums6_exact[1][0,0]

#    if (not quick):
#        for thing in [dif6_exact]:
#            thing = thing.subs(logargs).ratsimp().subs(expargs)
#            thing = thing.subs(logargs).cancel().subs(expargs)

    coupling6_exact = hyd6_exact/vel6_exact
    invfano6_exact = 2*vel6_exact/dif6_exact

    response6_exact = -diff(vel6_exact,f)
    tmech6_exact = response6_exact / dif6_exact

    for thing in vel6_exact, hyd6_exact, coupling6_exact:
        thing = thing.subs(logargs).ratsimp().subs(expargs)
        thing = thing.subs(logargs).cancel().subs(expargs)
        thing = thing.subs(logargs).simplify().subs(expargs)


##########################################
#  Altaner,Wachtel,Vollmer               #
##########################################

    if(quick):
        cums4_exact = cumulants.getCumulants(model4State,[(0,2),(1,2)], kinesin4_exact)
    else:
        cums4_exact = cumulants.getCumulants(model4State,[(0,2),(1,2)], kinesin4_exact, logargs, expargs) 

    vel4_exact = cums4_exact[0][0]
    hyd4_exact = vel4_exact + 2 * cums4_exact[0][1]
    coupling4_exact = hyd4_exact/vel4_exact

    for thing in vel4_exact, hyd4_exact, coupling4_exact:
        thing = thing.subs(logargs).ratsimp().subs(expargs)
        thing = thing.subs(logargs).cancel().subs(expargs)
        thing = thing.subs(logargs).simplify().subs(expargs)

    dif4_exact = .5 * cums4_exact[1][0,0]

    invfano4_exact =  vel4_exact/(2*dif4_exact)
    response4_exact = -diff(vel4_exact,f)
    tmech4_exact = response4_exact / dif4_exact

    #only uncomment this if you have too much time
    #for thing in dif4_exact, invfano6_exact, response4_exact:
    #    thing = thing.subs(logargs).ratsimp().subs(expargs)
    #    thing = thing.subs(logargs).cancel().subs(expargs)

##########################################
#  Comparison of 4-state and 6-state     #
##########################################

    vel_relerr_exact = vel4_exact/vel6_exact - 1
    hyd_relerr_exact = hyd4_exact/hyd6_exact - 1

    for thing in vel_relerr_exact, hyd_relerr_exact:
        thing = thing.subs(logargs).ratsimp().subs(expargs)
        thing = thing.subs(logargs).cancel().subs(expargs)

    dif_relerr_exact = dif4_exact/dif6_exact - 1

##########################################
#  Output of everything                  #
##########################################

    return [ lambdify( (x,y), N(thing).subs(dummy), "numpy" )\
        for thing in ( \
        vel6_exact, hyd6_exact, dif6_exact, coupling6_exact, invfano6_exact, tmech6_exact, \
        vel4_exact, hyd4_exact, dif4_exact, coupling4_exact, invfano4_exact, tmech4_exact, \
        vel_relerr_exact, hyd_relerr_exact, dif_relerr_exact ) ]
