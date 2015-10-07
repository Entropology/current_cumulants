#import sympy as sy
from sympy import zeros, exp, sqrt, log, symbols, S, Rational

# Define some symbols

# The symbols 'w' represent the transition rates for the Markov models

w12, w21, w13, w31, w14, w41, w24, w42, w23, w32, w34, w43,\
w16, w61, w25, w52, w45, w54, w56, w65 = \
    symbols(
         """w_{12}, w_{21}, w_{13}, w_{31}, w_{14}, w_{41},
            w_{24}, w_{42}, w_{23}, w_{32}, w_{34}, w_{43},
            w_{16}, w_{61}, w_{25}, w_{52}, w_{54}, w_{45},
            w_{56}, w_{65}""", real=True, positive=True)

# f and mu are the thermodynamic forces driving the kinesin models
f, mu = symbols("f, \Delta\mu", real=True)

# x and y are coordinates for 2D plots
x, y = symbols("x, y", real=True)

# u and v are used as exponentials of f and mu, thus only positive,
# which is important for symbolic simplification
u, v = symbols("u, v", positive=True)

# In the following, a model is a dictionary

# this model is a general 4-state model with all 6 possible edges.
# can be useful for quick checks, even though not used anymore
gen4State = {(0,1): w12, (1,0): w21,\
               (0,2): w13, (2,0): w31,\
               (0,3): w14, (3,0): w41,\
               (1,3): w24, (3,1): w42,\
               (1,2): w23, (2,1): w32,\
               (2,3): w34, (3,2): w43 \
               }


############################################
#         Kinesin models based on          #
#     Liepelt, Lipowsky, PRL 98 (2007)     #
############################################


# The 6 state model already including the topology
# for kinesin by Liepelt+Lipowsky.

model6State = {(0,1): w12, (1,0): w21,\
               (0,5): w16, (5,0): w61,\
               (1,2): w23, (2,1): w32,\
               (1,4): w25, (4,1): w52,\
               (2,3): w34, (3,2): w43,\
               (3,4): w45, (4,3): w54,\
               (4,5): w56, (5,4): w65\
               }

# A 4 state model already including the topology
# for kinesin as described in Altaner+Vollmer2014.

model4State = {(0,1): w12, (1,0): w21,\
               (0,2): w13, (2,0): w31,\
               (0,3): w14, (3,0): w41,\
               #(1,3): w24, (3,1): w42,\
               (1,2): w23, (2,1): w32,\
               (2,3): w34, (3,2): w43 \
               }

##########################################
## A numerical version of the kinesin4 model
##########################################

# Equilibrium constant hydrolysis reaction
K = 4.9e11

## Load dependence parameters
xi1=0.15
xi2=0.25
theta = 0.65

# First order rates
k13 = 3e5
k31 = 0.24

k14 = 100.0
k41 = 2.0
k43 = 2.51608e-6 # fitted rate
k34 = K*k43*k14*k31/(k41*k13) #should be ~49.3
k32 = (k31/k13)**2*k14 #should be 6.4e-11
k23 = k41
k21 = k43
k12 = k34

## Mechanical transition
w13p = k13*exp(-theta*f)
w31p = k31*exp((1-theta)*f)


## Chemical transitions
w12p = k12*2/(1+exp(xi1*f)) #ADP, P attachment
w21p = k21*2/(1+exp(xi1*f))
w23p = k23*2/(1+exp(xi2*f))/K*exp(mu) #ATP attachment
w32p = k32*2/(1+exp(xi2*f))
w34p = k34*2/(1+exp(xi1*f)) #ADP, P attachment
w43p = k43*2/(1+exp(xi1*f))
w41p = k41*2/(1+exp(xi2*f))/K*exp(mu) #ATP attachment
w14p = k14*2/(1+exp(xi2*f))

# substitution list
kinesin4_numeric = [(w12,w12p),(w21,w21p),(w13,w13p),(w31,w31p),\
                    (w14,w14p),(w41,w41p),(w23,w23p),(w32,w32p),\
                    (w34,w34p),(w43,w43p),(w24,0),(w42,0)]



##########################################
## Analytical version of the kinesin4 Model
##########################################


# Equilibrium constant hydrolysis reaction
K = S(490000000000) # 4.9e11

## Load dependence parameters
xi1=Rational(15, 100) #0.15
xi2=Rational(25, 100) #0.25
theta = Rational(65, 100) #0.65

# First order rates
k13 = S(300000)  #SympyFloat
k31 = Rational(24,100) #0.24

k14 = S(100)
k41 = S(2)
k43 = Rational(251608,10**11) # fitted rate = 2.51608e-6
k34 = K*k43*k14*k31/(k41*k13) #should be ~49.3
k32 = (k31/k13)**2*k14 #should be 6.4e-11
k23 = k41
k21 = k43
k12 = k34

## Mechanical transition
w13p = k13*exp(-theta*f)
w31p = k31*exp((1-theta)*f)


## Chemical transitions
w12p = k12*2/(1+exp(xi1*f)) #ADP, P attachment
w21p = k21*2/(1+exp(xi1*f))
w23p = k23*2/(1+exp(xi2*f))/K*exp(mu) #ATP attachment
w32p = k32*2/(1+exp(xi2*f))
w34p = k34*2/(1+exp(xi1*f)) #ADP, P attachment
w43p = k43*2/(1+exp(xi1*f))
w41p = k41*2/(1+exp(xi2*f))/K*exp(mu) #ATP attachment
w14p = k14*2/(1+exp(xi2*f))

# substitution list
kinesin4_exact = [(w12,w12p),(w21,w21p),(w13,w13p),(w31,w31p),\
                  (w14,w14p),(w41,w41p),(w23,w23p),(w32,w32p),\
                  (w34,w34p),(w43,w43p),(w24,0),(w42,0)]


##########################################
## Numerical version of the kinesin6 model
##########################################


# Equilibrium constant hydrolysis reaction
K = 4.9e11

## Load dependence parameters
xi16= 0.15
xi12= 0.25
theta = 0.65

# First order rates
k25 = 3e5  #SympyFloat
k52 = 0.24

k56 = 100.
k65 = 5./196.
k16 = 2./100.
k12 = 2.
k21 = k56
k23 = k56
k34 = k56
k32 = k65
k43 = k16
k45 = k12
k54 = k21*(k52/k25)**2
k61 = k56

## Mechanical transition
w25p = k25*exp(-theta*f)
w52p = k52*exp((1-theta)*f)

## Chemical transitions
w12p = k12*2/(1+exp(xi12*f))/K*exp(mu) # ATP attachment
w21p = k21*2/(1+exp(xi12*f))
w16p = k16*2/(1+exp(xi16*f)) # P attachment
w61p = k61*2/(1+exp(xi16*f))

w45p = k45*2/(1+exp(xi12*f))/K*exp(mu) # ATP attachment
w54p = k54*2/(1+exp(xi12*f))
w43p = k43*2/(1+exp(xi16*f)) # P attachment
w34p = k34*2/(1+exp(xi16*f))

w23p = k23*2/(1+exp(xi16*f))
w32p = k32*2/(1+exp(xi16*f)) # ADP attachment
w56p = k56*2/(1+exp(xi16*f))
w65p = k65*2/(1+exp(xi16*f)) # ADP attachment


# substitution list
kinesin6_numeric=[(w12,w12p),(w21,w21p),(w16,w16p),(w61,w61p),\
                  (w45,w45p),(w54,w54p),(w23,w23p),(w32,w32p),\
                  (w34,w34p),(w43,w43p),(w25,w25p),(w52,w52p),\
                  (w56,w56p),(w65,w65p)]

##########################################
## Analytical version of the kinesin6 model
##########################################


# Equilibrium constant hydrolysis reaction
K = S(490000000000) # 4.9e11

## Load dependence parameters
xi16=Rational(15, 100) #0.15
xi12=Rational(25, 100) #0.25
theta = Rational(65, 100) #0.65

# First order rates
k25 = S(300000)  #SympyFloat
k52 = Rational(24,100) #0.24

k56 = S(100)
k65 = Rational(5,196)
k16 = Rational(2,100)
k12 = S(2)
k21 = k56
k23 = k56
k34 = k56
k32 = k65
k43 = k16
k45 = k12
k54 = k21*(k52/k25)**2
k61 = k56

## Mechanical transition
w25p = k25*exp(-theta*f)
w52p = k52*exp((1-theta)*f)

## Chemical transitions
w12p = k12*2/(1+exp(xi12*f))/K*exp(mu) # ATP attachment
w21p = k21*2/(1+exp(xi12*f))
w16p = k16*2/(1+exp(xi16*f)) # P attachment
w61p = k61*2/(1+exp(xi16*f))

w45p = k45*2/(1+exp(xi12*f))/K*exp(mu) # ATP attachment
w54p = k54*2/(1+exp(xi12*f))
w43p = k43*2/(1+exp(xi16*f)) # P attachment
w34p = k34*2/(1+exp(xi16*f))

w23p = k23*2/(1+exp(xi16*f))
w32p = k32*2/(1+exp(xi16*f)) # ADP attachment
w56p = k56*2/(1+exp(xi16*f))
w65p = k65*2/(1+exp(xi16*f)) # ADP attachment


# substitution list
kinesin6_exact = [(w12,w12p),(w21,w21p),(w16,w16p),(w61,w61p),\
                  (w45,w45p),(w54,w54p),(w23,w23p),(w32,w32p),\
                  (w34,w34p),(w43,w43p),(w25,w25p),(w52,w52p),\
                  (w56,w56p),(w65,w65p)]

########################################
# 
########################################

# this substitution transforms our expressions to rational functions
# of u and v for efficient simplification
logargs = [ (f,20*log(u)), (mu,log(v)) ]

# this substitution transforms the rational functions back to
# the physiological parameters
expargs = [ (u,exp(f/20)), (v,exp(mu)) ]


# this substitution transforms our expressions to rational functions
# of u and v for efficient simplification
alt_logargs=[ (f,10*log(u)), (mu,log(v)) ]

# this substitution transforms the rational functions back
# to the physiological parameters
alt_expargs=[ (u,exp(f/10)), (v,exp(mu)) ]




############################################
#         Kinesin models based on          #
#   Lau, Lacoste, Mallick PRL 99 (2007)    #
############################################


# Numerical model parameters
(e, a, aa, om, omm, tap, tam, tbp, tbm) =\
        symbols("""
        epsilon,
        alpha,
        alpha',
        omega,
        omega',
        Theta_a^+,
        Theta_a^-,
        Theta_b^+,
        Theta_b^-
        """ , positive=True)

valsubs = [(e,10.81), (a,0.57), (aa,1.3e-6), (om,3.5) , (omm,108.15),\
           (tap,0.25), (tam,1.83), (tbp,0.08), (tbm,-0.16)]

# Symbolic transition rates

wlBm = a*exp(-tbm*f)
wlBn = om*exp(-tbm*f)
wrAp = a*exp(-e + mu + tap*f)
wrAn = om*exp(-e + tap*f)
wlAp = aa*exp(-e + mu - tam *f)
wlAn = omm*exp(-e - tam*f)
wrBm = aa*exp(tbp*f)
wrBn = omm*exp(tbp*f)

# effective rates
wrA = (wrAp + wrAn)
wlA = (wlAp + wlAn)
wrB = (wrBm + wrBn)
wlB = (wlBm + wlBn)

# tilted matrix
l,g = symbols("lambda, gamma")
TiltedWLa = zeros(2)
TiltedWLa[0,0] = -wrA - wlA
TiltedWLa[0,1] = exp(l)*(wlBm*exp(g) + wlBn) + exp(-l)*(wrBm*exp(g) + wrBn)
TiltedWLa[1,0] = exp(l)*(wlAp*exp(-g) + wlAn) + exp(-l)*(wrAp*exp(-g) + wrAn)
TiltedWLa[1,1] = -wrB - wlB

# Now the conventions regarding the definition
# of the parameters agree between the models
TiltedWLa = TiltedWLa.subs(f,-f/2) 

trace = TiltedWLa.trace()
det = TiltedWLa.det()
scgfLa = trace/2 + sqrt(trace**2/4-det) #positive sign give largest EV
