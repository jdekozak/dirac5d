# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, './sympy')

#computer algebra system
from sympy import *
from sympy.matrices import *
from sympy.galgebra.ga import Ga
from sympy.galgebra.printer import *

########################################################################
#ALGEBRA & DEFINITIONS
########################################################################
#Clifford(1,4)
#Flat space, no metric, just signature
#All constants = 1
metric=[1
        ,-1
        ,-1
        ,-1
        ,-1]
#Dimensions
variables = (t, x, y, z, w) = symbols('t x y z w', real=True)
myBasis='gamma_t gamma_x gamma_y gamma_z gamma_w'
#Algebra
sp5d = Ga(myBasis, g=metric, coords=variables,norm=True)
(gamma_t, gamma_x, gamma_y, gamma_z, gamma_w) = sp5d.mv()
(grad, rgrad) = sp5d.grads()
