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

#Imaginary unit
imag=gamma_w
imag.texLabel='i'
#Associative Hyperbolic Quaternions
ihquat=gamma_t
jhquat=gamma_t*gamma_x*gamma_y*gamma_z*gamma_w
khquat=gamma_x*gamma_y*gamma_z*gamma_w
ihquat.texLabel='\\mathbf{i}'
jhquat.texLabel='\\mathbf{j}'
khquat.texLabel='\\mathbf{k}'
#Quaternions
iquat=gamma_y*gamma_z
jquat=gamma_z*gamma_x
kquat=gamma_x*gamma_y
iquat.texLabel='\\boldsymbol{\\mathit{i}}'
jquat.texLabel='\\boldsymbol{\\mathit{j}}'
kquat.texLabel='\\boldsymbol{\\mathit{k}}'
