# -*- coding: utf-8 -*-
#computer algebra system
from sympy import *
from sympy.matrices import *
from sympy.galgebra.GA import *
import sympy.galgebra.latex_ex as tex

set_main(sys.modules[__name__])

outputTex=True

########################################################################
#ALGEBRA & DEFINITIONS
########################################################################
#Clifford(1,4)
metric = '1  0  0  0  0,'+\
         '0 -1  0  0  0,'+\
         '0  0 -1  0  0,'+\
         '0  0  0 -1  0,'+\
         '0  0  0  0 -1'

#Dimensions
vars = make_symbols('t x y z w')
myBasis='gamma_t gamma_x gamma_y gamma_z gamma_w'
MV.setup(myBasis,metric,True,vars)
if outputTex:
    tex.Format('1 1 1 1')
else:
    MV.set_str_format(0)
#Imaginary unit
imag=gamma_w
#CoQuaternions
icquat=gamma_t
jcquat=gamma_t*gamma_x*gamma_y*gamma_z
kcquat=gamma_x*gamma_y*gamma_z
#Quaternions
iquat=gamma_y*gamma_z
jquat=gamma_z*gamma_x
kquat=gamma_x*gamma_y

def CheckProperties(i,j,k, ilabel, jlabel, klabel):
    print ilabel + '=', i
    print jlabel + '=', j
    print klabel + '=', k
    print ilabel + '^2=', i*i
    print jlabel + '^2=', j*j
    print klabel + '^2=', k*k
    print ilabel + jlabel + '=', i*j
    print jlabel + ilabel + '=', j*i
    print jlabel + klabel + '=', j*k
    print klabel + jlabel + '=', k*j
    print klabel + ilabel + '=', k*i
    print ilabel + klabel + '=', i*k
    print ilabel + jlabel + klabel + '=', i*j*k

def CheckGlobals():
    '''Check algebra definitions'''
    print 'CoQuaternions'
    CheckProperties(icquat,jcquat,kcquat, 'i_{cq}', 'j_{cq}', 'k_{cq}')
    print
    print 'Quaternions'
    CheckProperties( iquat, jquat, kquat, 'i_{q}', 'j_{q}', 'k_{q}')
    print
    print 'Imaginary Unit'
    print 'i = ',imag
    print 'i^2=', imag*imag
    print
    print 'PseudoScalar'
    print 'j_{cq}\\ i = ', jcquat*imag
    print 'j_{cq}\\ i\\ j_{cq}\\ i = ', jcquat*imag*jcquat*imag


#Computation with the gradient
def ShowGrad(af):
    """Dirac equation is the geometric algebra gradient"""
    i = 1
    for f in af:
        print '{\\nabla}f_'+str(i)+' = ', f.grad()
        i += 1

def CheckComputedGrad(Grad):
    i=1
    for grad in Grad:
        print '{\\nabla}f_'+str(i)+'^2=',grad*grad
        i=i+1

def CheckNil(Grad,A_i):
    print
    i=1
    for (g,a) in zip(Grad,A_i):
        sq=a*g
        sq.expand()
        print 'A_'+str(i)+' * {\\nabla}f_'+str(i)+'=', sq
        i+=1

def relativityEnergyConservation(aMultiVector):
    return aMultiVector.subs(-m**2-p_x**2-p_y**2-p_z**2,-E**2).subs(m**2+p_x**2+p_y**2+p_z**2,E**2)

########################################################################
#PHYSICS
########################################################################
parms = make_symbols('m E p_x p_y p_z e phi Ax Ay Az')
r = [x, y, z]
rquat = [iquat, jquat, kquat]
#vecteur d'onde (propagation spatiale)
p = [p_x, p_y, p_z]
psigma = S(0)
for (dim, var) in zip(p, r):
    psigma += var * dim
pmulti = S(0)
for (dim, var) in zip(p, rquat):
    pmulti += var * dim

r = [  x,  y,  z,      w]
#vecteur d'onde (propagation spatiale) en milieu absorbant ou amplificateur
k = [ p_x, p_y, p_z, imag*m]
kdotr = S(0)
for (dim, var) in zip(k, r):
    kdotr += var * dim
#Fonction de l'exponentielle
f1 = -imag * (E * t - kdotr)

k = [-p_x, -p_y, -p_z, -imag * m]
kdotr = S(0)
for (dim, var) in zip(k, r):
    kdotr += var * dim
f2 = -imag * (E * t - kdotr)
k = [ p_x, p_y, p_z, -imag * m]
kdotr = S(0)
for (dim, var) in zip(k, r):
    kdotr += var * dim
f3 = -imag * (E * t - kdotr)
k = [-p_x, -p_y, -p_z, imag * m]
kdotr = S(0)
for (dim, var) in zip(k, r):
    kdotr += var * dim
f4 = -imag * (E * t - kdotr)
f = [S(1)*f1, S(1)*f2, S(1)*f3, S(1)*f4]

#Constantes 'simples' x4
K1=jcquat*imag*(kcquat*E-jcquat*m+icquat*pmulti)
K1label='K_1 = j_{cq} i (k_{cq} E - j_{cq} m + i_{cq} \\mathbf{p})'
K2=jcquat*imag*(kcquat*E+jcquat*m-icquat*pmulti)
K2label='K_2 = j_{cq} i (k_{cq} E + j_{cq} m - i_{cq} \\mathbf{p})'
K3=jcquat*imag*(kcquat*E+jcquat*m+icquat*pmulti)
K3label='K_3 = j_{cq} i (k_{cq} E + j_{cq} m + i_{cq} \\mathbf{p})'
K4=jcquat*imag*(kcquat*E-jcquat*m-icquat*pmulti)
K4label='K_4 = j_{cq} i (k_{cq} E - j_{cq} m - i_{cq} \\mathbf{p})'
K=[K1, K2, K3, K4]

#Constantes 'combinees' x16
#premiere forme
A1a=(K1-K4*icquat-K3*jcquat+K2*kcquat)
A1alabel='A_a^1 = K_1-K_4 i_{cq}-K_3 j_{cq}+K_2 k_{cq}'
A2a=(K2-K3*icquat-K4*jcquat+K1*kcquat)
A2alabel='A_a^2 = K_2-K_3 i_{cq}-K_4 j_{cq}+K_1 k_{cq}'
A3a=(K3-K2*icquat-K1*jcquat+K4*kcquat)
A3alabel='A_a^3 = K_3-K_2 i_{cq}-K_1 j_{cq}+K_4 k_{cq}'
A4a=(K4-K1*icquat-K2*jcquat+K3*kcquat)
A4alabel='A_a^4 = K_4-K_1 i_{cq}-K_2 j_{cq}+K_3 k_{cq}'
#deuxieme forme
A1b=(K1+K4*icquat+K3*jcquat+K2*kcquat)
A1blabel='A_b^1 = K_1+K_4 i_{cq}+K_3 j_{cq}+K_2 k_{cq}'
A2b=(K2+K3*icquat+K4*jcquat+K1*kcquat)
A2blabel='A_b^2 = K_2+K_3 i_{cq}+K_4 j_{cq}+K_1 k_{cq}'
A3b=(K3+K2*icquat+K1*jcquat+K4*kcquat)
A3blabel='A_b^3 = K_3+K_2 i_{cq}+K_1 j_{cq}+K_4 k_{cq}'
A4b=(K4+K1*icquat+K2*jcquat+K3*kcquat)
A4blabel='A_b^4 = K_4+K_1 i_{cq}+K_2 j_{cq}+K_3 k_{cq}'
#troisieme forme
A1c=(-K1+K4*icquat-K3*jcquat+K2*kcquat)
A1clabel='A_c^1 = -K_1+K_4 i_{cq}-K_3 j_{cq}+K_2 k_{cq}'
A2c=(-K2+K3*icquat-K4*jcquat+K1*kcquat)
A2clabel='A_c^2 = -K_2+K_3 i_{cq}-K_4 j_{cq}+K_1 k_{cq}'
A3c=(-K3+K2*icquat-K1*jcquat+K4*kcquat)
A3clabel='A_c^3 = -K_3+K_2 i_{cq}-K_1 j_{cq}+K_4 k_{cq}'
A4c=(-K4+K1*icquat-K2*jcquat+K3*kcquat)
A4clabel='A_c^4 = -K_4+K_1 i_{cq}-K_2 j_{cq}+K_3 k_{cq}'
#quatrieme forme
A1d=(-K1-K4*icquat+K3*jcquat+K2*kcquat)
A1dlabel='A_d^1 = -K_1-K_4 i_{cq}+K_3 j_{cq}+K_2 k_{cq}'
A2d=(-K2-K3*icquat+K4*jcquat+K1*kcquat)
A2dlabel='A_d^2 = -K_2-K_3 i_{cq}+K_4 j_{cq}+K_1 k_{cq}'
A3d=(-K3-K2*icquat+K1*jcquat+K4*kcquat)
A3dlabel='A_d^3 = -K_3-K_2 i_{cq}+K_1 j_{cq}+K_4 k_{cq}'
A4d=(-K4-K1*icquat+K2*jcquat+K3*kcquat)
A4dlabel='A_d^4 = -K_4-K_1 i_{cq}+K_2 j_{cq}+K_3 k_{cq}'
#
Aa=[A1a, A2a, A3a, A4a]
Ab=[A1b, A2b, A3b, A4b]
Ac=[A1c, A2c, A3c, A4c]
Ad=[A1d, A2d, A3d, A4d]
#
A1=[A1a, A1b, A1c, A1d]
A2=[A2a, A2b, A2c, A2d]
A3=[A3a, A3b, A3c, A3d]
A4=[A4a, A4b, A4c, A4d]
#Equation d'onde (Onde Plane Progressive Monochromatique)
#psi=A1a*exp(f1)

########################################################################
#MAIN DIRAC
########################################################################
print('Algebra is Clifford(1,4) over the reals : $Cl_{1,4}(\\mathbb{R})$')
print('The five dimensions are $t$, $x$, $y$, $z$, $w$')
print('It defines two sets of quaternions with one imaginary unit and a pseudoscalar')
print('')
CheckGlobals()
print('')
print('Gradient')
print('{\\nabla}=({\gamma}_t {\\partial}/{\\partial t}+{\gamma}_x {\\partial}/{\\partial x}+{\gamma_y} {\\partial}/{\\partial y}+{\gamma_z} {\\partial}/{\\partial z}+{\gamma_w} {\\partial}/{\\partial w})')
print
print('Wavefunction : $A$ is a constant and $f$ is a function of $t$, $x$, $y$, $z$, $w$')
print('{\psi}=A e^f')
print('{\\nabla}{\psi}=A ({\\nabla}f) e^f')
print
print('The following symbols are defined : (a positive value is NOT required)')
print('$E$ is for energy, $E \\in \\mathbb{R}$')
print('$m$ is for mass, $m \\in \\mathbb{R}$')
print('$\\mathbf{p}$ is the momentum. $p_x, p_y, p_z \\in \\mathbb{R}$')
print
print('Exponential function $f$')
i=1
for aF in f:
    print 'f_'+str(i)+'=', aF
    i += 1
print
print('Gradient for $f$')
ShowGrad(f)
print
print('Square of the gradient')
Grad=[f1.grad(),
      f2.grad(),
      f3.grad(),
      f4.grad()]
CheckComputedGrad(Grad)
print('')
print('Dirac')
print('0=({\gamma}_0 {\\partial}/{\\partial t}+{\gamma}_1 {\\partial}/{\\partial x}+{\gamma}_2 {\\partial}/{\\partial y}+{\gamma}_3 {\\partial}/{\\partial z}+i\\ m) {\psi}')
print
print('With the above gradients, identify the Dirac algebra aka gamma matrices')
#print('http://en.wikipedia.org/wiki/Gamma_matrices#Normalisation')
gamma_0 =  gamma_t*gamma_w
gamma_1 = -gamma_x*gamma_w
gamma_2 = -gamma_y*gamma_w
gamma_3 = -gamma_z*gamma_w
gamma_5 = imag*gamma_0*gamma_1*gamma_2*gamma_3
print '{\gamma}_0 = ', gamma_0
print '{\gamma}_0^2 = ', gamma_0*gamma_0
print '{\gamma}_1 = ', gamma_1
print '{\gamma}_1^2 = ', gamma_1*gamma_1
print '{\gamma}_2 = ', gamma_2
print '{\gamma}_2^2 = ', gamma_2*gamma_2
print '{\gamma}_3 = ', gamma_3
print '{\gamma}_3^2 = ', gamma_3*gamma_3
print '{\gamma}_5 = ', gamma_5
print '{\gamma}_5^2 = ', gamma_5*gamma_5
print
print('Simple Constants (exactly the gradient)')
print '\\mathbf{p}=', pmulti
i=1
for (k,label) in zip(K,[K1label, K2label, K3label, K4label]):
    print label
    k.convert_from_blades()
    print 'K_'+str(i)+'=', k
    i+=1

print
if outputTex:
    tex.Format('1 1 1 2')
else:
    MV.set_str_format(1)
print('Mixed Constants (built from simple constants and the coquaternions) see details at the end')
for label in [[A1alabel, A2alabel, A3alabel, A4alabel],
              [A1blabel, A2blabel, A3blabel, A4blabel],
              [A1clabel, A2clabel, A3clabel, A4clabel],
              [A1dlabel, A2dlabel, A3dlabel, A4dlabel]
              ]:
    for lab in label:
        print lab
print 'Symmetry $A^1$'
print '[A_a^1, A_b^1, A_c^1, A_d^1]\\ with\\ i_{cq} \\Rightarrow [A_c^4, A_d^4, A_a^4, A_b^4]'
for (label, product) in [('i_{cq} A_a^1 i_{cq} - A_c^4 = ', icquat*A1a*icquat - A4c),
                         ('i_{cq} A_b^1 i_{cq} - A_d^4 = ', icquat*A1b*icquat - A4d),
                         ('i_{cq} A_c^1 i_{cq} - A_a^4 = ', icquat*A1c*icquat - A4a),
                         ('i_{cq} A_d^1 i_{cq} - A_b^4 = ', icquat*A1d*icquat - A4b)
                         ]:
    if product == 0:
        print label + '= 0'
print '[A_a^1, A_b^1, A_c^1, A_d^1]\\ with\\ j_{cq} \\Rightarrow [-A_d^3, -A_c^3, -A_b^3, -A_a^3]'
for (label, product) in [('j_{cq} A_a^1 j_{cq} + A_d^3 = ', jcquat*A1a*jcquat + A3d),
                         ('j_{cq} A_b^1 j_{cq} + A_c^3 = ', jcquat*A1b*jcquat + A3c),
                         ('j_{cq} A_c^1 j_{cq} + A_b^3 = ', jcquat*A1c*jcquat + A3b),
                         ('j_{cq} A_d^1 j_{cq} + A_a^3 = ', jcquat*A1d*jcquat + A3a)
                         ]:
    if product == 0:
        print label + '= 0'
print '[A_a^1, A_b^1, A_c^1, A_d^1]\\ with\\ k_{cq} \\Rightarrow [A_b^2, A_a^2, A_d^2, A_c^2]'
for (label, product) in [('k_{cq} A_a^1 k_{cq} - A_b^2 = ', kcquat*A1a*kcquat - A2b),
                         ('k_{cq} A_b^1 k_{cq} - A_a^2 = ', kcquat*A1b*kcquat - A2a),
                         ('k_{cq} A_c^1 k_{cq} - A_d^2 = ', kcquat*A1c*kcquat - A2d),
                         ('k_{cq} A_d^1 k_{cq} - A_c^2 = ', kcquat*A1d*kcquat - A2c)
                         ]:
    if product == 0:
        print label + '= 0'
print 'Symmetry $A^2$'
print '[A_a^2, A_b^2, A_c^2, A_d^2]\\ with\\ i_{cq} \\Rightarrow [A_c^3, A_d^3, A_a^3, A_b^3]'
for (label, product) in [('i_{cq} A_a^2 i_{cq} - A_b^3 = ', icquat*A2a*icquat - A3c),
                         ('i_{cq} A_b^2 i_{cq} - A_a^3 = ', icquat*A2b*icquat - A3d),
                         ('i_{cq} A_c^2 i_{cq} - A_d^3 = ', icquat*A2c*icquat - A3a),
                         ('i_{cq} A_d^2 i_{cq} - A_c^3 = ', icquat*A2d*icquat - A3b),
                         ]:
    if product == 0:
        print label + '= 0'
print '[A_a^2, A_b^2, A_c^2, A_d^2]\\ with\\ j_{cq} \\Rightarrow [-A_d^4, -A_c^4, -A_b^4, -A_a^4]'
for (label, product) in [('j_{cq} A_a^2 j_{cq} + A_d^4 = ', jcquat*A2a*jcquat + A4d),
                         ('j_{cq} A_b^2 j_{cq} + A_c^4 = ', jcquat*A2b*jcquat + A4c),
                         ('j_{cq} A_c^2 j_{cq} + A_b^4 = ', jcquat*A2c*jcquat + A4b),
                         ('j_{cq} A_d^2 j_{cq} + A_a^4 = ', jcquat*A2d*jcquat + A4a)
                         ]:
    if product == 0:
        print label + '= 0'
print '[A_a^2, A_b^2, A_c^2, A_d^2]\\ with\\ k_{cq} \\Rightarrow [A_b^1, A_a^1, A_d^1, A_c^1]'
for (label, product) in [('k_{cq} A_a^2 k_{cq} - A_b^1 = ', kcquat*A2a*kcquat - A1b),
                         ('k_{cq} A_b^2 k_{cq} - A_a^1 = ', kcquat*A2b*kcquat - A1a),
                         ('k_{cq} A_c^2 k_{cq} - A_d^1 = ', kcquat*A2c*kcquat - A1d),
                         ('k_{cq} A_d^2 k_{cq} - A_c^1 = ', kcquat*A2d*kcquat - A1c)
                         ]:
    if product == 0:
        print label + '= 0'
print 'Symmetry $A^3$'
print '[A_a^3, A_b^3, A_c^3, A_d^3]\\ with\\ i_{cq} \\Rightarrow [A_c^2, A_d^2, A_a^2, A_b^2]'
for (label, product) in [('i_{cq} A_a^3 i_{cq} - A_b^2 = ', icquat*A3a*icquat - A2c),
                         ('i_{cq} A_b^3 i_{cq} - A_a^2 = ', icquat*A3b*icquat - A2d),
                         ('i_{cq} A_c^3 i_{cq} - A_d^2 = ', icquat*A3c*icquat - A2a),
                         ('i_{cq} A_d^3 i_{cq} - A_c^2 = ', icquat*A3d*icquat - A2b),
                         ]:
    if product == 0:
        print label + '= 0'
print '[A_a^3, A_b^3, A_c^3, A_d^3]\\ with\\ j_{cq} \\Rightarrow [-A_d^1, -A_c^1, -A_b^1, -A_a^1]'
for (label, product) in [('j_{cq} A_a^3 j_{cq} + A_d^1 = ', jcquat*A3a*jcquat + A1d),
                         ('j_{cq} A_b^3 j_{cq} + A_c^1 = ', jcquat*A3b*jcquat + A1c),
                         ('j_{cq} A_c^3 j_{cq} + A_b^1 = ', jcquat*A3c*jcquat + A1b),
                         ('j_{cq} A_d^3 j_{cq} + A_a^1 = ', jcquat*A3d*jcquat + A1a)
                         ]:
    if product == 0:
        print label + '= 0'
print '[A_a^3, A_b^3, A_c^3, A_d^3]\\ with\\ k_{cq} \\Rightarrow [A_b^4, A_a^4, A_d^4, A_c^4]'
for (label, product) in [('k_{cq} A_a^3 k_{cq} - A_b^4 = ', kcquat*A3a*kcquat - A4b),
                         ('k_{cq} A_b^3 k_{cq} - A_a^4 = ', kcquat*A3b*kcquat - A4a),
                         ('k_{cq} A_c^3 k_{cq} - A_d^4 = ', kcquat*A3c*kcquat - A4d),
                         ('k_{cq} A_d^3 k_{cq} - A_c^4 = ', kcquat*A3d*kcquat - A4c)
                         ]:
    if product == 0:
        print label + '= 0'
print 'Symmetry $A^4$'
print '[A_a^4, A_b^4, A_c^4, A_d^4]\\ with\\ i_{cq} \\Rightarrow [A_c^1, A_d^1, A_a^1, A_b^1]'
for (label, product) in [('i_{cq} A_a^4 i_{cq} - A_b^1 = ', icquat*A4a*icquat - A1c),
                         ('i_{cq} A_b^4 i_{cq} - A_a^1 = ', icquat*A4b*icquat - A1d),
                         ('i_{cq} A_c^4 i_{cq} - A_d^1 = ', icquat*A4c*icquat - A1a),
                         ('i_{cq} A_d^4 i_{cq} - A_c^1 = ', icquat*A4d*icquat - A1b),
                         ]:
    if product == 0:
        print label + '= 0'
print '[A_a^4, A_b^4, A_c^4, A_d^4]\\ with\\ j_{cq} \\Rightarrow [-A_d^2, -A_c^2, -A_b^2, -A_a^2]'
for (label, product) in [('j_{cq} A_a^4 j_{cq} + A_d^2 = ', jcquat*A4a*jcquat + A2d),
                         ('j_{cq} A_b^4 j_{cq} + A_c^2 = ', jcquat*A4b*jcquat + A2c),
                         ('j_{cq} A_c^4 j_{cq} + A_b^2 = ', jcquat*A4c*jcquat + A2b),
                         ('j_{cq} A_d^4 j_{cq} + A_a^2 = ', jcquat*A4d*jcquat + A2a)
                         ]:
    if product == 0:
        print label + '= 0'
print '[A_a^4, A_b^4, A_c^4, A_d^4]\\ with\\ k_{cq} \\Rightarrow [A_b^3, A_a^3, A_d^3, A_c^3]'
for (label, product) in [('k_{cq} A_a^4 k_{cq} - A_b^3 = ', kcquat*A4a*kcquat - A3b),
                         ('k_{cq} A_b^4 k_{cq} - A_a^3 = ', kcquat*A4b*kcquat - A3a),
                         ('k_{cq} A_c^4 k_{cq} - A_d^3 = ', kcquat*A4c*kcquat - A3d),
                         ('k_{cq} A_d^4 k_{cq} - A_c^3 = ', kcquat*A4d*kcquat - A3c)
                         ]:
    if product == 0:
        print label + '= 0'
print 'Symmetry'
for (label, product) in [('i_{cq} K_1 i_{cq} + K_4 = ', icquat*K1*icquat + K4),
                         ('i_{cq} K_2 i_{cq} + K_3 = ', icquat*K2*icquat + K3),
                         ('i_{cq} K_3 i_{cq} + K_2 = ', icquat*K3*icquat + K2),
                         ('i_{cq} K_4 i_{cq} + K_1 = ', icquat*K4*icquat + K1)
                         ]:
    if product == 0:
        print label + '= 0'
print 'Symmetry'
for (label, product) in [('j_{cq} K_1 j_{cq} - K_3 = ', jcquat*K1*jcquat - K3),
                         ('j_{cq} K_2 j_{cq} - K_4 = ', jcquat*K2*jcquat - K4),
                         ('j_{cq} K_3 j_{cq} - K_1 = ', jcquat*K3*jcquat - K1),
                         ('j_{cq} K_4 j_{cq} - K_2 = ', jcquat*K4*jcquat - K2)
                         ]:
    if product == 0:
        print label + '= 0'
print 'Symmetry'
for (label, product) in [('k_{cq} K_1 k_{cq} - K_2 = ', kcquat*K1*kcquat - K2),
                         ('k_{cq} K_2 k_{cq} - K_1 = ', kcquat*K2*kcquat - K1),
                         ('k_{cq} K_3 k_{cq} - K_4 = ', kcquat*K3*kcquat - K4),
                         ('k_{cq} K_4 k_{cq} - K_3 = ', kcquat*K4*kcquat - K3)
                         ]:
    if product == 0:
        print label + '= 0'
print 'Symmetry (tilde is the reverse)'
for (label, product) in [('\\tilde{K}_1 + K_3 = ', K1.rev() + K3),
                         ('\\tilde{K}_2 + K_4 = ', K2.rev() + K4),
                         ('\\tilde{K}_3 + K_1 = ', K3.rev() + K1),
                         ('\\tilde{K}_4 + K_2 = ', K4.rev() + K2)
                         ]:
    if product == 0:
        print label + '= 0'
#print '**********************************************************************'
#for (label, product) in [('K_1 + (K_2 - 2 E i_{cq} i) = ', K1 + (K2 - 2*E*imag*icquat)),
#                         ('K_2 + (K_1 - 2 E i_{cq} i) = ', K2 + (K1 - 2*E*imag*icquat)),
#                         ('K_3 + (K_4 - 2 E i_{cq} i) = ', K3 + (K4 - 2*E*imag*icquat)),
#                         ('K_4 + (K_3 - 2 E i_{cq} i) = ', K4 + (K3 - 2*E*imag*icquat))
#                         ]:
#    if product == 0:
#        print label + '= 0'

print('Details about mixed constants (built from simple constants and the coquaternions)')
for (A,label) in [(Aa, [A1alabel, A2alabel, A3alabel, A4alabel]),
                  (Ab, [A1blabel, A2blabel, A3blabel, A4blabel]),
                  (Ac, [A1clabel, A2clabel, A3clabel, A4clabel]),
                  (Ad, [A1dlabel, A2dlabel, A3dlabel, A4dlabel])]:
    i = 1
    for (a, lab, k) in zip(A,label, K):
        a.convert_from_blades()
        print lab
        print lab + '=', a
        product=a*k
        product.expand()
        print 'Product\\ by\\ K_'+str(i)+'=', product
        print '**********************************************************************'
        i += 1

if outputTex:
    tex.xdvi(filename='evq.tex', debug=True)

#res=S(1)/8*(A1a*A3a + A1b*A3b + A1c*A3c + A1d*A3a)
#res=S(1)/8*(A1a*jcquat*A1d*jcquat + A1b*jcquat*A1c*jcquat + A1c*jcquat*A1b*jcquat + A1d*jcquat*A1a*jcquat)
#res.expand()
#print res
