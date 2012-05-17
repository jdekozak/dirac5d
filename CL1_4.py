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
#Flat space, no metric, just signature
#All constants = 1
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
imag.texLabel='i'
#CoQuaternions
icquat=gamma_t*gamma_w
jcquat=gamma_t*gamma_x*gamma_y*gamma_z
kcquat=gamma_x*gamma_y*gamma_z*gamma_w
icquat.texLabel='\\mathbf{i}'
jcquat.texLabel='\\mathbf{j}'
kcquat.texLabel='\\mathbf{k}'
#Quaternions
iquat=gamma_y*gamma_z
jquat=gamma_z*gamma_x
kquat=gamma_x*gamma_y
iquat.texLabel='\\boldsymbol{\\mathit{i}}'
jquat.texLabel='\\boldsymbol{\\mathit{j}}'
kquat.texLabel='\\boldsymbol{\\mathit{k}}'
#PseudoScalar
#MV.I=jcquat*imag

def CheckProperties(i,j,k, ilabel, jlabel, klabel):
    i.convert_from_blades()
    j.convert_from_blades()
    k.convert_from_blades()
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
    print 'CoQuaternions $(1,'+icquat.texLabel+','+jcquat.texLabel+','+kcquat.texLabel+')$ http://en.wikipedia.org/wiki/Split-quaternion'
    CheckProperties(icquat,jcquat,kcquat,
                    icquat.texLabel,jcquat.texLabel,kcquat.texLabel)
    print
    print 'Quaternions $(1,'+iquat.texLabel+','+jquat.texLabel+','+kquat.texLabel+')$ http://en.wikipedia.org/wiki/Quaternion'
    CheckProperties(iquat,jquat,kquat,
                    iquat.texLabel,jquat.texLabel,kquat.texLabel)
    print
    print 'Imaginary unit $' +imag.texLabel+'$'
    print imag.texLabel+'=',imag
    print imag.texLabel+'^2=', imag*imag
    print

#idempotents with icquat, jcquat, kcquat and jcquat*imag
i=1
U={}
for (signs, symbols) in zip([( -1, +1, -1),( -1, -1, +1),( +1, -1, -1),( +1, +1, +1), ( +1, -1, +1),( +1, +1, -1),( -1, +1, +1),( -1, -1, -1)],
                            [('-','+','-'),('-','-','+'),('+','-','-'),('+','+','+'), ('+','-','+'),('+','+','-'),('-','+','+'),('-','-','-')]):
    U[i]=0.5*(1+jcquat*imag)*0.5*(1+signs[0]*icquat+signs[1]*jcquat+signs[2]*kcquat)
    U[i].texLabel='\\frac{1}{2}(1+'+jcquat.texLabel+imag.texLabel+')\\frac{1}{2}(1'+symbols[0]+icquat.texLabel+symbols[1]+jcquat.texLabel+symbols[2]+kcquat.texLabel+')'
    i+=1
    U[i]=0.5*(1-jcquat*imag)*0.5*(1+signs[0]*icquat+signs[1]*jcquat+signs[2]*kcquat)
    U[i].texLabel='\\frac{1}{2}(1-'+jcquat.texLabel+imag.texLabel+')\\frac{1}{2}(1'+symbols[0]+icquat.texLabel+symbols[1]+jcquat.texLabel+symbols[2]+kcquat.texLabel+')'
    i+=1

def Idempotents():
    print('Idempotents in $Cl_{1,4}(\\mathbb{R})$')
    for index in range(1,17):
        print('U_{'+str(index)+'} : '+U[index].texLabel+'\Rightarrow '+'U_{'+str(index)+'}U_{'+str(index)+'} - U_{'+str(index)+'} = ' + str(U[index]*U[index]-U[index]))

########################################################################
#PHYSICS
########################################################################
parms = make_symbols('q m E p_x p_y p_z E_x E_y E_z B_x B_y B_z J_t J_x J_y J_z')

r = [y, z, x]
rquat = [iquat, jquat, kquat]

pv =[p_x, p_y, p_z]
Ev =[E_x, E_y, E_z]
Bv =[B_x, B_y, B_z]

p = S(0)
El = S(0)
B = S(0)

for (dim, var) in zip(pv, rquat):
    p += var * dim
for (dim, var) in zip(Ev, rquat):
    El += var * dim
for (dim, var) in zip(Bv, rquat):
    B += var * dim

p.texLabel='\\mathbf{p}'
El.texLabel='\\mathbf{E}'
B.texLabel='\\mathbf{B}'

J=J_t*gamma_t + J_x*gamma_x + J_y*gamma_y + J_z*gamma_z
J.texLabel='\\mathbf{J}'

r.append(w)
#potentiels A
ADm=S(1)*(-imag*(E*t -p_x*x -p_y*y -p_z*z) -m*w )
AD=-imag*(E*t -p_x*x -p_y*y -p_z*z)
AM=(E_x*x+E_y*y+E_z*z)*gamma_t-B_y*gamma_x*z-B_z*gamma_y*x-B_x*gamma_z*y
A=AD+AM

#Solutions particulieres
i=1
K={}
for (signs,symbols) in zip([(+1, +1), (+1, -1), (-1, +1), (-1, -1)],
                           [('+','+'),('+','-'),('-','+'),('-','-')]):
    k=(-icquat*E+signs[0]*imag*m+signs[1]*kcquat*p)
    k.texLabel='(-'+icquat.texLabel+'E'+symbols[0]+imag.texLabel+'m'+symbols[1]+kcquat.texLabel+p.texLabel+')'
    K[i] = k
    i += 1
texLabel='-'+K[2].texLabel+kcquat.texLabel
K[2]=-K[2]*kcquat
K[2].texLabel=texLabel
texLabel='-'+K[3].texLabel+jcquat.texLabel
K[3]=-K[3]*jcquat
K[3].texLabel=texLabel
texLabel=K[4].texLabel+icquat.texLabel
K[4]=K[4]*icquat
K[4].texLabel=texLabel
texLabel=K[1].texLabel+jcquat.texLabel+imag.texLabel
K[5]=K[1]*jcquat*imag
K[5].texLabel=texLabel
texLabel='-'+'(-'+icquat.texLabel+'E-'+imag.texLabel+'m+'+kcquat.texLabel+p.texLabel+')'+imag.texLabel
K[6]=-(-icquat*E-imag*m+kcquat*p)*imag
K[6].texLabel=texLabel
texLabel='(-'+icquat.texLabel+'E+'+imag.texLabel+'m-'+kcquat.texLabel+p.texLabel+')'+icquat.texLabel+imag.texLabel
K[7]=(-icquat*E+imag*m-kcquat*p)*icquat*imag
K[7].texLabel=texLabel
texLabel='-'+'(-'+icquat.texLabel+'E-'+imag.texLabel+'m-'+kcquat.texLabel+p.texLabel+')'+kcquat.texLabel+imag.texLabel
K[8]=-(-icquat*E-imag*m-kcquat*p)*kcquat*imag
K[8].texLabel=texLabel


########################################################################
#MAIN DIRAC
########################################################################
print('ALGEBRA')
print('')
print('Algebra playground is Clifford $Cl_{1,4}(\\mathbb{R})$, hence signature is (+ - - - -). The five dimensions are $(t, x, y, z, w)$.')
# 1 scalar
# 5 vectors
# 10 bivectors
# 10 tri vectors
# 5 pseudo vectors = quadri vectors
# 1 pseudo scalar
print('It defines two sets of quaternions with one imaginary unit.')
CheckGlobals()
print('Gradient definition')
print('{\\nabla}=({\gamma}_t \\frac{\\partial}{\\partial t}+{\gamma}_x \\frac{\\partial}{\\partial x}+{\gamma_y} \\frac{\\partial}{\\partial y}+{\gamma_z} \\frac{\\partial}{\\partial z}+{\gamma_w} \\frac{\\partial}{\\partial w})')
print('http://www.mrao.cam.ac.uk/~clifford/ptIIIcourse/course99/handouts/hout07.ps.gz')
Idempotents()

print('\\end{equation*}\\newpage\\begin{equation*}')
print('PHYSICS')
print('The following symbols are defined :')
print('Energy $E \\in \\mathbb{R}$')
print('Mass $m \\in \\mathbb{R}$')
print('Momentum $'+p.texLabel+'$ is defined with $p_x, p_y, p_z \\in \\mathbb{R}$')
print p.texLabel+'=p_x'+iquat.texLabel+'+p_y'+jquat.texLabel+'+p_z'+kquat.texLabel
print p.texLabel+'=', p
print('Electric field $'+El.texLabel+'$ is defined with $E_x, E_y, E_z \\in \\mathbb{R}$')
print El.texLabel+'=E_x'+iquat.texLabel+'+E_y'+jquat.texLabel+'+E_z'+kquat.texLabel
print El.texLabel+'=', El
print('Magnetic field $'+B.texLabel+'$ is defined with $B_x, B_y, B_z \\in \\mathbb{R}$')
print B.texLabel+'=B_x'+iquat.texLabel+'+B_y'+jquat.texLabel+'+B_z'+kquat.texLabel
print B.texLabel+'=', B
print('Wave : $K$ is a constant and $f$ is a function of $(t, x, y, z, w)$')
print('{\psi}=Ke^f')

print('First derivative')
print('{\\nabla}{\psi}={\\nabla}(Ke^f)')
print('{\\nabla}{\psi}=({\\nabla}K)e^f+K{\\nabla}(e^f)')
print('{\\nabla}K=0')
print('{\\nabla}{\psi}=K({\\nabla}f)e^f')

print ('$f$ is the exponential function')
print A


print('\\end{equation*}\\newpage\\begin{equation*}')
print('DIRAC')
print('http://en.wikipedia.org/wiki/Dirac\_equation')
print('0=({\gamma}_0 \\frac{\\partial}{\\partial t}+{\gamma}_1 \\frac{\\partial}{\\partial x}+{\gamma}_2 \\frac{\\partial}{\\partial y}+{\gamma}_3 \\frac{\\partial}{\\partial z}+'+imag.texLabel+'m) {\psi}')

print('Substitution in Dirac equation')
print('-'+imag.texLabel+'m\psi = (-'+icquat.texLabel+'\\frac{\\partial}{\\partial t}-'+iquat.texLabel+kcquat.texLabel+'\\frac{\\partial}{\\partial x}-'+jquat.texLabel+kcquat.texLabel+'\\frac{\\partial}{\\partial y}-'+kquat.texLabel+kcquat.texLabel+'\\frac{\\partial}{\\partial z}) {\psi}')

#print('Factorization and multiplication')
#print('0 = '+jcquat.texLabel+'(-'+kcquat.texLabel+'\\frac{\\partial}{\\partial t}+'+iquat.texLabel+icquat.texLabel+'\\frac{\\partial}{\\partial x}+'+jquat.texLabel+icquat.texLabel+'\\frac{\\partial}{\\partial y}+'+kquat.texLabel+icquat.texLabel+'\\frac{\\partial}{\\partial z}-'+jcquat.texLabel+imag.texLabel+'m) {\psi}')

if outputTex:
    tex.Format('1 1 1 1')
else:
    MV.set_str_format(1)

print('With the above gradients, identify the Dirac algebra aka gamma matrices')
#print('http://en.wikipedia.org/wiki/Gamma_matrices#Normalisation')
gamma_0 = -gamma_t*gamma_w
gamma_1 = gamma_x*gamma_w
gamma_2 = gamma_y*gamma_w
gamma_3 = gamma_z*gamma_w
gamma_5 = imag*gamma_0*gamma_1*gamma_2*gamma_3
print '{\gamma}_0 = ', gamma_0
print '-'+icquat.texLabel+' = ', -icquat
print '{\gamma}_0^2 = ', gamma_0*gamma_0
print '{\gamma}_1 = ', gamma_1
print '-'+iquat.texLabel+kcquat.texLabel+' = ', -iquat*kcquat
print '{\gamma}_1^2 = ', gamma_1*gamma_1
print '{\gamma}_2 = ', gamma_2
print '-'+jquat.texLabel+kcquat.texLabel+' = ', -jquat*kcquat
print '{\gamma}_2^2 = ', gamma_2*gamma_2
print '{\gamma}_3 = ', gamma_3
print '-'+kquat.texLabel+kcquat.texLabel+' = ', -kquat*kcquat
print '{\gamma}_3^2 = ', gamma_3*gamma_3
print '{\gamma}_5 = ', gamma_5
print '-'+imag.texLabel+jcquat.texLabel+' = ', -imag*jcquat
print '{\gamma}_5^2 = ', gamma_5*gamma_5

print('Exponential function with energy and momentum only, electric and magnetic fields are null $f$')
print 'f = ', AD
print('Gradient for $f$')
print '{\\nabla}f = ', AD.grad()
print('Square of the gradient is a Lorentz invariant')
print '{\\nabla}f^2 = ', AD.grad()*AD.grad()

print('Particular solutions $K$')
for (kkey,k) in K.items():
    klabel = 'K_'+str(kkey)
    print klabel+'='+k.texLabel
print('Check first derivative is null $K({\\nabla}f + i m)$')
for (kkey,k) in K.items():
    klabel = 'K_'+str(kkey)
    print klabel+'({\\nabla}f+ i m) =', k * (AD.grad()+imag*m)


print('\\end{equation*}\\newpage\\begin{equation*}')
print('MAXWELL')
print('Exponential function with electric and magnetic fields only $f$')
print 'f = ', AM
print('Gradient for $f$')
print '{\\nabla}f = ', AM.grad()
print('Square of the gradient is a Lorentz invariant')
print '{\\nabla}f^2 = ', AM.grad()*AM.grad()
#print ('${\\nabla}F_{Maxwell}$ is the current $J$')
#print 'J = ', F_Maxwell.grad()

print('\\end{equation*}\\newpage\\begin{equation*}')
print('DIRAC MAXWELL')
print ('${\\nabla}f$ is the electromagnetic field $F$ with a Dirac component')
if outputTex:
    tex.Format('1 1 1 3')
else:
    MV.set_str_format(1)
F_DiracMaxwell = A.grad()
print 'F_{DiracMaxwell} = ', F_DiracMaxwell

F_Maxwell = B-jcquat*El
F_Maxwell.texLabel = B.texLabel + '-' + jcquat.texLabel + El.texLabel 
print 'F_{Maxwell} = ' + F_Maxwell.texLabel 
print 'F_{Maxwell} = ', F_Maxwell

F_Dirac = AD.grad()
F_Dirac.texLabel = '('+kcquat.texLabel+p.texLabel+'-'+icquat.texLabel+'E)'
print 'F_{Dirac} = '+F_Dirac.texLabel
print 'F_{Dirac} = ', F_Dirac

#print('\\end{equation*}\\newpage\\begin{equation*}')


if outputTex:
    tex.xdvi(filename='CL1_4.tex', debug=True)
