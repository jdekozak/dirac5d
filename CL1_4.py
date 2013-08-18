# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, './sympy')

#computer algebra system
from sympy import *
from sympy.matrices import *
from sympy.galgebra.ga import *
from sympy.galgebra.printing import *

outputTex=True
#outputTex=False

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
vars = t, x, y, z, w = symbols('t x y z w')
myBasis='gamma_t gamma_x gamma_y gamma_z gamma_w'
gamma_t, gamma_x, gamma_y, gamma_z, gamma_w, grad = MV.setup(myBasis,metric,coords=vars)
if outputTex:
    Format()
#Imaginary unit
imag=gamma_w
imag.texLabel='i'
#CoQuaternions
icquat=gamma_t*gamma_w
jcquat=gamma_t*gamma_x*gamma_y*gamma_z*gamma_w
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
#MV.I=jcquat

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
    print '#CoQuaternions $('+imag.texLabel+','+icquat.texLabel+','+jcquat.texLabel+','+kcquat.texLabel+')$ http://en.wikipedia.org/wiki/Split-quaternion ?\\newline'
    CheckProperties(icquat,jcquat,kcquat,
                    icquat.texLabel,jcquat.texLabel,kcquat.texLabel)
    print
    print '#Quaternions $(1,'+iquat.texLabel+','+jquat.texLabel+','+kquat.texLabel+')$ http://en.wikipedia.org/wiki/Quaternion\\newline'
    CheckProperties(iquat,jquat,kquat,
                    iquat.texLabel,jquat.texLabel,kquat.texLabel)
    print
    print '#Imaginary unit $' +imag.texLabel+'$\\newline'
    print imag.texLabel+'=',imag
    print imag.texLabel+'^2=', imag*imag
    print

#idempotents with icquat, jcquat, kcquat and jcquat*imag
i=1
U={}
for (signs, the_symbols) in zip([( -1, +1, -1),( -1, -1, +1),( +1, -1, -1),( +1, +1, +1), ( +1, -1, +1),( +1, +1, -1),( -1, +1, +1),( -1, -1, -1)],
                            [('-','+','-'),('-','-','+'),('+','-','-'),('+','+','+'), ('+','-','+'),('+','+','-'),('-','+','+'),('-','-','-')]):
    U[i]=(1+jcquat)*(1+signs[0]*icquat+signs[1]*imag+signs[2]*kcquat)*(S(1)/4)
    U[i].texLabel='\\frac{1}{2}(1+'+jcquat.texLabel+')\\frac{1}{2}(1'+the_symbols[0]+icquat.texLabel+the_symbols[1]+imag.texLabel+the_symbols[2]+kcquat.texLabel+')'
    i+=1
    U[i]=(1-jcquat)*(1+signs[0]*icquat+signs[1]*imag+signs[2]*kcquat)*(S(1)/4)
    U[i].texLabel='\\frac{1}{2}(1-'+jcquat.texLabel+')\\frac{1}{2}(1'+the_symbols[0]+icquat.texLabel+the_symbols[1]+imag.texLabel+the_symbols[2]+kcquat.texLabel+')'
    i+=1

def Idempotents():
    print('#Idempotents in $Cl_{1,4}(\\mathbb{R})$')
    for index in range(1,17):
        print('U_{'+str(index)+'} : '+U[index].texLabel+'\Rightarrow '+'U_{'+str(index)+'}U_{'+str(index)+'} - U_{'+str(index)+'} = ' + str(U[index]*U[index]-U[index]))

########################################################################
#PHYSICS
########################################################################
q, m, E, p_x, p_y, p_z, E_x, E_y, E_z, B_x, B_y, B_z, phi, A_x, A_y, A_z = symbols('q m E p_x p_y p_z E_x E_y E_z B_x B_y B_z phi A_x A_y A_z')

R = x*gamma_x+y*gamma_y+z*gamma_z
rquat = [iquat, jquat, kquat]

pv =[p_x, p_y, p_z]
Ev =[E_x, E_y, E_z]
Bv =[B_x, B_y, B_z]
Av =[A_x, A_y, A_z]

p = S(0)
El= S(0)
B = S(0)
A = S(0)

for (dim, var) in zip(pv, rquat):
    p += var * dim
for (dim, var) in zip(Ev, rquat):
    El += var * dim
for (dim, var) in zip(Bv, rquat):
    B += var * dim
for (dim, var) in zip(Av, rquat):
    A += var * dim

p.texLabel='\\mathbf{p}'
El.texLabel='\\mathbf{E}'
B.texLabel='\\mathbf{B}'
A.texLabel='\\mathbf{A}'

#potentiels A
ADm=S(1)*(-imag*(E*t -p_x*x -p_y*y -p_z*z) -m*w )
AD=-imag*(E*t -p_x*x -p_y*y -p_z*z)
AM=(E_x*x+E_y*y+E_z*z)*gamma_t+(-B_y*z)*gamma_x+(-B_z*x)*gamma_y+(-B_x*y)*gamma_z
ADM=AD+AM

#Solutions particulieres
i=1
K={}
for (signs,the_symbols) in zip([(+1, +1), (+1, -1), (+1, +1), (-1, -1)],
                           [('+','+'),('+','-'),('+','+'),('-','-')]):
    k=(-icquat*E+signs[0]*imag*m+signs[1]*kcquat*p)
    k.texLabel='(-'+icquat.texLabel+'E'+the_symbols[0]+imag.texLabel+'m'+the_symbols[1]+kcquat.texLabel+p.texLabel+')'
    K[i] = k
    i += 1
texLabel='-'+K[2].texLabel+kcquat.texLabel
K[2]=-K[2]*kcquat
K[2].texLabel=texLabel
texLabel='-'+K[3].texLabel+jcquat.texLabel
K[3]=K[3]*jcquat
K[3].texLabel=texLabel
texLabel=K[4].texLabel+icquat.texLabel
K[4]=K[4]*icquat
K[4].texLabel=texLabel
texLabel='(-'+icquat.texLabel+'E-'+imag.texLabel+'m+'+kcquat.texLabel+p.texLabel+')'+jcquat.texLabel+imag.texLabel
K[5]=(-icquat*E-imag*m+kcquat*p)*jcquat*imag
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

def energy(self):
    return self.subs({E**2/4-m**2/4-p_x**2/4-p_y**2/4-p_z**2/4:0,-E**2/4+m**2/4+p_x**2/4+p_y**2/4+p_z**2/4:0})
MV.energy=energy

########################################################################
#MAIN DIRAC
########################################################################
print('#ALGEBRA\\newline')
print('#\\newline')
print('#Algebra playground is Clifford $Cl_{1,4}(\\mathbb{R})$, hence signature is $(+ - - - -)$. The five dimensions are $(t, x, y, z, w)$. It defines two sets of quaternions with one imaginary unit.\\newline')
# 1 scalar
# 5 vectors
# 10 bivectors
# 10 tri vectors
# 5 pseudo vectors = quadri vectors
# 1 pseudo scalar
CheckGlobals()
print('#Gradient definition\\newline')
print('{\\nabla}=({\gamma}_t \\frac{\\partial}{\\partial t}+{\gamma}_x \\frac{\\partial}{\\partial x}+{\gamma_y} \\frac{\\partial}{\\partial y}+{\gamma_z} \\frac{\\partial}{\\partial z}+{\gamma_w} \\frac{\\partial}{\\partial w})')
print('#http://www.mrao.cam.ac.uk/~clifford/ptIIIcourse/course99/handouts/hout07.ps.gz\\newline')
Idempotents()

print('\\end{equation*}\\newpage\\begin{equation*}')
print('#PHYSICS\\newline')
print('#The following symbols are defined :\\newline')
print('#Energy $E \\in \\mathbb{R}$\\newline')
print('#Mass $m \\in \\mathbb{R}$\\newline')
print('#Momentum $'+p.texLabel+'$ is defined with $p_x, p_y, p_z \\in \\mathbb{R}$\\newline')
print p.texLabel+'=p_x'+iquat.texLabel+'+p_y'+jquat.texLabel+'+p_z'+kquat.texLabel
print p.texLabel+'=', p
print('#Electric field $'+El.texLabel+'$ is defined with $E_x, E_y, E_z \\in \\mathbb{R}$\\newline')
print El.texLabel+'=E_x'+iquat.texLabel+'+E_y'+jquat.texLabel+'+E_z'+kquat.texLabel
print El.texLabel+'=', El
print('#Magnetic field $'+B.texLabel+'$ is defined with $B_x, B_y, B_z \\in \\mathbb{R}$\\newline')
print B.texLabel+'=B_x'+iquat.texLabel+'+B_y'+jquat.texLabel+'+B_z'+kquat.texLabel
print B.texLabel+'=', B
print('#Wave : $K$ is a constant and $f$ is a function of $(t, x, y, z, w)$\\newline')
print('{\psi}=Ke^f')

print('#First derivative\\newline')
print('{\\nabla}{\psi}={\\nabla}(Ke^f)')
print('{\\nabla}{\psi}=({\\nabla}K)e^f+K{\\nabla}(e^f)')
print('{\\nabla}K=0')
print('{\\nabla}{\psi_L}=({\\nabla}f)Ke^f')
print('{\\nabla}{\psi_R}=Ke^f({\\nabla}f)')
print('{\\nabla}{\psi_R}=K({\\nabla}f)e^f')
print('#Find solutions that fullfills left and right multiplication, (start with right multiplication)\\newline')

print ('#$f$ is the exponential function\\newline')
print ADM


print('\\end{equation*}\\newpage\\begin{equation*}')
print('#DIRAC\\newline')
print('#http://en.wikipedia.org/wiki/Dirac\_equation\\newline')
print('0=({\gamma}_0 \\frac{\\partial}{\\partial t}+{\gamma}_1 \\frac{\\partial}{\\partial x}+{\gamma}_2 \\frac{\\partial}{\\partial y}+{\gamma}_3 \\frac{\\partial}{\\partial z}+'+imag.texLabel+'m) {\psi}')

print('#Substitution in Dirac equation')
print('-'+imag.texLabel+'\\frac{\\partial}{\\partial w}\psi = (-'+icquat.texLabel+'\\frac{\\partial}{\\partial t}-'+iquat.texLabel+kcquat.texLabel+'\\frac{\\partial}{\\partial x}-'+jquat.texLabel+kcquat.texLabel+'\\frac{\\partial}{\\partial y}-'+kquat.texLabel+kcquat.texLabel+'\\frac{\\partial}{\\partial z}) {\psi}')
#print('Factorization and multiplication')
#print('0 = '+jcquat.texLabel+'(-'+kcquat.texLabel+'\\frac{\\partial}{\\partial t}+'+iquat.texLabel+icquat.texLabel+'\\frac{\\partial}{\\partial x}+'+jquat.texLabel+icquat.texLabel+'\\frac{\\partial}{\\partial y}+'+kquat.texLabel+icquat.texLabel+'\\frac{\\partial}{\\partial z}-'+jcquat.texLabel+imag.texLabel+'m) {\psi}')

print('#With the above gradients, identify the Dirac algebra aka gamma matrices\\newline')
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
print '-'+jcquat.texLabel+' = ', -jcquat
print '{\gamma}_5^2 = ', gamma_5*gamma_5

print('#Exponential function with energy and momentum only, electric and magnetic fields are null $f$\\newline')
print 'f = ', ADm
print('#Gradient for $f$\\newline')
print '{\\nabla}f = ', grad*ADm
print('#Square of the gradient is a Lorentz invariant\\newline')
print '{\\nabla}f^2 = ', (grad*ADm)*(grad*ADm)

print('#Particular solutions $K$\\newline')
for (kkey,k) in K.items():
    klabel = 'K_'+str(kkey)
    print klabel+'='+k.texLabel
print('#Check first derivative is null $K{\\nabla}f$\\newline')
for (kkey,k) in K.items():
    klabel = 'K_'+str(kkey)
    print klabel+'{\\nabla}f =', k * (grad*ADm)

print('#Right multiplication, left multiplication and combinations of the eight particular solutions\\newline')

print 'U_1 = ', U[1]

res=U[1]*K[1]*(grad*ADm)
res.expand()
print 'U_1K_1{\\nabla}f = ', res

res=(grad*ADm)*K[1]*U[1]
res.expand()
print '{\\nabla}f K_1U_1 = ', res

psi=-K[1]*U[1]*K[1]
psi.expand()
psi=psi.energy()

print('\psi_L = -e^{-f}K_1U_1')
print('\psi_R = U_1K_1e^f')
print('\psi = \psi_L\psi_R')
print('\psi : -e^{-f}K_1U_1U_1K_1e^f = -e^{-f}K_1U_1K_1e^f')
print('#IBOZOO UU : $-K_1U_1K_1 = \Psi$\\newline')
print '\Psi = ', psi
print('\\end{equation*}\\newpage\\begin{equation*}')
print('#Dirac observables handout 11 chapter 3.1\\newline')
i=1
O={}
#for (aBasis,aBasisLabel) in [(1,''),(imag,imag.texLabel),(icquat,icquat.texLabel),(kcquat,kcquat.texLabel)]:
#    O[i]=psi*aBasis*psi.rev()
#    O[i].texLabel = '{\psi}'+aBasisLabel+'{\\tilde{\psi}}'
#    O[i].expand()
#    print 'O_'+str(i)+'='+ O[i].texLabel
#    print 'O_'+str(i)+'=', O[i]
#    print('\\end{equation*}\\newpage\\begin{equation*}')
#    i+=1


print('\\end{equation*}\\newpage\\begin{equation*}')
print('#MAXWELL\\newline')
print('#Exponential function with electric and magnetic fields only $f$\\newline')
print 'f = ', AM
print('#Gradient for $f$\\newline')
print '{\\nabla}f = ', grad*AM
print('#Square of the gradient is a Lorentz invariant\\newline')
print '{\\nabla}f^2 = ', (grad*AM)*(grad*AM)
#print ('{\\nabla}F_{Maxwell} is the current J')
#print 'J = ', F_Maxwell.grad()

print('\\end{equation*}\\newpage\\begin{equation*}')
print('#DIRAC MAXWELL\\newline')
print ('#${\\nabla}f$ is the electromagnetic field $F$ with a Dirac component\\newline')

F_DiracMaxwell = grad*ADM
print 'F_{DiracMaxwell} = ', F_DiracMaxwell

F_Maxwell = -imag*(imag*B-jcquat*El)
F_Maxwell.texLabel = '-'+imag.texLabel+'('+imag.texLabel+B.texLabel + '-' + jcquat.texLabel + El.texLabel+')'
print 'F_{Maxwell} = ' + F_Maxwell.texLabel 
print 'F_{Maxwell} = ', F_Maxwell

F_Dirac = grad*AD
F_Dirac.texLabel = '('+kcquat.texLabel+p.texLabel+'-'+icquat.texLabel+'E)'
print 'F_{Dirac} = '+F_Dirac.texLabel
print 'F_{Dirac} = ', F_Dirac

print('\\end{equation*}\\newpage\\begin{equation*}')


print('#DIRAC COUPLING\\newline')
AD=-imag*(E*t -p_x*x -p_y*y -p_z*z)+q*w*(phi*gamma_t+A_x*gamma_x+A_y*gamma_y+A_z*gamma_z)
print('-('+imag.texLabel+'m+'+icquat.texLabel+'q \phi+'+kcquat.texLabel+'q'+A.texLabel+')\psi = (-'+icquat.texLabel+'\\frac{\\partial}{\\partial t}-'+iquat.texLabel+kcquat.texLabel+'\\frac{\\partial}{\\partial x}-'+jquat.texLabel+kcquat.texLabel+'\\frac{\\partial}{\\partial y}-'+kquat.texLabel+kcquat.texLabel+'\\frac{\\partial}{\\partial z}) {\psi}')
print('#Exponential function with energy and momentum only, and a potential that depends on $w$\\newline')
print('#Electric and magnetic field are null\\newline')
print 'f = ', AD
print('#Gradient for $f$')
print '{\\nabla}f = ', grad*AD

ADM=AD+AM
print('#Full exponential function\\newline')
print 'f = ', ADM

print('#Gradient for $f$\\newline')
print '{\\nabla}f = ', grad*ADM

print('\\end{equation*}\\newpage\\begin{equation*}')

if outputTex:
    #xdvi(filename='CL1_4.tex', debug=True)
    xdvi()
