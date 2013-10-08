# -*- coding: utf-8 -*-
from sp5d import *

outputTex=True
#outputTex=False

if outputTex:
    Format()

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


def CheckProperties(i,j,k, ilabel, jlabel, klabel):
    i.Fmt(fmt=1, title=ilabel)
    j.Fmt(fmt=1, title=jlabel)
    k.Fmt(fmt=1, title=klabel)
    (i*i).Fmt(fmt=1, title='%'+ilabel+'^2')
    (j*j).Fmt(fmt=1, title='%'+jlabel+'^2')
    (k*k).Fmt(fmt=1, title='%'+klabel+'^2')
    (i*j).Fmt(fmt=1, title=ilabel + jlabel)
    (j*i).Fmt(fmt=1, title=jlabel + ilabel)
    (j*k).Fmt(fmt=1, title=jlabel + klabel)
    (k*j).Fmt(fmt=1, title=klabel + jlabel)
    (k*i).Fmt(fmt=1, title=klabel + ilabel)
    (i*k).Fmt(fmt=1, title=ilabel + klabel)
    (i*j*k).Fmt(fmt=1, title=ilabel + jlabel + klabel)

def CheckGlobals():
    '''Check algebra definitions'''
    print '#Associative Hyperbolic Quaternions $(1,'+ihquat.texLabel+','+jhquat.texLabel+','+khquat.texLabel+')$\\newline'
    CheckProperties(ihquat,jhquat,khquat,
                    ihquat.texLabel,jhquat.texLabel,khquat.texLabel)
    print
    print '#Quaternions $(1,'+iquat.texLabel+','+jquat.texLabel+','+kquat.texLabel+')$ http://en.wikipedia.org/wiki/Quaternion\\newline'
    CheckProperties(iquat,jquat,kquat,
                    iquat.texLabel,jquat.texLabel,kquat.texLabel)
    print
    print '#Imaginary unit $' +imag.texLabel+'$\\newline'
    imag.Fmt(fmt=1, title=imag.texLabel)
    (imag*imag).Fmt(fmt=1, title='%'+imag.texLabel+'^2')
    print

#idempotents with ihquat, jhquat, khquat and jhquat*imag
i=1
U={}
for (signs, the_symbols) in zip([( -1, +1, -1),( -1, -1, +1),( +1, -1, -1),( +1, +1, +1), ( +1, -1, +1),( +1, +1, -1),( -1, +1, +1),( -1, -1, -1)],
                            [('-','+','-'),('-','-','+'),('+','-','-'),('+','+','+'), ('+','-','+'),('+','+','-'),('-','+','+'),('-','-','-')]):
    U[i]=(1+jhquat)*(1+signs[0]*ihquat*imag+signs[1]*imag+signs[2]*khquat)*(S(1)/4)
    U[i].texLabel='\\frac{1}{2}(1+'+jhquat.texLabel+')\\frac{1}{2}(1'+the_symbols[0]+ihquat.texLabel+imag.texLabel+the_symbols[1]+imag.texLabel+the_symbols[2]+khquat.texLabel+')'
    i+=1
    U[i]=(1-jhquat)*(1+signs[0]*ihquat*imag+signs[1]*imag+signs[2]*khquat)*(S(1)/4)
    U[i].texLabel='\\frac{1}{2}(1-'+jhquat.texLabel+')\\frac{1}{2}(1'+the_symbols[0]+ihquat.texLabel+imag.texLabel+the_symbols[1]+imag.texLabel+the_symbols[2]+khquat.texLabel+')'
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
AD=       -imag*(E*t -p_x*x -p_y*y -p_z*z)
AM=(E_x*x+E_y*y+E_z*z)*gamma_t+(-B_y*z)*gamma_x+(-B_z*x)*gamma_y+(-B_x*y)*gamma_z
ADM=ADm+AM

#Solutions particulieres
i=1
K={}
for (signs,the_symbols) in zip([(+1, +1), (+1, -1), (+1, +1), (-1, -1)],
                           [('+','+'),('+','-'),('+','+'),('-','-')]):
    k=(-ihquat*imag*E+signs[0]*imag*m+signs[1]*khquat*p)
    k.texLabel='(-'+ihquat.texLabel+imag.texLabel+'E'+the_symbols[0]+imag.texLabel+'m'+the_symbols[1]+khquat.texLabel+p.texLabel+')'
    K[i] = k
    i += 1
texLabel='-'+K[2].texLabel+khquat.texLabel
K[2]=-K[2]*khquat
K[2].texLabel=texLabel
texLabel=K[3].texLabel+jhquat.texLabel
K[3]=K[3]*jhquat
K[3].texLabel=texLabel
texLabel=K[4].texLabel+ihquat.texLabel+imag.texLabel
K[4]=K[4]*ihquat*imag
K[4].texLabel=texLabel
texLabel='(-'+ihquat.texLabel+imag.texLabel+'E-'+imag.texLabel+'m+'+khquat.texLabel+p.texLabel+')'+jhquat.texLabel+imag.texLabel
K[5]=(-ihquat*imag*E-imag*m+khquat*p)*jhquat*imag
K[5].texLabel=texLabel
texLabel='-'+'(-'+ihquat.texLabel+imag.texLabel+'E-'+imag.texLabel+'m+'+khquat.texLabel+p.texLabel+')'+imag.texLabel
K[6]=-(-ihquat*imag*E-imag*m+khquat*p)*imag
K[6].texLabel=texLabel
texLabel='(-'+ihquat.texLabel+imag.texLabel+'E+'+imag.texLabel+'m-'+khquat.texLabel+p.texLabel+')'+ihquat.texLabel
K[7]=(-ihquat*imag*E+imag*m-khquat*p)*ihquat
K[7].texLabel=texLabel
texLabel='-'+'(-'+ihquat.texLabel+imag.texLabel+'E-'+imag.texLabel+'m-'+khquat.texLabel+p.texLabel+')'+khquat.texLabel+imag.texLabel
K[8]=-(-ihquat*imag*E-imag*m-khquat*p)*khquat*imag
K[8].texLabel=texLabel

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
p.Fmt(fmt=1, title=p.texLabel)
print('#Electric field $'+El.texLabel+'$ is defined with $E_x, E_y, E_z \\in \\mathbb{R}$\\newline')
print El.texLabel+'=E_x'+iquat.texLabel+'+E_y'+jquat.texLabel+'+E_z'+kquat.texLabel
El.Fmt(fmt=1, title=El.texLabel)
print('#Magnetic field $'+B.texLabel+'$ is defined with $B_x, B_y, B_z \\in \\mathbb{R}$\\newline')
print B.texLabel+'=B_x'+iquat.texLabel+'+B_y'+jquat.texLabel+'+B_z'+kquat.texLabel
B.Fmt(fmt=1, title=B.texLabel)
print('#Wavefunction : $K$ is a constant and $f$ is a function of $(t, x, y, z, w)$\\newline')
print('{\psi}=Ke^f')

print('#First derivative\\newline')
print('{\\nabla}{\psi}={\\nabla}(Ke^f)')
print('{\\nabla}{\psi}=({\\nabla}K)e^f+K{\\nabla}(e^f)')
print('{\\nabla}K=0')
print('{\\nabla}{\psi_L}=({\\nabla}f)Ke^f')
print('{\\nabla}{\psi_R}=Ke^f({\\nabla}f)')
print('{\\nabla}{\psi_R}=K({\\nabla}f)e^f')
print('K({\\nabla}f) = ({\\nabla}f)K = 0')
print('#Find solutions that fullfills left and right multiplication\\newline')

print ('#$f$ is the exponential function\\newline')
ADM.Fmt(fmt=1,title='f')


print('\\end{equation*}\\newpage\\begin{equation*}')
print('#DIRAC\\newline')
print('#http://en.wikipedia.org/wiki/Dirac\_equation\\newline')
print('0=({\gamma}_0 \\frac{\\partial}{\\partial t}+{\gamma}_1 \\frac{\\partial}{\\partial x}+{\gamma}_2 \\frac{\\partial}{\\partial y}+{\gamma}_3 \\frac{\\partial}{\\partial z}+'+imag.texLabel+'m) {\psi}')

#print('Factorization and multiplication')
#print('0 = '+jhquat.texLabel+'(-'+khquat.texLabel+'\\frac{\\partial}{\\partial t}+'+iquat.texLabel+ihquat.texLabel+'\\frac{\\partial}{\\partial x}+'+jquat.texLabel+ihquat.texLabel+'\\frac{\\partial}{\\partial y}+'+kquat.texLabel+ihquat.texLabel+'\\frac{\\partial}{\\partial z}-'+jhquat.texLabel+imag.texLabel+'m) {\psi}')

print('#Define the Dirac algebra (aka gamma matrices) with the two sets of quaternions\\newline')
#print('http://en.wikipedia.org/wiki/Gamma_matrices#Normalisation')
gamma_0 = -gamma_t*gamma_w
gamma_1 = gamma_x*gamma_w
gamma_2 = gamma_y*gamma_w
gamma_3 = gamma_z*gamma_w
gamma_5 = imag*gamma_0*gamma_1*gamma_2*gamma_3
print '{\gamma}_0 = ', gamma_0
print '-'+ihquat.texLabel+imag.texLabel+' = ', -ihquat*imag
print '%{\gamma}_0^2 = ', gamma_0*gamma_0
print '{\gamma}_1 = ', gamma_1
print '-'+iquat.texLabel+khquat.texLabel+' = ', -iquat*khquat
print '%{\gamma}_1^2 = ', gamma_1*gamma_1
print '{\gamma}_2 = ', gamma_2
print '-'+jquat.texLabel+khquat.texLabel+' = ', -jquat*khquat
print '%{\gamma}_2^2 = ', gamma_2*gamma_2
print '{\gamma}_3 = ', gamma_3
print '-'+kquat.texLabel+khquat.texLabel+' = ', -kquat*khquat
print '%{\gamma}_3^2 = ', gamma_3*gamma_3
print '{\gamma}_5 = ', gamma_5
print '-'+jhquat.texLabel+' = ', -jhquat
print '%{\gamma}_5^2 = ', gamma_5*gamma_5

print('#Substitution in Dirac equation')
print('-'+imag.texLabel+'\\frac{\\partial}{\\partial w}\psi = (-'+ihquat.texLabel+imag.texLabel+'\\frac{\\partial}{\\partial t}-'+iquat.texLabel+khquat.texLabel+'\\frac{\\partial}{\\partial x}-'+jquat.texLabel+khquat.texLabel+'\\frac{\\partial}{\\partial y}-'+kquat.texLabel+khquat.texLabel+'\\frac{\\partial}{\\partial z}) {\psi}')

print('#Exponential function with energy and momentum only, electric and magnetic fields are null $f$\\newline')
ADm.Fmt(fmt=1, title='f')
print('#Gradient for $f$\\newline')
(grad*ADm).Fmt(fmt=1, title='{\\nabla}f')
print('#Square of the gradient is a Lorentz invariant\\newline')
((grad*ADm)*(grad*ADm)).Fmt(fmt=1, title='%{\\nabla}f^2')

print('#Particular solutions $K$\\newline')
for (kkey,k) in K.items():
    klabel = 'K_'+str(kkey)
    print klabel+'='+k.texLabel 
print('#Check first derivative is null $K{\\nabla}f$\\newline')
for (kkey,k) in K.items():
    klabel = 'K_'+str(kkey)
    (k * (grad*ADm)).Fmt(fmt=1, title=klabel+'{\\nabla}f')
print('#Check first derivative is null ${\\nabla}fK$\\newline')
for (kkey,k) in K.items():
    klabel = 'K_'+str(kkey)
    ((grad*ADm)*k).Fmt(fmt=1, title='{\\nabla}f'+klabel)

print('#Right multiplication, left multiplication\\newline')
print('')

U[1].Fmt(fmt=1, title='U_1')
(jhquat*U[1] - U[1]).Fmt(fmt=1, title=jhquat.texLabel+'U_1 - U_1')
(U[1]*jhquat - U[1]).Fmt(fmt=1, title='U_1'+jhquat.texLabel+' - U_1')
(K[3] - K[1]*jhquat).Fmt(fmt=1, title='K_3 - K_1'+jhquat.texLabel)
(K[3] - jhquat*K[1]).Fmt(fmt=1, title='K_3 - '+jhquat.texLabel+'K_1')

res=U[1]*K[1]*(grad*ADm)
#res.expand()

res.Fmt(fmt=2, title='U_1K_1{\\nabla}f')

res=(grad*ADm)*K[1]*U[1]
#res.expand()

res.Fmt(fmt=2, title='{\\nabla}f K_1U_1')

(K[1]*U[1] + U[1]*K[1]).Fmt(fmt=2, title='U_1K_1 + K_1U_1')

print('#16 idempotents, 16 solutions ?')
print('\psi = U_1K_1e^f')
print('\\end{equation*}\\newpage\\begin{equation*}')

print('#Dirac observables handout 11 chapter 3.1\\newline')
print('%\psi^{\dag} = -e^{-f}K_1U_1')
print('%\\rho = \psi^{\dag}\psi = -e^{-f}K_1U_1U_1K_1e^f')
print('%\\rho = \psi^{\dag}\psi = -e^{-f}K_1U_1K_1e^f')
print('#this is the definition of a rotation :')
print('%e^{-f}e^f=1')
print('-K_1U_1K_1 = \Phi')
phi=-K[1]*U[1]*K[1]
#phi.expand()

phi.Fmt(fmt=2, title='\Phi')
i=1
O={}
#for (aBasis,aBasisLabel) in [(1,''),(imag,imag.texLabel),(ihquat,ihquat.texLabel),(khquat,khquat.texLabel)]:
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
AM.Fmt(fmt=1, title='f')
print('#Gradient for $f$\\newline')
(grad*AM).Fmt(fmt=1, title='{\\nabla}f')
print('#Square of the gradient is a Lorentz invariant\\newline')
((grad*AM)*(grad*AM)).Fmt(fmt=1, title='%{\\nabla}f^2')

print('\\end{equation*}\\newpage\\begin{equation*}')
print('#DIRAC MAXWELL\\newline')
print ('#${\\nabla}f$ is the electromagnetic field $F$ with a Dirac component\\newline')
ADM.Fmt(fmt=1, title='f')

F_DiracMaxwell = grad*ADM
F_DiracMaxwell.Fmt(fmt=1, title='F_{DiracMaxwell} : {\\nabla}f')

F_Maxwell = -imag*(imag*B-jhquat*El)
F_Maxwell.texLabel = '-'+imag.texLabel+'('+imag.texLabel+B.texLabel + '-' + jhquat.texLabel + El.texLabel+')'
print 'F_{Maxwell} = ' + F_Maxwell.texLabel 
F_Maxwell.Fmt(fmt=1, title='F_{Maxwell}')

F_Dirac = grad*ADm
F_Dirac.texLabel = '('+khquat.texLabel+p.texLabel+'-'+ihquat.texLabel+imag.texLabel+'E)'
print 'F_{Dirac} = '+F_Dirac.texLabel
F_Dirac.Fmt(fmt=1, title='F_{Dirac}')

print('#Square of the gradient is a Lorentz invariant\\newline')
((grad*ADM)*(grad*ADM)).Fmt(fmt=2, title='%{\\nabla}f^2')


print('\\end{equation*}\\newpage\\begin{equation*}')


print('#DIRAC COUPLING\\newline')
AD=-imag*(E*t -p_x*x -p_y*y -p_z*z)+q*w*(phi*gamma_t+A_x*gamma_x+A_y*gamma_y+A_z*gamma_z)
print('-('+imag.texLabel+'m+'+ihquat.texLabel+imag.texLabel+'q \phi+'+khquat.texLabel+'q'+A.texLabel+')\psi = (-'+ihquat.texLabel+imag.texLabel+'\\frac{\\partial}{\\partial t}-'+iquat.texLabel+khquat.texLabel+'\\frac{\\partial}{\\partial x}-'+jquat.texLabel+khquat.texLabel+'\\frac{\\partial}{\\partial y}-'+kquat.texLabel+khquat.texLabel+'\\frac{\\partial}{\\partial z}) {\psi}')
print('#Exponential function with energy and momentum only, and a potential that depends on $w$\\newline')
print('#Electric and magnetic field are null\\newline')
AD.Fmt(fmt=1, title='f')
print('#Gradient for $f$')
(grad*AD).Fmt(fmt=1, title='{\\nabla}f')

ADM=AD+AM
print('#Full exponential function\\newline')
ADM.Fmt(fmt=1, title='f')

print('#Gradient for $f$\\newline')
(grad*ADM).Fmt(fmt=3, title='{\\nabla}f')

print('\\end{equation*}\\newpage\\begin{equation*}')

if outputTex:
    xpdf()