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
icquat=gamma_t
jcquat=gamma_t*gamma_x*gamma_y*gamma_z
kcquat=gamma_x*gamma_y*gamma_z
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

#Clear equality
def energy(self):
    return self.subs(-0.5*m**2-0.5*p_x**2-0.5*p_y**2-0.5*p_z**2,-0.5*E**2).subs(0.5*m**2+0.5*p_x**2+0.5*p_y**2+0.5*p_z**2,0.5*E**2).subs(-m**2-p_x**2-p_y**2-p_z**2,-E**2).subs(m**2+p_x**2+p_y**2+p_z**2,E**2)
MV.energy=energy

########################################################################
#PHYSICS
########################################################################
parms = make_symbols('m E p_x p_y p_z q A_x A_y A_z')
r = [x, y, z]
rquat = [iquat, jquat, kquat]
#vecteur d'onde (propagation spatiale)
p = [p_x, p_y, p_z]
A = [A_x, A_y, A_z]
psigma = S(0)
for (dim, var) in zip(p, r):
    psigma += var * dim
pmulti = S(0)
for (dim, var) in zip(p, rquat):
    pmulti += var * dim
pmulti.texLabel='\\mathbf{p}'
Asigma = S(0)
for (dim, var) in zip(A, r):
    Asigma += var * dim
Amulti = S(0)
for (dim, var) in zip(A, rquat):
    Amulti += var * dim
Amulti.texLabel='\\mathbf{A}'

r.append(w)
#vecteur d'onde (propagation spatiale) en milieu absorbant ou amplificateur
F=S(1)*(-imag*(E*t -p_x*x -p_y*y -p_z*z) -m*w )
#Solutions particulieres
i=1
K={}
for (signs,symbols) in zip([(-1, +1), (+1, -1), (+1, +1), (-1, -1)],
                           [('-','+'),('+','-'),('+','+'),('-','-')]):
    k=jcquat*imag*(kcquat*E+signs[0]*jcquat*m+signs[1]*icquat*pmulti)
    k.texLabel=jcquat.texLabel+imag.texLabel+'('+kcquat.texLabel+'E'+symbols[0]+jcquat.texLabel+'m'+symbols[1]+icquat.texLabel+pmulti.texLabel+')'
    K[i] = k
    i += 1
texLabel=K[2].texLabel+kcquat.texLabel
K[2]=K[2]*kcquat
K[2].texLabel=texLabel
texLabel='-'+K[3].texLabel+jcquat.texLabel
K[3]=-K[3]*jcquat
K[3].texLabel=texLabel
texLabel='-'+K[4].texLabel+icquat.texLabel
K[4]=-K[4]*icquat
K[4].texLabel=texLabel
texLabel=K[1].texLabel+jcquat.texLabel+imag.texLabel
K[5]=K[1]*jcquat*imag
K[5].texLabel=texLabel
texLabel='-'+jcquat.texLabel+imag.texLabel+'('+kcquat.texLabel+'E+'+jcquat.texLabel+'m+'+icquat.texLabel+pmulti.texLabel+')'+imag.texLabel
K[6]=-jcquat*imag*(kcquat*E+jcquat*m+icquat*pmulti)*imag
K[6].texLabel=texLabel
texLabel=jcquat.texLabel+imag.texLabel+'('+kcquat.texLabel+'E+'+jcquat.texLabel+'m-'+icquat.texLabel+pmulti.texLabel+')'+icquat.texLabel+imag.texLabel
K[7]=jcquat*imag*(kcquat*E+jcquat*m-icquat*pmulti)*icquat*imag
K[7].texLabel=texLabel
texLabel='-'+jcquat.texLabel+imag.texLabel+'('+kcquat.texLabel+'E-'+jcquat.texLabel+'m-'+icquat.texLabel+pmulti.texLabel+')'+kcquat.texLabel+imag.texLabel
K[8]=-jcquat*imag*(kcquat*E-jcquat*m-icquat*pmulti)*kcquat*imag
K[8].texLabel=texLabel

#idempotents with icquat, jcquat, kcquat and jcquat*imag
i=1
idem={}
for (signs, symbols) in zip([( -1, +1, -1),( -1, -1, +1),( +1, -1, -1),( +1, +1, +1), ( +1, -1, +1),( +1, +1, -1),( -1, +1, +1),( -1, -1, -1)],
                            [('-','+','-'),('-','-','+'),('+','-','-'),('+','+','+'), ('+','-','+'),('+','+','-'),('-','+','+'),('-','-','-')]):
    idem[i]=0.5*(1+jcquat*imag)*0.5*(1+signs[0]*icquat+signs[1]*jcquat+signs[2]*kcquat)
    idem[i].texLabel='\\frac{1}{2}(1+'+jcquat.texLabel+imag.texLabel+')\\frac{1}{2}(1'+symbols[0]+icquat.texLabel+symbols[1]+jcquat.texLabel+symbols[2]+kcquat.texLabel+')'
    i+=1
    idem[i]=0.5*(1-jcquat*imag)*0.5*(1+signs[0]*icquat+signs[1]*jcquat+signs[2]*kcquat)
    idem[i].texLabel='\\frac{1}{2}(1-'+jcquat.texLabel+imag.texLabel+')\\frac{1}{2}(1'+symbols[0]+icquat.texLabel+symbols[1]+jcquat.texLabel+symbols[2]+kcquat.texLabel+')'
    i+=1

def Idempotents():
    print('Some idempotents in $Cl_{1,4}(\\mathbb{R})$')
    for index in range(1,17):
        print('Id_{'+str(index)+'} = '+idem[index].texLabel)
        print('Id_{'+str(index)+'}Id_{'+str(index)+'} - Id_{'+str(index)+'} = ' + str(idem[index]*idem[index]-idem[index]))

########################################################################
#MAIN DIRAC
########################################################################
print('Algebra is $Cl_{1,4}(\\mathbb{R})$')
print('The five dimensions are $(t, x, y, z, w)$.')
print('It defines two sets of quaternions with one imaginary unit.')
print('Signature is (+ - - - -)')
CheckGlobals()
print('Gradient definition')
print('{\\nabla}=({\gamma}_t \\frac{\\partial}{\\partial t}+{\gamma}_x \\frac{\\partial}{\\partial x}+{\gamma_y} \\frac{\\partial}{\\partial y}+{\gamma_z} \\frac{\\partial}{\\partial z}+{\gamma_w} \\frac{\\partial}{\\partial w})')
print('http://www.mrao.cam.ac.uk/~clifford/ptIIIcourse/course99/handouts/hout07.ps.gz')

print('Wave : $K$ is a constant and $f$ is a function of $(t, x, y, z, w)$')
print('{\psi}=Ke^f')

print('First derivative')
print('{\\nabla}{\psi}={\\nabla}(Ke^f)')
print('{\\nabla}{\psi}=({\\nabla}K)e^f+K{\\nabla}(e^f)')
print('{\\nabla}K=0')
print('{\\nabla}{\psi}=K({\\nabla}f)e^f')

print('The following symbols are defined :')
print('$E \\in \\mathbb{R}$')
print('$m \\in \\mathbb{R}$')
print('$'+pmulti.texLabel+'$ is defined with $p_x, p_y, p_z \\in \\mathbb{R}$')
print pmulti.texLabel+'=p_x'+iquat.texLabel+'+p_y'+jquat.texLabel+'+p_z'+kquat.texLabel
print pmulti.texLabel+'=', pmulti

print('Exponential function $f$')
print 'f = ', F
print('Gradient for $f$')
print '{\\nabla}f = ', F.grad()
print('Square of the gradient')
print '{\\nabla}f^2 = ', F.grad()*F.grad()

print('Dirac')
print('http://en.wikipedia.org/wiki/Dirac\_equation')
print('0=({\gamma}_0 \\frac{\\partial}{\\partial t}+{\gamma}_1 \\frac{\\partial}{\\partial x}+{\gamma}_2 \\frac{\\partial}{\\partial y}+{\gamma}_3 \\frac{\\partial}{\\partial z}+'+imag.texLabel+'m) {\psi}')

print('With the above gradients, identify the Dirac algebra aka gamma matrices')
#print('http://en.wikipedia.org/wiki/Gamma_matrices#Normalisation')
gamma_0 =  gamma_t*gamma_w
gamma_1 = -gamma_x*gamma_w
gamma_2 = -gamma_y*gamma_w
gamma_3 = -gamma_z*gamma_w
gamma_5 = imag*gamma_0*gamma_1*gamma_2*gamma_3
print '{\gamma}_0 = ', gamma_0
print '-'+imag.texLabel+icquat.texLabel+' = ', -imag*icquat
print '{\gamma}_0^2 = ', gamma_0*gamma_0
print '{\gamma}_1 = ', gamma_1
print '-'+imag.texLabel+iquat.texLabel+kcquat.texLabel+' = ', -imag*iquat*kcquat
print '{\gamma}_1^2 = ', gamma_1*gamma_1
print '{\gamma}_2 = ', gamma_2
print '-'+imag.texLabel+jquat.texLabel+kcquat.texLabel+' = ', -imag*jquat*kcquat
print '{\gamma}_2^2 = ', gamma_2*gamma_2
print '{\gamma}_3 = ', gamma_3
print '-'+imag.texLabel+kquat.texLabel+kcquat.texLabel+' = ', -imag*kquat*kcquat
print '{\gamma}_3^2 = ', gamma_3*gamma_3
print '{\gamma}_5 = ', gamma_5
print '-'+imag.texLabel+jcquat.texLabel+' = ', -imag*jcquat
print '{\gamma}_5^2 = ', gamma_5*gamma_5

print('Substitution in Dirac equation')
print('0 = (-'+imag.texLabel+icquat.texLabel+'\\frac{\\partial}{\\partial t}-'+imag.texLabel+iquat.texLabel+kcquat.texLabel+'\\frac{\\partial}{\\partial x}-'+imag.texLabel+jquat.texLabel+kcquat.texLabel+'\\frac{\\partial}{\\partial y}-'+imag.texLabel+kquat.texLabel+kcquat.texLabel+'\\frac{\\partial}{\\partial z}+'+imag.texLabel+'m) {\psi}')

print('\\end{equation*}\\newpage\\begin{equation*}')
print('Factorization and multiplication')
print('0 = '+jcquat.texLabel+imag.texLabel+'('+kcquat.texLabel+'\\frac{\\partial}{\\partial t}+'+iquat.texLabel+icquat.texLabel+'\\frac{\\partial}{\\partial x}+'+jquat.texLabel+icquat.texLabel+'\\frac{\\partial}{\\partial y}+'+kquat.texLabel+icquat.texLabel+'\\frac{\\partial}{\\partial z}-'+jcquat.texLabel+'m) {\psi}')

print('Particular solutions $K$')
for (kkey,k) in K.items():
    klabel = 'K_'+str(kkey)
    print klabel+'='+k.texLabel
#    k.convert_from_blades()
#    print klabel+'=', k
print('Check first derivative is null $K{\\nabla}f$')
for (kkey,k) in K.items():
    klabel = 'K_'+str(kkey)
    print klabel+'{\\nabla}f=', k * F.grad()
print('8 previous solutions can be combined linearly together and factorized, this gives 16 idempotents !')
Idempotents()

#if outputTex:
#    tex.Format('1 1 1 2')
#else:
#    MV.set_str_format(1)


if outputTex:
    tex.xdvi(filename='CL1_4.tex', debug=True)
