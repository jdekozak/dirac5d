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

#Computation with the gradient
def CheckComputedGrad(Grad):
    i=1
    for grad in Grad:
        print '{\\nabla}f_'+str(i)+'^2=',grad*grad
        i=i+1

def CheckNil(Grad,L_i):
    print
    i=1
    for (g,a) in zip(Grad,L_i):
        sq=a*g
        sq.expand()
        print 'L_'+str(i)+' * {\\nabla}f_'+str(i)+'=', sq
        i+=1

def energy(self):
    return self.subs(-m**2-p_x**2-p_y**2-p_z**2,-E**2).subs(m**2+p_x**2+p_y**2+p_z**2,E**2)
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
i=1
F={}
for signs in [(+1, +1), (-1, -1), (+1, -1), (-1, +1)]:
    kprop=[signs[0]*p_x,signs[0]*p_y,signs[0]*p_z,signs[1]*imag*m]
    func=S(0)
    for (dim,var) in zip(kprop,r):
        func = func + S(1)*var*dim
    func=-imag*(E*t-func)
    F[i] = S(1)*func
    i += 1
#Solutions 'simples' x4
i=1
K={}
for (signs,symbols) in zip([(-1, +1), (+1, -1), (+1, +1), (-1, -1)],
                           [('-','+'),('+','-'),('+','+'),('-','-')]):
    k=jcquat*imag*(kcquat*E+signs[0]*jcquat*m+signs[1]*icquat*pmulti)
    k.texLabel=jcquat.texLabel+imag.texLabel+'('+kcquat.texLabel+'E'+symbols[0]+jcquat.texLabel+'m'+symbols[1]+icquat.texLabel+pmulti.texLabel+')'
    K[i] = k
    i += 1
#Solutions 'combinees' x16
i=1
L ={}
R ={}
L1={}
R1={}
L2={}
R2={}
L3={}
R3={}
L4={}
R4={}
for (signs, symbols) in zip([(+1,  -1, -1), (+1, +1, +1), (-1, +1, -1), (-1, -1, +1)],
                            [('+','-','-'),('+','+','+'),('-','+','-'),('-','-','+')]):
    l=( signs[0]*K[1]+signs[1]*K[4]*icquat+signs[2]*K[3]*jcquat+K[2]*kcquat)
    r=( signs[0]*K[1]+signs[1]*icquat*K[4]+signs[2]*jcquat*K[3]+kcquat*K[2])
    l.texLabel=symbols[0]+'K_1'+symbols[1]+'K_4'+icquat.texLabel+symbols[2]+'K_3'+jcquat.texLabel+'+K_2'+kcquat.texLabel
    r.texLabel=symbols[0]+'K_1'+symbols[1]+icquat.texLabel+'K_4'+symbols[2]+jcquat.texLabel+'K_3'+'+'+kcquat.texLabel+'K_2'
    L1[i]=l
    R1[i]=r
    l=( signs[0]*K[2]+signs[1]*K[3]*icquat+signs[2]*K[4]*jcquat+K[1]*kcquat)
    r=( signs[0]*K[2]+signs[1]*icquat*K[3]+signs[2]*jcquat*K[4]+kcquat*K[1])
    l.texLabel=symbols[0]+'K_2'+symbols[1]+'K_3'+icquat.texLabel+symbols[2]+'K_4'+jcquat.texLabel+'+K_1'+kcquat.texLabel
    r.texLabel=symbols[0]+'K_2'+symbols[1]+icquat.texLabel+'K_3'+symbols[2]+jcquat.texLabel+'K_4'+'+'+kcquat.texLabel+'K_1'
    L2[i]=l
    R2[i]=r
    l=( signs[0]*K[3]+signs[1]*K[2]*icquat+signs[2]*K[1]*jcquat+K[4]*kcquat)
    r=( signs[0]*K[3]+signs[1]*icquat*K[2]+signs[2]*jcquat*K[1]+kcquat*K[4])
    l.texLabel=symbols[0]+'K_3'+symbols[1]+'K_2'+icquat.texLabel+symbols[2]+'K_1'+jcquat.texLabel+'+K_4'+kcquat.texLabel
    r.texLabel=symbols[0]+'K_3'+symbols[1]+icquat.texLabel+'K_2'+symbols[2]+jcquat.texLabel+'K_1'+'+'+kcquat.texLabel+'K_4'
    L3[i]=l
    R3[i]=r
    l=( signs[0]*K[4]+signs[1]*K[1]*icquat+signs[2]*K[2]*jcquat+K[3]*kcquat)
    r=( signs[0]*K[4]+signs[1]*icquat*K[1]+signs[2]*jcquat*K[2]+kcquat*K[3])
    l.texLabel=symbols[0]+'K_4'+symbols[1]+'K_1'+icquat.texLabel+symbols[2]+'K_2'+jcquat.texLabel+'+K_3'+kcquat.texLabel
    r.texLabel=symbols[0]+'K_4'+symbols[1]+icquat.texLabel+'K_1'+symbols[2]+jcquat.texLabel+'K_2'+'+'+kcquat.texLabel+'K_3'
    L4[i]=l
    R4[i]=r
    i+=1
L[1]=L1
R[1]=R1
L[2]=L2
R[2]=R2
L[3]=L3
R[3]=R3
L[4]=L4
R[4]=R4

#Equation d'onde (Onde Plane Progressive Monochromatique)
#psi=A1a*exp(f1)

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

print('Left and right derivatives')
print('http://geocalc.clas.asu.edu/pdf/NFMPchapt1.pdf see chapter 1-3 Differentiation by vectors')
print('{\\nabla}{\psi}={\\nabla}(Ke^f)')
print('{\\nabla}{\psi}=({\\nabla}K)e^f+K{\\nabla}(e^f)')
print('{\\nabla}K=0')
print('{\\nabla}{\psi_L}=K({\\nabla}f)e^f')
print('{\\nabla}{\psi_R}=({\\nabla}f)Ke^f')

print('The following symbols are defined :')
print('$E \\in \\mathbb{R}$')
print('$m \\in \\mathbb{R}$')
print('$'+pmulti.texLabel+'$ is defined with $p_x, p_y, p_z \\in \\mathbb{R}$')
print pmulti.texLabel+'=p_x'+iquat.texLabel+'+p_y'+jquat.texLabel+'+p_z'+kquat.texLabel
print pmulti.texLabel+'=', pmulti

print('Exponential function $f$')
i=1
for f in F.values():
    print 'f_'+str(i)+'=', f
    i += 1

print('Gradient for $f$')
i=1
for f in F.values():
    print '{\\nabla}f_'+str(i)+' = ', f.grad()
    i += 1

print('Square of the gradient')
Grad=[F[1].grad(),
      F[2].grad(),
      F[3].grad(),
      F[4].grad()]
CheckComputedGrad(Grad)

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

print('Simple Constants $K$ (exactly the gradient)')
for (kkey,k) in K.items():
    klabel = 'K_'+str(kkey)
    print klabel+'='+k.texLabel
#    k.convert_from_blades()
#    print klabel+'=', k

if outputTex:
    tex.Format('1 1 1 2')
else:
    MV.set_str_format(1)
print '\\end{equation*}\\newpage\\begin{equation*}'
print('Mixed Constants $A$ (built from simple constants and the coquaternions)')
for (column,l) in L.items():
    print '$L_'+str(column)+'$'
    theLeftVector  = '\\begin{pmatrix}'
    theRightVector = '\\begin{pmatrix}'
    for (row,element) in l.items():
        theLeftVector +='L_'+str(column)+'^'+str(row)+'\\\\'
        theRightVector+=element.texLabel+'\\\\'
    theLeftVector += '\\end{pmatrix}'
    theRightVector+= '\\end{pmatrix}'
    print theLeftVector+'='+theRightVector
for (column,r) in R.items():
    print '$R_'+str(column)+'$'
    theLeftVector  = '\\begin{pmatrix}'
    theRightVector = '\\begin{pmatrix}'
    for (row,element) in r.items():
        theLeftVector +='R_'+str(column)+'^'+str(row)+'\\\\'
        theRightVector+=element.texLabel+'\\\\'
    theLeftVector += '\\end{pmatrix}'
    theRightVector+= '\\end{pmatrix}'
    print theLeftVector+'='+theRightVector

def symmetryLeftRight(p_rotor, p_fromColumn, p_toColumn, p_toRows, p_sign):
    theLeftVector  = p_rotor.texLabel + '\\begin{pmatrix}'
    theTargetVector= '\\begin{pmatrix}'
    theRightVector = '\\begin{pmatrix}'
    for (lkey, lvalue) in L[p_fromColumn].items():
        theLeftVector += 'L_'+str(p_fromColumn)+'^'+str(lkey)+'\\\\'
        theTargetVector += 'L_'+str(p_toColumn)+'^'+str(p_toRows[lkey-1])+'\\\\'
        product = p_rotor*lvalue*p_rotor + L[p_toColumn][p_toRows[lkey-1]]*p_sign
        if product == 0:
            theRightVector += '0'
        theRightVector += '\\\\'
    theLeftVector += '\\end{pmatrix}' + p_rotor.texLabel
    theTargetVector+='\\end{pmatrix}'
    theRightVector+= '\\end{pmatrix}'
    if p_sign > 0:
        print theLeftVector+'+'+theTargetVector+'='+theRightVector
    else:
        print theLeftVector+'-'+theTargetVector+'='+theRightVector

def symmetryLeft(p_rotor, p_column, p_toRows, p_signs):
    theLeftVector  = p_rotor.texLabel + '\\begin{pmatrix}'
    theTargetVector= '\\begin{pmatrix}'
    theRightVector = '\\begin{pmatrix}'
    for (lkey, lvalue) in L[p_column].items():
        theLeftVector += 'L_'+str(p_column)+'^'+str(lkey)+'\\\\'
        if p_signs[lkey-1] > 0:
            theTargetVector += 'L_'+str(p_column)+'^'+str(p_toRows[lkey-1])+'\\\\'
        else:
            theTargetVector += '-L_'+str(p_column)+'^'+str(p_toRows[lkey-1])+'\\\\'
        product = p_rotor*lvalue - L[p_column][p_toRows[lkey-1]]*p_signs[lkey-1]
        if product == 0:
            theRightVector += '0'
        theRightVector += '\\\\'
    theLeftVector += '\\end{pmatrix}'
    theTargetVector+='\\end{pmatrix}'
    theRightVector+= '\\end{pmatrix}'
    print theLeftVector+'-'+theTargetVector+'='+theRightVector

def symmetryRight(p_rotor, p_fromColumn, p_toColumn, p_toRows, p_signs):
    theLeftVector  = '\\begin{pmatrix}'
    theTargetVector= '\\begin{pmatrix}'
    theRightVector = '\\begin{pmatrix}'
    for (lkey, lvalue) in L[p_fromColumn].items():
        theLeftVector += 'L_'+str(p_fromColumn)+'^'+str(lkey)+'\\\\'
        if p_signs[lkey-1] > 0:
            theTargetVector += 'L_'+str(p_toColumn)+'^'+str(p_toRows[lkey-1])+'\\\\'
        else:
            theTargetVector += '-L_'+str(p_toColumn)+'^'+str(p_toRows[lkey-1])+'\\\\'
        product = lvalue*p_rotor - L[p_toColumn][p_toRows[lkey-1]]*p_signs[lkey-1]
        if product == 0:
            theRightVector += '0'
        theRightVector += '\\\\'
    theLeftVector += '\\end{pmatrix}'+p_rotor.texLabel
    theTargetVector+='\\end{pmatrix}'
    theRightVector+= '\\end{pmatrix}'
    print theLeftVector+'-'+theTargetVector+'='+theRightVector

def symmetries():
    for (source, destination) in [(1, [4, 3, 2]),
                                  (2, [3, 4, 1]),
                                  (3, [2, 1, 4]),
                                  (4, [1, 2, 3])]:
        print '\\end{equation*}\\newpage\\begin{equation*}'
        print 'Symmetry $L_'+str(source)+'$'
        symmetryLeftRight(icquat, source, destination[0], [3, 4, 1, 2], -1)
        symmetryLeftRight(jcquat, source, destination[1], [4, 3, 2, 1],  1)
        symmetryLeftRight(kcquat, source, destination[2], [2, 1, 4, 3], -1)
        symmetryLeftRight(imag,   source, destination[1], [4, 3, 2, 1],  1)
        
        symmetryLeft(icquat, source, [1, 2, 3, 4], [1, -1, 1, -1])
        symmetryLeft(jcquat, source, [2, 1, 4, 3], [-1, 1, 1, -1])
        symmetryLeft(kcquat, source, [2, 1, 4, 3], [1, 1, -1, -1])
        
        symmetryRight(icquat, source, destination[0], [3, 4, 1, 2], [1, -1, 1, -1])
        symmetryRight(jcquat, source, destination[1], [3, 4, 1, 2], [-1, 1, 1, -1])
        symmetryRight(kcquat, source, destination[2], [1, 2, 3, 4], [1, 1, -1, -1])

print('$L_i^jK_i=0$')
print('$K_iR_i^j=0$')
print('$L_i^jR_i^j=0$')
def NullDetails():
    for index in range(1,5):
        for (id,element) in L[index].items():
            res=element*K[index]
            res.expand()
            res=res.energy()
            if res==0:
                print 'L_'+str(index)+'^'+str(id)+'K_'+str(index)+'=0'
            else:
                print 'L_'+str(index)+'^'+str(id)+'K_'+str(index)+'=',res
        for (id,element) in R[index].items():
            res=K[index]*element
            res.expand()
            res=res.energy()
            if res==0:
                print 'K_'+str(index)+'R_'+str(index)+'^'+str(id)+'=0'
            else:
                print 'K_'+str(index)+'R_'+str(index)+'^'+str(id)+'=',res
            for lrow in range(1,5):
                res=L[index][lrow]*element*S(1)/2
                res.expand()
                res=res.energy()
                if res==0:
                    print 'L_'+str(index)+'^'+str(lrow)+'R_'+str(index)+'^'+str(id)+'=0'
                else:
                    print 'L_'+str(index)+'^'+str(lrow)+'R_'+str(index)+'^'+str(id)+'=', res
#NullDetails()
#symmetries()
if outputTex:
    tex.Format('1 1 1 2')
else:
    MV.set_str_format(1)

print '\\end{equation*}\\newpage\\begin{equation*}'
print('Details about full solutions (built from mixed constants and the imaginary unit)')
psiL={}
for kindex in range(1,5):
#for kindex in range(1,2):
    index=1
    print '$\\psi_{L'+str(kindex)+'}$'
    temppsi={}
    for (akey, avalue) in L[kindex].items():
        for (lkey, lvalue) in L[kindex].items():
            psiLabel='\\psi_{L'+str(kindex)+'}^{'+str(index)+'}'
            print psiLabel+'\\equiv(L_'+str(kindex)+'^'+str(akey)+'+'+imag.texLabel+'L_'+str(kindex)+'^'+str(lkey)+')='+avalue.texLabel+'+'+imag.texLabel+'('+lvalue.texLabel+')'
            aPsi=(avalue+imag*lvalue)
            aPsi.convert_from_blades()
            aPsi.texLabel='(L_'+str(kindex)+'^'+str(akey)+'+'+imag.texLabel+'L_'+str(kindex)+'^'+str(lkey)+')'
            #print psiLabel+'=', aPsi
            temppsi[index]=aPsi
            product = aPsi*K[kindex]
            product.expand()
            product = product.energy()
            if product == 0:
                print psiLabel+'K_'+str(kindex)+'=', 0
            else:
                print psiLabel+'K_'+str(kindex)+'=', product
            index += 1
            psiLabel='\\psi_{L'+str(kindex)+'}^{'+str(index)+'}'
            print psiLabel+'\\equiv(L_'+str(kindex)+'^'+str(akey)+'-'+imag.texLabel+'L_'+str(kindex)+'^'+str(lkey)+')='+avalue.texLabel+'-'+imag.texLabel+'('+lvalue.texLabel+')'
            aPsi=(avalue-imag*lvalue)
            aPsi.convert_from_blades()
            aPsi.texLabel='(L_'+str(kindex)+'^'+str(akey)+'-'+imag.texLabel+'L_'+str(kindex)+'^'+str(lkey)+')'
            #print psiLabel+'=', aPsi
            temppsi[index]=aPsi
            product = aPsi*K[kindex]
            product.expand()
            product = product.energy()
            if product == 0:
                print psiLabel+'K_'+str(kindex)+'=', 0
            else:
                print psiLabel+'K_'+str(kindex)+'=', product
            index += 1
    psiL[kindex]=temppsi
print '\\end{equation*}\\newpage\\begin{equation*}'
psiR={}
for kindex in range(1,5):
#for kindex in range(1,2):
    index=1
    print '$\\psi_{R'+str(kindex)+'}$'
    temppsi={}
    for (akey, avalue) in R[kindex].items():
        for (lkey, lvalue) in R[kindex].items():
            psiLabel='\\psi_{R'+str(kindex)+'}^{'+str(index)+'}'
            print psiLabel+'\\equiv(R_'+str(kindex)+'^'+str(akey)+'+R_'+str(kindex)+'^'+str(lkey)+imag.texLabel+')='+avalue.texLabel+'+('+lvalue.texLabel+')'+imag.texLabel
            aPsi=(avalue+lvalue*imag)
            aPsi.convert_from_blades()
            aPsi.texLabel='(R_'+str(kindex)+'^'+str(akey)+'+R_'+str(kindex)+'^'+str(lkey)+imag.texLabel+')'
            #print psiLabel+'=', aPsi
            temppsi[index]=aPsi
            product = K[kindex]*aPsi
            product.expand()
            product = product.energy()
            if product == 0:
                print 'K_'+str(kindex)+psiLabel+'=0'
            else:
                print 'K_'+str(kindex)+psiLabel+'=', product
            index += 1
            psiLabel='\\psi_{R'+str(kindex)+'}^{'+str(index)+'}'
            print psiLabel+'\\equiv(R_'+str(kindex)+'^'+str(akey)+'-R_'+str(kindex)+'^'+str(lkey)+imag.texLabel+')='+avalue.texLabel+'-'+'('+lvalue.texLabel+')'+imag.texLabel
            aPsi=(avalue-lvalue*imag)
            aPsi.convert_from_blades()
            aPsi.texLabel='(R_'+str(kindex)+'^'+str(akey)+'-'+'R_'+str(kindex)+'^'+str(lkey)+imag.texLabel+')'
            #print psiLabel+'=', aPsi
            temppsi[index]=aPsi
            product = K[kindex]*aPsi
            product.expand()
            product = product.energy()
            if product == 0:
                print 'K_'+str(kindex)+psiLabel+'=0'
            else:
                print 'K_'+str(kindex)+psiLabel+'=', product
            index += 1
    psiR[kindex]=temppsi

if outputTex:
    tex.Format('1 1 1 1')
else:
    MV.set_str_format(0)
#for kindex in range(1,5):
for kindex in range(1,2):
    print('Symmetry with $\\psi_{L'+str(kindex)+'}$ solutions')
    print [1, 2, 3, 4, 9, 10, 11, 12]
    for groupe in [[ 1,  2],
                   [11, 12],
                   [ 3, 10],
                   [ 9,  4]]:
        text = imag.texLabel+'\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}+\\psi_{L'+str(kindex)+'}^{'+str(groupe[1])+'}='
        print text, imag*psiL[kindex][groupe[0]]+psiL[kindex][groupe[1]]
    for groupe in [[1, 2, 11, 12],
                   [3, 10],
                   [4, 9]]:
        print groupe
        if len(groupe) == 2:
            text = icquat.texLabel+'\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}-\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}='
            print text, icquat*psiL[kindex][groupe[0]]-psiL[kindex][groupe[0]]
            text = jcquat.texLabel+'\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}+\\psi_{L'+str(kindex)+'}^{'+str(groupe[1])+'}='
            print text, jcquat*psiL[kindex][groupe[0]]+psiL[kindex][groupe[1]]
            text = kcquat.texLabel+'\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}-\\psi_{L'+str(kindex)+'}^{'+str(groupe[1])+'}='
            print text, kcquat*psiL[kindex][groupe[0]]-psiL[kindex][groupe[1]]
        else:
            text = icquat.texLabel+'\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}-\\psi_{L'+str(kindex)+'}^{'+str(groupe[1])+'}='
            print text, icquat*psiL[kindex][groupe[0]]-psiL[kindex][groupe[1]]
            text = jcquat.texLabel+'\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}+\\psi_{L'+str(kindex)+'}^{'+str(groupe[2])+'}='
            print text, jcquat*psiL[kindex][groupe[0]]+psiL[kindex][groupe[2]]
            text = kcquat.texLabel+'\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}-\\psi_{L'+str(kindex)+'}^{'+str(groupe[3])+'}='
            print text, kcquat*psiL[kindex][groupe[0]]-psiL[kindex][groupe[3]]
    print [5, 6, 15, 16, 17, 18, 27, 28]
    for groupe in [[ 5, 18],
                   [17,  6],
                   [15, 28],
                   [27, 16]]:
        text = imag.texLabel+'\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}+\\psi_{L'+str(kindex)+'}^{'+str(groupe[1])+'}='
        print text, imag*psiL[kindex][groupe[0]]+psiL[kindex][groupe[1]]
    print [7, 8, 13, 14, 19, 20, 25, 26]
    for groupe in [[ 7, 26],
                   [25,  8],
                   [13, 20],
                   [19, 14]]:
        text = imag.texLabel+'\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}+\\psi_{L'+str(kindex)+'}^{'+str(groupe[1])+'}='
        print text, imag*psiL[kindex][groupe[0]]+psiL[kindex][groupe[1]]
    for groupe in [[5, 6, 16, 15],
                   [7, 13],
                   [8, 14]]:
        print groupe
        if len(groupe) == 2:
            text = icquat.texLabel+'\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}-\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}='
            print text, icquat*psiL[kindex][groupe[0]]-psiL[kindex][groupe[0]]
            text = jcquat.texLabel+'\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}+\\psi_{'+str(kindex)+'}^{'+str(groupe[1])+'}='
            print text, jcquat*psiL[kindex][groupe[0]]+psiL[kindex][groupe[1]]
            text = kcquat.texLabel+'\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}-\\psi_{L'+str(kindex)+'}^{'+str(groupe[1])+'}='
            print text, kcquat*psiL[kindex][groupe[0]]-psiL[kindex][groupe[1]]
        else:
            text = icquat.texLabel+'\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}-\\psi_{L'+str(kindex)+'}^{'+str(groupe[1])+'}='
            print text, icquat*psiL[kindex][groupe[0]]-psiL[kindex][groupe[1]]
            text = jcquat.texLabel+'\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}+\\psi_{L'+str(kindex)+'}^{'+str(groupe[2])+'}='
            print text, jcquat*psiL[kindex][groupe[0]]+psiL[kindex][groupe[2]]
            text = kcquat.texLabel+'\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}-\\psi_{L'+str(kindex)+'}^{'+str(groupe[3])+'}='
            print text ,kcquat*psiL[kindex][groupe[0]]-psiL[kindex][groupe[3]]
    for groupe in [[17, 18, 28, 27],
                   [19, 25],
                   [20, 26]]:
        print groupe
        if len(groupe) == 2:
            text = icquat.texLabel+'\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}-\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}='
            print text, icquat*psiL[kindex][groupe[0]]-psiL[kindex][groupe[0]]
            text = jcquat.texLabel+'\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}-\\psi_{L'+str(kindex)+'}^{'+str(groupe[1])+'}='
            print text, jcquat*psiL[kindex][groupe[0]]-psiL[kindex][groupe[1]]
            text = kcquat.texLabel+'\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}+\\psi_{L'+str(kindex)+'}^{'+str(groupe[1])+'}='
            print text, kcquat*psiL[kindex][groupe[0]]+psiL[kindex][groupe[1]]
        else:
            text = icquat.texLabel+'\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}-\\psi_{L'+str(kindex)+'}^{'+str(groupe[1])+'}='
            print text, icquat*psiL[kindex][groupe[0]]-psiL[kindex][groupe[1]]
            text = jcquat.texLabel+'\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}-\\psi_{L'+str(kindex)+'}^{'+str(groupe[2])+'}='
            print text, jcquat*psiL[kindex][groupe[0]]-psiL[kindex][groupe[2]]
            text = kcquat.texLabel+'\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}+\\psi_{L'+str(kindex)+'}^{'+str(groupe[3])+'}='
            print text, kcquat*psiL[kindex][groupe[0]]+psiL[kindex][groupe[3]]
    print [21, 22, 23, 24, 29, 30, 31, 32]
    for groupe in [[21, 22],
                   [31, 32],
                   [23, 30],
                   [29, 24]]:
        text = imag.texLabel+'\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}+\\psi_{L'+str(kindex)+'}^{'+str(groupe[1])+'}='
        print text, imag*psiL[kindex][groupe[0]]+psiL[kindex][groupe[1]]
    for groupe in [[21, 22, 31, 32],
                   [23, 30],
                   [24, 29]]:
        print groupe
        if len(groupe) == 2:
            text = icquat.texLabel+'\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}-\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}='
            print text, icquat*psiL[kindex][groupe[0]]-psiL[kindex][groupe[0]]
            text = jcquat.texLabel+'\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}-\\psi_{L'+str(kindex)+'}^{'+str(groupe[1])+'}='
            print text, jcquat*psiL[kindex][groupe[0]]-psiL[kindex][groupe[1]]
            text = kcquat.texLabel+'\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}+\\psi_{L'+str(kindex)+'}^{'+str(groupe[1])+'}='
            print text, kcquat*psiL[kindex][groupe[0]]+psiL[kindex][groupe[1]]
        else:
            text = icquat.texLabel+'\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}-\\psi_{L'+str(kindex)+'}^{'+str(groupe[1])+'}='
            print text, icquat*psiL[kindex][groupe[0]]-psiL[kindex][groupe[1]]
            text = jcquat.texLabel+'\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}-\\psi_{L'+str(kindex)+'}^{'+str(groupe[2])+'}='
            print text, jcquat*psiL[kindex][groupe[0]]-psiL[kindex][groupe[2]]
            text = kcquat.texLabel+'\\psi_{L'+str(kindex)+'}^{'+str(groupe[0])+'}+\\psi_{L'+str(kindex)+'}^{'+str(groupe[3])+'}='
            print text, kcquat*psiL[kindex][groupe[0]]+psiL[kindex][groupe[3]]

for (fromColumn, toColumns) in [(1, [4, 3, 2])
#                                ,(2, [3, 4, 1])
#                                ,(3, [2, 1, 4])
#                                ,(4, [1, 2, 3])
                                ]:
    print 'Rotations $\\psi_{L'+str(fromColumn)+'}$'
    for (fromRow,toRow) in [( 1, [22, 31, 12]),
                            ( 2, [21, 32, 11]),
                            ( 3, [24, 29, 10]),
                            ( 4, [23, 30,  9]),
                            ( 5, [18, 27, 16]),
                            ( 6, [17, 28, 15]),
                            ( 7, [20, 25, 14]),
                            ( 8, [19, 26, 13]),
                            ( 9, [30, 23,  4]),
                            (10, [29, 24,  3]),
                            (11, [32, 21,  2]),
                            (12, [31, 22,  1]),
                            (13, [26, 19,  8]),
                            (14, [25, 20,  7]),
                            (15, [28, 17,  6]),
                            (16, [27, 18,  5]),
                            (17, [ 6, 15, 28]),
                            (18, [ 5, 16, 27]),
                            (19, [ 8, 13, 26]),
                            (20, [ 7, 14, 25]),
                            (21, [ 2, 11, 32]),
                            (22, [ 1, 12, 31]),
                            (23, [ 4,  9, 30]),
                            (24, [ 3, 10, 29]),
                            (25, [14,  7, 20]),
                            (26, [13,  8, 19]),
                            (27, [16,  5, 18]),
                            (28, [15,  6, 17]),
                            (29, [10,  3, 24]),
                            (30, [ 9,  4, 23]),
                            (31, [12,  1, 22]),
                            (32, [11,  2, 21])
                            ]:
        psiLabel='\\psi_{L'+str(fromColumn)+'}^{'+str(fromRow)+'}'
        text = icquat.texLabel+psiLabel+icquat.texLabel+'-\\psi_{L'+str(toColumns[0])+'}^{'+str(toRow[0])+'}='
        print text, icquat*psiL[fromColumn][fromRow]*icquat-psiL[toColumns[0]][toRow[0]]
        text = jcquat.texLabel+psiLabel+jcquat.texLabel+'+\\psi_{L'+str(toColumns[1])+'}^{'+str(toRow[1])+'}='
        print text, jcquat*psiL[fromColumn][fromRow]*jcquat+psiL[toColumns[1]][toRow[1]]
        text = kcquat.texLabel+psiLabel+kcquat.texLabel+'-\\psi_{L'+str(toColumns[2])+'}^{'+str(toRow[2])+'}='
        print text, kcquat*psiL[fromColumn][fromRow]*kcquat-psiL[toColumns[2]][toRow[2]]

if outputTex:
    tex.xdvi(filename='evq.tex', debug=True)
