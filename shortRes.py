#################################################################################
# V1.0 April 16, 2025
# Author:
# Mario González-Sánchez. Contact: mario.gonzalez.sanchez@uva.es
# GitHub repository: https://github.com/mgonzalezsanchez/shortRes

# This package provides functions to compute the short resolution of a weighted
# homogeneous ideal (whose variables are in Noether position) in SageMath. 
# For simplicial toric rings of dimension 3, we provide an algorithm to compute 
# the shifts in a graded free resolution and a pruning algorithm to get the
# shifts in the short resolution

# ****************************************************************************
#    Copyright (C) 2025 Mario González-Sánchez <mario.gonzalez.sanchez@uva.es>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
# ****************************************************************************

# Required imports
import numpy

#################################################################################
# AUXILIARY FUNCTIONS (only for internal use of the main functions)
#################################################################################

def sortLex(r,K,x,n,B1_u):
    s = PolynomialRing(K,[x[n-3],x[n-2]], order = 'lex');
    L = [s(b) for b in B1_u];
    L.sort();
    return L

#################################################################################
# MAIN FUNCTIONS   
#################################################################################

def toricIdeal(A,K):
    r"""
    Computes the toric ideal determined by the columns of the matrix A (using 
    elimination) over the polynomial ring with coefficients in K.
    

    OUTPUT:
    
    The polynomial ring r (equipped with weights) and the toric ideal I

    EXAMPLES::
    
        sage: A = numpy.array([[1,6],[2,5],[6,1],[7,0],[0,7]]);
        sage: r,I = toricIdeal(A,QQ);
        sage: I
    
        Ideal (x2^2 - x1*x3) of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

    """
    n = A.shape[0]; m = A.shape[1]; 
    w = [*A.sum(axis=1)]; g = gcd([*w]);
    if g > 1:
        w = [f//g for f in w];
    wT = w + m*[1];
    if wT == (n+m)*[1]:
        r1 = PolynomialRing(K, ['x%s' % j for j in range(1,n+m+1)], order = 'degrevlex');
    else:
        r1 = PolynomialRing(K, ['x%s' % j for j in range(1,n+m+1)], order = TermOrder('wdegrevlex',tuple(wT)));
    x = r1.gens();
    IA = ideal([x[i] - r1.monomial(*n*[0]+[*A[i,:]]) for i in range(n)]);
    I = IA.elimination_ideal([*x[n:n+m]]);
    # Now we map I to a ring r with n variables
    if wT == (n+m)*[1]:
        r = PolynomialRing(K, ['x%s' % j for j in range(1,n+1)], order = 'degrevlex');
    else:
        r = PolynomialRing(K, ['x%s' % j for j in range(1,n+1)], order = TermOrder('wdegrevlex',tuple(w)));
    x = r.gens();
    phi = r1.hom(list(x)+m*[1], r);
    I = ideal([phi(g) for g in I.gens()]);
    return r,I

def shortRes(r,I,printlevel=0):
    r"""
    Computes the short resolution of the ring r/I, where I is a weighted homogeneous
    ideal. The weights have to be included in the ring, and the variables must be in
    Noether position
    
    OUTPUT: 
    
    A list containing the resolution. 
    When printlevel=1, the function prints the Betti diagram and the regularity 
    of r/I when it is standard graded.
    0 if the ideal is not homogeneous or the variables are not in Noether position

    EXAMPLES:: 
 
        sage: A = numpy.array([[7,2,3],[1,8,3],[3,8,1],[12,0,0],[0,12,0],[0,0,12]]);
        sage: r,I = toricIdeal(A,QQ);
        sage: res = shortRes(r,I,1);
        
                   0     1     2
        ------------------------
            0:     1     -     -
            1:     3     -     -
            2:     6     1     -
            3:    10     3     -
            4:    15     6     -
            5:    21    10     -
            6:    26    15     -
            7:    29    20     -
            8:    32    26     1
            9:    29    26     2
           10:    20    19     2
           11:     9     9     1
           12:     2     2     -
           13:     1     1     -
        ------------------------
        total:   204   138     6
        reg(R/I) = 13

    """
    if I.is_homogeneous() == 0:
        print('I is not homogeneous')
        return 0
    # We change the order just in case it's not wdegrevlex
    w = tuple([x.degree() for x in r.gens()]);
    if w == tuple(len(r.gens())*[1]):
        r = r.change_ring(order = TermOrder('degrevlex'));
    else:
        r = r.change_ring(order = TermOrder('wdegrevlex',tuple(w)));
    K = r.base_ring();
    x = r.gens(); n = len(x); d = I.dimension();
    last_var = x[n-d:n];
    G = I.groebner_basis();
    gen_inI = [g.lm() for g in G]; inI = ideal(gen_inI); 
    # 0th step of the resolution: computation of B0
    inI_ext = ideal(gen_inI+list(last_var)); 
    kdim = inI_ext.vector_space_dimension();
    if kdim == Infinity:
        print('dim(R/I) = '+str(d))
        print('The variables are not in Noether position')
        return 0
    B0 = inI_ext.normal_basis(algorithm = 'singular'); b0 = len(B0);
    # 1st step of the resolution: computation of B1'
    aux_vect = [*x[0:n-d],*d*[1]];
    XinI = ideal([f(aux_vect) for f in gen_inI]); 
    H = [g for g in B0 if g in XinI];
    B1 = [];
    for u in H:
        Iu = inI.quotient(ideal(u));
        Iu = Iu.elimination_ideal(*[x[0:n-d]]);
        B1_u = Iu.groebner_basis();
        B1 = B1 + [u*g for g in B1_u];
    # 1st step of the resolution: compute ker(psi_0)
    if w == tuple(len(r.gens())*[1]):
        s = PolynomialRing(K, last_var, order = TermOrder('degrevlex'));
    else:    
        s = PolynomialRing(K, last_var, order = TermOrder('wdegrevlex',tuple(w[n-d:n])));
    hp = [u-u.reduce(G) for u in B1];
    h = [];
    for f in hp:
        hf = b0*[0];
        l = len(f.monomials());
        while l > 0:
            inf = f.lm();
            Xinf = inf(aux_vect);
            i = B0.index(Xinf);
            hf[i] = hf[i] + f.lt()/Xinf;
            f = f - f.lt();
            l = len(f.monomials());
        hf = [s(g) for g in hf]; 
        h.append(hf);
    # Construct the module generated by these vectors h_i in Singular
    singular.setring(singular(s));
    syz1 = [singular.vector(g) for g in h];
    singular.lib('gradedModules.lib');
    degs = [g.degree() for g in B0];
    Syz1 = singular.grobj(singular.module(*syz1), vector(degs));
    # Finally, compute the graded resolution with singular and print the Betti diagram (if printlevel=1)
    res_Syz1 = singular.grres(Syz1,0,1);
    if printlevel == 1:
        display(singular.print(res_Syz1.betti(),' "betti" '))
        if w == tuple(len(r.gens())*[1]):
            print('reg(R/I) = '+str(singular.betti(res_Syz1).nrows()-1))
    return res_Syz1

def schreyerResDim3(r,I):
    r"""
    Computes sets of monomials whose degrees are the shifts in a graded free resolution 
    of r/I as A-module, where r/I is a simplicial toric ring of dimension 3. 
    The variables must be in Noether position.

    OUTPUT: 
    
    Lists B0, B1, B2 containing the monomials whose degrees are the shifts in the 
    resolution.
    0 if the ideal is not homogeneous, not toric (prime and binomial), the dimension 
    of r/I is not 3, or the variables are not in Noether position.

    EXAMPLES:: 
    
        sage: A = numpy.array([[7,2,3],[1,8,3],[3,8,1],[12,0,0],[0,12,0],[0,0,12]]);
        sage: r,I = toricIdeal(A,QQ);
        sage: B0,B1,B2 = schreyerResDim3(r,I); 
        sage: [len(B0),len(B1),len(B2)] 
    
        [204, 174, 42]
    
    """
    K = r.base_ring();
    x = r.gens(); n = len(x); d = I.dimension();
    if I.is_homogeneous() == 0:
        print('I is not homogeneous')
        return 0
    if d != 3:
        print('dim(R/I) != 3')
        return 0
    if I.is_prime() == 0:
        print('I is not prime')
        return 0
    G = I.groebner_basis();
    lG = [len(g.monomials()) for g in G];
    if lG != len(lG)*[2]:
        print('I is not binomial')
    gen_inI = [g.lm() for g in G]; inI = ideal(gen_inI);
    # Computation of B0
    last_var = x[n-d:n]; 
    inI_ext = ideal(gen_inI+list(last_var));
    kdim = inI_ext.vector_space_dimension();
    if kdim == Infinity:
        print('The variables are not in Noether position')
        return 0
    B0 = inI_ext.normal_basis(algorithm = 'singular'); b0 = len(B0);
    # Computation of B1' and B2'
    aux_vect = [*x[0:n-d],*d*[1]];
    XinI = ideal([f(aux_vect) for f in gen_inI]); 
    H = [g for g in B0 if g in XinI];
    B1 = []; B2 = [];
    for u in H:
        Iu = inI.quotient(ideal(u));
        Iu = Iu.elimination_ideal(*[x[0:n-d]]);
        B1_u = Iu.groebner_basis();
        # We sort the elements in B1 w.r.t. lex order x_{n-2}>x_{n-1}
        B1_u = sortLex(r,K,x,n,B1_u);
        b1_u = len(B1_u);
        B1 = B1 + [u*g for g in B1_u];
        if b1_u > 1:
            B2_u = [];
            for i in range(b1_u-1):
                B2_u.append(lcm(B1_u[i],B1_u[i+1]));
            B2 = B2 + [u*g for g in B2_u];
    return B0,B1,B2

def pruningDim3(r,I,B0,B1,B2):
    r"""
    Computes the sets of monomials whose degrees are the shifts in the short resolution of 
    r/I, a simplicial toric ring of dimension 3. The function takes as input the output of 
    schreyerResDim3, which computes the shifts in a (non-necessarily minimal) graded free 
    resolution of r/I as A-module

    OUTPUT: 
    
    Lists B0, B1_min c B1, B2_min c B2 containing the monomials whose degrees are the 
    shifts in the short resolution of r/I.

    EXAMPLES:: 
    
        sage: A = numpy.array([[7,2,3],[1,8,3],[3,8,1],[12,0,0],[0,12,0],[0,0,12]]);
        sage: r,I = toricIdeal(A,QQ);
        sage: B0,B1,B2 = shortResDim3(r,I); 
        sage: B0m,B1m,B2m = pruningDim3(r,I,B0,B1,B2)
        sage: [len(B0m),len(B1m),len(B2m)]
    
        [204, 138, 6]

    """
    x = r.gens(); n = len(x);
    G = I.groebner_basis();
    H1 = ideal([*G,x[n-3]]);
    GH1 = H1.groebner_basis();
    B1_min = [];
    for g in B1:
        f = g.exponents()[0];
        if f[n-3] == 0 and f[n-2] >= 2:
            gtail = g.reduce(G);
            ftail = gtail.exponents()[0];
            if ftail[n-3] > 0: # Element in C1
                gn2 = r(g/x[n-2]);
                if gn2.reduce(GH1) != 0:
                    B1_min.append(g);
            else: # Element in C2
                gtailn1 = r(gtail/x[n-1]);
                if (gn2.reduce(GH1) != 0) or (gtailn1.reduce(GH1) != 0):
                    B1_min.append(g);
        else:
            B1_min.append(g);
    B2_min = [];
    H2 = ideal([*G,x[n-1]]);
    GH2 = H2.groebner_basis();
    for g in B2:
        f = g.exponents()[0];
        gn2 = r(g/x[n-2]);
        if gn2.reduce(GH2) == 0:
            B2_min.append(g);
    return B0,B1_min,B2_min
