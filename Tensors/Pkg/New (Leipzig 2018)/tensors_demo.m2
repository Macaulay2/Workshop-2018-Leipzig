restart

--~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
--    TENSORS PROJECT - LEIPZIG, JUNE 7TH, 2018
--~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
installPackage "Tensors"
check "Tensors"

-- We have a new Type 'TensorSpace'
V = tensorSpace(QQ, X, {3,3,3})
peek V

-- We have a new Type 'Tensor'
T1 = makeTensor(27:1_QQ, V)
peek T1

T2 = 3*V_(0,0,0) + 2*V_(1,1,1) + 1/2*V_(1,2,1)
T2_(0,0,0)
T2_(2,2,2)

-- We can make the usual operations between tensors
3*T1 
peek oo

3*T1-T2
peek oo

-- We can use different coefficients
W = tensorSpace(QQ[a,b,c], Y, {2,2,2})
T3 = makeTensor(8:1_QQ, W)

-- We can produce bigger tensor spaces by using tensor products...
V ** W
W = tensorSpace(QQ, symbol Y, {2,2,2})

V ** W
peek oo

-- ...and we can do the tensor product of two tensors 
T1 = V_(0,0,0)
T2 = W_(0,0,0)

T1 ** T2
(T1 ** T2)#tensorSpace == V ** W

T1 ^** 3
-- ...names of new variables can be improved... we'll do it! Stay tuned.

-- We can look at slices of tensors
V = tensorSpace(QQ, symbol X, {2,3,4})
T = makeTensor(1..24, V);
M1 = slice(T, {,,3}); matrix entries M1
M2 = slice(T, {0,,}); matrix entries M2

-- and can do contraction of tensors, e.g. matrix product.
M3 = contraction(M1, M2, {1}, {0})
matrix entries M3 == (matrix entries M1) * (matrix entries M2)

-------------------------------------------------------------------------------------------------
-- EXERCISE 
-- Consider the tensor space QQ[3x3x3]:
--     - find the equations of the space of rank 1 tensors, 
--    	      i.e., the Segre embedding of PP^2 x PP^2 x PP^2 inside PP^26 = PP(QQ[3x3x3])
--     - find the equations of the space of symmetric rank 1 tensors, 
--    	      i.e., the Veronese variety of plane cubics in the linear subspace PP(Sym^3(QQ[3]))
-------------------------------------------------------------------------------------------------

-- We use the natural actions of:
--    1. GL(3) x GL(3) x GL(3):     (g1,g2,g3).(v1*v2*v3) = g1.v1 * g2.v2 * g3.v3;
--    2. GL(3):	       	       	    g.(v1*v2*v3) = (g,g,g).(v1*v2*v3)
a = symbol a; b = symbol b; c = symbol c;
ringT = QQ[a_0..a_8,b_0..b_8,c_0..c_8,Z_(0,0,0)..Z_(2,2,2)]
G1 = sub(genericMatrix(QQ[a_0..a_8],3,3),ringT)
G2 = sub(genericMatrix(QQ[b_0..b_8],3,3),ringT)
G3 = sub(genericMatrix(QQ[c_0..c_8],3,3),ringT)

-- We construct the tensor space...
use ringT
P222 = tensorSpace(ringT,symbol X,{3,3,3})

-- ...and a rank 1 (and symmetric) tensor
T = P222_(0,0,0)

-- We construct the generic tensor
generic222 = makeTensor(Z_(0,0,0)..Z_(2,2,2), P222)

-- We construct the orbit of the tensor with respect to the two actions
orbitT = glAction({G1,G2,G3},T)	    	    
orbitTsym = glAction(G1,T)	      	      	    

-- We find the equations of the orbits with elimination
ring222 = QQ[Z_(0,0,0)..Z_(2,2,2)]

-- Veronese
Isym = ideal (generic222 - orbitTsym)#coeff
sub(eliminate(Isym,toList(a_0..a_8)), ring222)
netList first entries gens oo

-- Segre
I = ideal (generic222 - orbitT)#coeff
sub(eliminate(I,toList(a_0..a_8 | b_0..b_8 | c_0..c_8)), ring222)
netList first entries gens oo
-------------------------------------------------------------------------------------------------
-- EXERCISE 
-- Consider the tensor space QQ[2x2x2x2]:
-- 
-------------------------------------------------------------------------------------------------
a = symbol a;
ringT = QQ[a_0..a_3,Z_(0,0,0,0)..Z_(1,1,1,1)]
G = sub(genericMatrix(QQ[a_0..a_3],2,2),ringT)

use ringT
P1111 = tensorSpace(ringT, symbol X, {2,2,2,2}) 

T = P1111_(0,0,0,0) + P1111_(1,1,1,1)
generic1111 = makeTensor(Z_(0,0,0,0)..Z_(1,1,1,1), P1111)

orbitTsym = glAction(G,T)
I = ideal (generic1111 - orbitTsym)#coeff

ring1111 = QQ[Z_(0,0,0,0)..Z_(1,1,1,1)]
sub(eliminate(I,toList(a_0..a_3)), ring1111)
netList first entries gens oo
-------------------------------------------------------------------------------------------------
-- EXERCISE 
-- Assume two tensors T1 and T2 in the space QQ[5x5x5] and Bernd guarantees you that there 
--     exists a unique element G in GL(5) such that G.T1 = T2. How can we compute it?
-------------------------------------------------------------------------------------------------
S = QQ[a_0..a_24];
G = genericMatrix(S,5,5)

V = tensorSpace(S, symbol X, {5,5,5})

t1 = {1/6, 1/4, 3/10, 1/3, 5/14, 1/6, 4/15, 1/3, 8/21, 5/12, 3/20, 1/4, 9/28, 3/8, 5/12, 2/15, 8/35, 3/10,
 16/45, 2/5, 5/42, 5/24, 5/18, 1/3, 25/66, 1/12, 2/15,1/6, 4/21, 5/24, 1/10, 1/6, 3/14, 1/4, 5/18, 1/10,
 6/35, 9/40, 4/15, 3/10, 2/21, 1/6, 2/9, 4/15, 10/33, 5/56, 10/63,3/14, 20/77, 25/84, 1/20, 1/12, 3/28, 1/8,
 5/36, 1/15, 4/35, 3/20, 8/45, 1/5, 1/14, 1/8, 1/6, 1/5, 5/22, 1/14, 8/63, 6/35, 16/77, 5/21, 5/72, 1/8,
 15/88, 5/24, 25/104, 1/30, 2/35, 3/40, 4/45, 1/10, 1/21, 1/12, 1/9, 2/15, 5/33, 3/56, 2/21, 9/70, 12/77, 5/28,
 1/18, 1/10, 3/22, 1/6, 5/26, 1/18, 10/99, 5/36, 20/117, 25/126, 1/42, 1/24, 1/18, 1/15, 5/66, 1/28, 4/63, 3/35,
 8/77, 5/42, 1/24, 3/40, 9/88, 1/8, 15/104, 2/45, 8/99, 1/9, 16/117, 10/63, 1/22, 1/12, 3/26, 1/7, 1/6};
t2 = {4775190420, 1238008758, 2619367102, 6286464106, 9028270177, 1431499812, 374499612, 788337868, 1877876116, 2695798084, 2665694614, 690571367,
 1463360546, 3512136023, 5056951682, 5727553416, 1474569334, 3135116908, 7564662774, 10890125698, 7119531710, 1800960792, 3880842820, 9480638704,
 13763752612, 995172450, 250598232, 543153951, 1329073106, 1943631652, 312378696, 79939860, 171110532, 413801588, 600387706, 581941703,
 147261930, 317851985, 775054978, 1129870989, 1201834586, 300225240, 654440012, 1610911884, 2362092740, 1595078258, 394098726, 865565139,
 2147727408, 3158369173, 2044300324, 520285677, 1117405280, 2715784181, 3942986554, 625615594, 161266698, 343262018, 826161842, 1194312148,
 1169140334, 298164317, 639518880, 1551586317, 2251831608, 2466335366, 622895070, 1345051146, 3287948790, 4786067768, 3189857811, 794136782,
 1732976622, 4279537392, 6263873656, 5643081938, 1463796838, 3095415147, 7427082040, 10660208142, 1686714860, 441467512, 928923072, 2212258600,
 3174573426, 3142942055, 814502386, 1725339897, 4140277890, 5959027177, 6771656846, 1744211336, 3706536960, 8941552620, 12865145080, 8398739600,
 2124724714, 4577686815, 11183684080, 16231267793, 5507131413, 1409953364, 3010760688, 7292754338, 10529347896, 1633846050, 422629268, 896963382,
 2154620492, 3107517924, 3077136845, 787614360, 1683654570, 4076874190, 5898688573, 6659462118, 1692103420, 3632491058, 8849336140, 12809319614,
 8441100992, 2107978403, 4585168771, 11306219093, 16493977500}

T1 = makeTensor(t1,V)
T2 = makeTensor(t2,V)

gT1 = glAction(G,T1);
I = ideal (T2-gT1)#coeff;

--now we try to solve the system exactly, using companion matrices:
time bb=sub(basis(S/I),S)--it takes about 25 minutes
L=for i in 0..24 list a_i=>0_QQ
for j in 0..24 do (
    compx=sub(contract(transpose bb,(bb_(0,0))*a_j%I),L);
    for i from 1 to (numcols bb-1) do compx=compx|sub(contract(transpose bb,(bb_(0,i))*a_j%I),L);
    fac_j=det (sub(compx,S)-a_j*sub(compx^0,S));
    print (factor fac_j)--third power of the solution a_j
    )
