restart
path = path|{"/Users/andrew/Dropbox/Macaulay2","/Users/andrew/SVN/WFU-2012/"};
debug needsPackage("Tensors")

{{*------------------
-- Making a tensor --
------------------*}};
R=QQ[a..z]
S=QQ[x,y]

{{* from a nested list *}};
T=makeTensor{
     {{a,b},{c,d}},
     {{e,f},{g,h}},
     {{i,j},{k,l}}}

{{* from a list and given dimensions *}};
makeTensor({3,2,2},{a,b,c,d,e,f,g,h,i,j,k,l})

{{* generic tensors *}};
T=genericTensor(R,{3,2,2})
U=genericTensor(R,12,{3,2,2})
genericTensor(R,m,{3,2,2})

{{* random tensors of numbers *}};
randomTensor(ZZ,{4,2,2})
randomTensor(QQ,{4,2,2})

{{* random tensors of homogeneous polynomials *}};
randomTensor(S,1,{4,2,2})

{{*------------------------
-- Operations on tensors --
------------------------*}};
3*T+U

{{*----------
-- Entries --
----------*}};
T
T_5 -- by ordinal
T_(1,0,1) -- by position

{{*-----------------------
-- Permutation action: ---
-----------------------*}};
T
T@{0,1,2}--same as T
T@{2,0,1}--different!

{{*The action is designed to be
equivariant with respct to the 
action on the keys: *}}
I=(0,1,0);perm={2,0,1};
T_I===(T@perm)_(toSequence I_perm)

{{*---------
-- Slices --
---------*}};

{{* to get slices of a tensor, use a 
list subscript with blank spots (nulls)
for unspecified indices: *}};
T
T_{,0,}
t=T
l={,0,}
T
T_{1,,}
T_{1,1,1}

{{*------------
-- Marginals --
------------*}};
T
marginalize(T,{0})
marginalize(T,{0,2})
marginalize(T,{0,1,2})
a+b+c+d+e+f+g

{{*-------------------------
-- Indexed Tensor Methods --
-------------------------*}};
R=QQ[a..z]
S=QQ[x,y]
i=symbol i;j=symbol j;k=symbol k;l=symbol l;

T=genericTensor(R,{3,2,2})
U=randomTensor(ZZ,{2,2})

tensor T_(0,i,j)
tensor sum(k,T_(i,j,k) * U_(k,l))

{{*how indexed tensors work *}}
t=T_(i,j,k)
u=U_(k,l)
peek t
peek u
t*u
sum(k,t*u)
tensor sum(k, t*u)
tensor sum(k, T_(i,j,k) * U_(k,l))
tensor (T_(i,j,k) * U_(k,l))
tensor sum(i, T_(i,j,k) * U_(k,l))
tensor einsum{T_(i,j,k), U_(j,k)}

{{*---------------------------
-- Index permutation action --
---------------------------*}};
x=symbol x;i=symbol i;j=symbol j;k=symbol k;
R'=QQ[x_0..x_60]
T=genericTensor(R',{3,4,5})
t=T_(i,j,k)
t'=t@(k,i,j)--index permutations are sequences, not lists
T'=tensor t'
s=(k,i,j)
T
T_(2,3,0)
T'_(0,2,3)

{{*The action is designed to be
equivariant with respect to the 
action on the index sequence: *}}
I=(1,3,4)
perm={2,0,1}
t_I===(t@(toSequence (i,j,k)_perm))_(toSequence I_perm)
(tensor t)@perm===tensor(t@(toSequence (i,j,k)_perm))

{{*
time T_(i1,i2) * T_(i3,i4) * T_(i5,i6) * T_(i7,i8)
time indexedTensorProduct{T_(i1,i2),T_(i3,i4),T_(i5,i6),T_(i7,i8)}
*}}


{{*--------------------
--Modules of tensors
--------------------*}};

{{* a tensor is an element
of a "tensor module", *}};
M=class T

{{* which is a module that remembers
it is a tensor product of smaller 
modules.*}};
M0=R^3 ** R^2 ** R^2 -- doesn't remember it's a tensor product
M=tensorModule(R,{3,2,2}) -- remembers

M'=tensorModule(R,{2,2,3})
new TensorModule from M -- doesn't know it's free
M==M0--they are equal as modules, 
M===M0--but not as hashtables,
M==M'--and tensor modules with different factors are different

{{* tensor products of tensor modules *}};
N=tensorModule(R,{2,2})
P=M**N
P_7
P_(0,0,1,1,1)


end

