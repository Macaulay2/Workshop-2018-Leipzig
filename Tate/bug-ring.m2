restart
E = ZZ/101[a,b, SkewCommutative =>true]
T1=(dual res(coker gens prune (ideal vars E)^2,LengthLimit=>11))[1]
T2=(dual res(coker gens prune (ideal vars E)^2,LengthLimit=>11))
T3=(res(coker gens prune (ideal vars E)^2,LengthLimit=>11))
I2 = prune (ideal vars E)^2
I3 = (ideal vars E)^2
E === ring I2 --
E === ring I3


T3a=(res(prune module ((ideal vars E)^2),LengthLimit=>11))
T4 = res prune coker vars E
T5 = res coker vars E
E === ring T1
E === ring T2
E === ring T3
E === ring T3a
E === ring T4
restart
E = ZZ/101[a,b]
T1=(dual res(coker gens prune (ideal vars E)^2,LengthLimit=>11))[1]
T2=(dual res(coker gens prune (ideal vars E)^2,LengthLimit=>11))
T3=(res(coker gens prune (ideal vars E)^2,LengthLimit=>11))
T4 = res prune coker vars E
T5 = res coker vars E
E === ring T1
E === ring T2
E === ring T3
E === ring T4
E === ring T5

end--
restart
load "bug-ring.m2"

code methods prune
