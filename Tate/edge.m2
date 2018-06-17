restart
needsPackage "EdgeIdeals"
--viewHelp EdgeIdeals
n= 10
MM = id_(ZZ^n)
S = QQ[a_1..a_n, Degrees=>apply(n, i->apply(n,j->MM_(i,j)))]
describe S
G = cycle(S,n)
I = edgeIdeal G
betti( F=res I)
netList degrees F_2
H = path(4,S)
