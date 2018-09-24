---
--- Presentation in Leipzig, Juni 8th, 2018
---

restart
debug needsPackage "TateOnProducts"

N = {2,2}
(S, E) = productOfProjectiveSpaces N

M0 = cokernel matrix {{x_(0,0)^2*x_(1,1), x_(0,0)^2*x_(0,1)^2, x_(0,0)^2*x_(1,0)^3, x_(0,0)^3*x_(1,0)^2, x_(0,1)^3*x_(1,0)^2, x_(0,1)^2*x_(1,0)*x_(1,1)^3}}
multigradedRegularity M0


M = cokernel matrix {{x_(1,0)^2*x_(1,1), x_(0,0)*x_(0,1)^3, x_(0,0)^2*x_(0,1)*x_(1,1)^2, x_(0,0)*x_(0,1)^2*x_(1,1)^2, x_(0,1)^3*x_(1,0)^2, x_(0,0)^3*x_(1,0)^3}}
multigradedRegularity M

---

R = regularity M
(low, high) = (-{sum N, sum N}, {R, R})

m = cohomologyMatrix(M, low, high)
findCorners(m, low, high)

ht = cohomologyHashTable(M, low, high);
findCorners ht
