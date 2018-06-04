-- Lists all moments of the univariate Gaussian
listOfMoments = (d,R) -> (
  S := R[t]/t^(d+1);
  use S;
  g := gens R;
  a := g_0*t + 1/2 * g_1^2 * t^2;
  b := exp(a);
  li := for i from 1 to d list i! * coefficient(t^i,b);
  use R;
  li
)

--Gaussian
momentIdeal = d->(
    R=QQ[mn,sd,m_0..m_d][t]/t^(d+1);
    series:=exp(mn*t+(1/2)*sd^2*t^2);
    I:=ideal for i from 1 to d list i!*coefficient(t^i,series)-m_i;
    eliminate({mn,sd},I)
    )

--Exponential mixture
--takes highest  degree d of moments and number of mixtures
--NEEDS TO EB HOMOGENISED, NEEDS TO FIX DOUBLE ELIMINATE
momentIdealExponential = (d,mix) ->(
    R:=QQ[lam_1..lam_mix,alp_1..alp_mix,m_0..m_d];
    I:=ideal (for i from 1 to d list -m_i+sum for j from 1 to mix list alp_j*lam_j^i)+ideal(-1+sum for i from 1 to mix list alp_i);
    eliminate for j from 1 to mix list alp_j,eliminate(for i from 1 to mix list lam_i ,I))
)

--Gaussian Mixtures
--takes as input the number of mixtures and the highest degree of moments appearing
--computes the homogeneous moment ideal by eliminating the means and standard deviations
momentIdealGaussian = (mix,d)->(
    R=QQ[mn_1..mn_mix,sd_1..sd_mix,m_0..m_d][t]/t^(d+1);
    series:=sum for i from 1 to mix list exp(mn_i*t+(1/2)*sd_i^2*t^2);
    I:=ideal for i from 1 to d list i!*coefficient(t^i,series)-m_i;
    homogenize(eliminate((for i from 1 to mix list mn_i)|(for i from 1 to mix list sd_i),I),m_0)
)