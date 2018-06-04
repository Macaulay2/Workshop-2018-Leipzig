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
momentIdealExponential = (d,mix) ->(
    R:=QQ[lam_1..lam_mix,alp_1..alp_mix,m_0..m_d];
    I:=ideal (for i from 1 to d list -m_i+sum for j from 1 to mix list alp_j*lam_j^i)+ideal(-1+sum for i from 1 to mix list alp_i);
    eliminate for j from 1 to mix list alp_j,eliminate(for i from 1 to mix list lam_i ,I))
)