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