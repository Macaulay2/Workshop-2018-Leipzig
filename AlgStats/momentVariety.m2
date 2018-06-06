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
    R=QQ[mn_1..mn_mix,sd_1..sd_mix,a_1..a_mix,m_0..m_d][t]/t^(d+1);
    series:=sum for i from 1 to mix list a_i*exp(mn_i*t+(1/2)*sd_i^2*t^2);
    I:=ideal for i from 1 to d list i!*coefficient(t^i,series)-m_i+ideal(-1+sum for i from 1 to mix list a_i);
    homogenize(eliminate((for i from 1 to mix list a_i)|(for i from 1 to mix list mn_i)|(for i from 1 to mix list sd_i),I),m_0)
)

-- This computes the moment map up to order d for an n-dimensional Gaussian. 
-- The parameters for the mean are x_1,..,x_n and for the covariance matrix are s_(i,j)

momentMapGaussians =  (n,d) -> (  
  par:=toList(x_1..x_n);
  for i from 1 to n do (for j from i to n do (par=append(par,s_(i,j))) );
  par=toSequence(par);
  R := QQ[par];
  mu := matrix({toList(x_1..x_n)});
  Sigma := genericSymmetricMatrix(R,s_(1,1),n);
     
  S := R[t_1..t_n]/((ideal(t_1..t_n))^(d+1));
  use S;
  a := vars(S)*transpose(mu) + (1/2) * vars(S)*Sigma*transpose(vars(S));
  MGF := exp(a_(0,0));
  
  
  (M,C):=coefficients(MGF);
  use R;
  C = mutableMatrix(C);
  for i from 0 to numColumns(M)-1 do (C_(i,0) = (((degree(M_(0,i)))_0)!)*C_(i,0));
  C = matrix(C);
  C=lift(C,R);
  return matrix({(reverse((entries(transpose(C)))_0))});   
)

--Gaussian Mixtures Test
--written to eliminate a_mix
momentIdealGaussianTest = (mix,d)->(
    R=QQ[mn_1..mn_mix,sd_1..sd_mix,a_1..a_mix,m_0..m_d][t]/t^(d+1);
    if mix == 1 then(
    	series:=sum for i from 1 to mix list a_i*exp(mn_i*t+(1/2)*sd_i^2*t^2);
    	I:=ideal for i from 1 to d list i!*coefficient(t^i,series)-m_i+ideal(a_1 - 1);
    	return homogenize(eliminate({a_1}|(for i from 1 to mix list mn_i)|(for i from 1 to mix list sd_i),I),m_0);)
    else( 
    	amix := 1 - sum for i from 1 to mix-1 list a_i;
    	series:=sum for i from 1 to mix-1 list a_i*exp(mn_i*t+(1/2)*sd_i^2*t^2) + amix*exp(mn_mix*t+(1/2)*sd_mix^2*t^2);
    	I:=ideal for i from 1 to d list i!*coefficient(t^i,series)-m_i;
    	return homogenize(eliminate((for i from 1 to mix-1 list a_i)|(for i from 1 to mix list mn_i)|(for i from 1 to mix list sd_i),I),m_0))
)