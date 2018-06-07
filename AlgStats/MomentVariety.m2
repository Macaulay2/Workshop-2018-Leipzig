--
-- PREAMBLE -----------------------------------------------------
-- -*- coding: utf-8 -*-
newPackage(
    "MomentVariety",
    Version => "0.1", 
    Date => "6 June 2018",
    Authors => {{Name => "Some people", 
	    Email => "your.name@email.edu", 
	    HomePage => "http://your.website.html"}},
    Headline => "a package that computes the moments of some distributions",
    AuxiliaryFiles => false,
    Reload => true,
    DebuggingMode => true
    )

-- EXPORT LIST --------------------------------------------------
export {
    "listOfMoments",
    "momentIdeal",
    "momentIdealExponential",
    "momentIdealGaussian",
    "momentMapGaussians",
    "momentVarietyGaussians",
    "momentIdealPoisson",
    "momentIdealMultinomial",
    "momentIdealMultinomialMixture",
    "formalLog",
    "cumulantIdealGaussian",
    "cumulantIdealExponential",
    "cumulantIdealPoisson",
    "cumulantsToMoments",
    "momentsToCumulants",
    "momentIdealToCumulants",
    "cumulantToMoments",
    "momentIdealToCumulantsMultivariate",
    "cumulantIdealToMomentsMultivariate"
    }

-- Lists all moments of the univariate Gaussian

GaussianMoments = method()
GaussianMoments (ZZ, Ring) := List => (d, R) -> (
  t := symbol t;
  S := R[t]/t^(d+1);
  use S;
  g := gens R;
  a := g_0*t + 1/2 * g_1^2 * t^2;
  b := exp(a);
  li := for i from 1 to d list i! * coefficient(t^i,b);
  use R;
  li
)

--Exponential mixture
--takes as input the number of mixtures, the highest degree d of moments and a ring
momentIdealExponential = method(Options => {K => QQ})
momentIdealExponential (ZZ, ZZ) := o -> (mix, d)->(
    l := local l;
    a := local a;
    m := local m;
    R := o.K[l_1..l_mix,a_1..a_mix,m_0..m_d];
    I := ideal (for i from 1 to d list -m_i+sum for j from 1 to mix list a_j*l_j^i*i!) + ideal(-1+sum for i from 1 to mix list a_i);
    I = homogenize(eliminate (toList(a_1..a_mix)|toList(l_1..l_mix) ,I),m_0);
    sub(I, o.K[m_0..m_d])
)

--------------------------------------------------------------------------------------

momentMapGaussians =  (n,d) -> (
  x := symbol x;    
  s := symbol s;
  par:=toList(x_1..x_n);
  for i from 1 to n do (for j from i to n do (par=append(par,s_(i,j))) );
  par=toSequence(par);
  R := QQ[par];
  mu := matrix({toList(x_1..x_n)});
  Sigma := genericSymmetricMatrix(R,s_(1,1),n);
  t := symbol t;   
  S := R[t_1..t_n]/((ideal(t_1..t_n))^(d+1));
  use S;
  a := vars(S)*transpose(mu) + (1/2) * vars(S)*Sigma*transpose(vars(S));
  MGF := exp(a_(0,0));
  
  
  (M,C):=coefficients(MGF);
  use R;
  C = mutableMatrix(C);
  lM :=  flatten (entries M);
  lexpM := flatten (apply(lM,mon->exponents(mon)));
  c := 1;
  for i from 0 to numColumns(M)-1 do (
      (for e in lexpM_i do c = c*(e!));
      -- (for m in ( (entries vars S)_0 ) do c = c*((degree(m,M_(0,i)))!));
      C_(i,0) = c*C_(i,0);
      c=1;
      );
  C = matrix(C);
  C=lift(C,R);
  m := symbol m;
  momvars := toSequence reverse (apply(lexpM,e->m_e));
  
  return (matrix({(reverse((entries(transpose(C)))_0))}),momvars);
     
)   	    	    	

-- This computes the homogeneous ideal of the moment variety.
momentVarietyGaussians = method()
momentVarietyGaussians (ZZ, ZZ) := Ideal => (n,d) -> (
    
  (C,momvars) := momentMapGaussians(n,d);   
  R := ring(C);
  k := coefficientRing(R);
    
  PPM := k[momvars];
  varmoms := gens PPM;
  f := map(R,PPM,C);
  I := kernel f;
  I = homogenize(I,varmoms_0);
  
  zeroes = for i from 0 to n - 1 list 0;
  full = for i from 0 to n - 1 list d;
  return sub(I,QQ[m_zeroes..m_full]);  
--  return I
   
)

-------------------------------------------------------------------------------------
--------------------------------------------------------------
-- This tries to compute mixtures of multivariate Gaussians

momentMapGaussiansMixtures = (n,d,k,KK) -> (
 
  x := symbol x;
  parx := {};
  parxtemp :={};
 
  s := symbol s;
  pars := {};
  parstemp := {};
 
  a := symbol a;
  para := {};
 
  -- This creates the variables for the moments.
  for i from 1 to k do (
      for j from 1 to n do (parxtemp = append(parxtemp,x_(i,j)););
      parx = append(parx,parxtemp);
      parxtemp = {};
      );
 
  -- This creates the variables for the covariances.
  for i from 1 to k do (
      for j from 1 to n do (
	  for h from j to n do (
	      parstemp = append(parstemp,s_(i,(j,h)));
	      );
	  );
      pars = append(pars,parstemp);
      parstemp = {};
      );
 
  -- This creates the variables for the mixture parameters.
  para = toList(a_1..a_k);
 
  -- This makes the ring with all the parameters.
  par = join(flatten(parx),flatten(pars),flatten(para));
  R := KK[par];
 
  -- This makes the generating function
  t := symbol t;
  auxvars = toList(t_1..t_n);
  auxring := R[auxvars];
  auxideal := ideal(vars(auxring));
  S = auxring/(auxideal^(d+1));
  use S;
 
  mu := 0;
  Sigma:= 0;
  MGF := 0;
  CGF := 0;
 
 
  for i from 1 to k do (
    mu = genericMatrix(R,x_(i,1),n,1);
    mu = promote(mu,S);
    Sigma= genericSymmetricMatrix(R,s_(i,(1,1)),n);
    Sigma = promote(Sigma,S);
    logarithm =  vars(S)*mu + (1/2) * vars(S)*Sigma*transpose(vars(S));
    MGF = MGF + (a_i)*exp(logarithm_(0,0));
    );
 
 
  (M,C):=coefficients(MGF);
  use R;
  C = mutableMatrix(C);
  lM :=  flatten (entries M);
  lexpM := flatten (apply(lM,mon->exponents(mon)));
  c := 1;
  for i from 0 to numColumns(M)-1 do (
      (for e in lexpM_i do c = c*(e!));
      -- (for m in ( (entries vars S)_0 ) do c = c*((degree(m,M_(0,i)))!));
      C_(i,0) = c*C_(i,0);
      c=1;
      );
  C = matrix(C);
  C=lift(C,R);
 
  momvars := toSequence reverse (apply(lexpM,e->m_e));
 
  return (matrix({(reverse((entries(transpose(C)))_0))}),momvars);
 
)


----------------------------------------------------------------
-- This computes the homogeneous ideal of the moment variety of multidimensional Gaussian mixtures

momentVarietyGaussiansMixtures = (n,d,k,KK) -> (
 
  (C,momvars) := momentMapGaussiansMixtures(n,d,k,KK);
  R := ring(C);
  k := coefficientRing(R);
 
  PPM := KK[momvars];
  varmoms := gens PPM;
  f := map(R,PPM,C);
  I := kernel f;
  I = homogenize(I,varmoms_0);
 
  return I;
 
)


-------------------------------------------------------------------------------------

--Poisson Mixtures
--takes as input the number of mixtures and the highest degree of moments appearing
--computes the homogeneous moment ideal 

momentIdealPoisson = method(Options => {K => QQ})
momentIdealPoisson (ZZ, ZZ) := o -> (mix, d)-> (
    lambda := symbol lambda;
    a := symbol a;
    m := symbol m;
    t := symbol t;
    R := o.K[lambda_1..lambda_mix,a_1..a_mix,m_0..m_d][t]/t^(d+1);
    use R;
    series := sum for i from 1 to mix list a_i*exp(lambda_i*(exp(t)-1));
    I := ideal for i from 1 to d list i!*coefficient(t^i,series)-m_i+ideal(-1+sum for i from 1 to mix list a_i);
    I = homogenize(eliminate((for i from 1 to mix list a_i)|(for i from 1 to mix list lambda_i),I),m_0);
    sub(I, o.K[m_0..m_d])
)

--Gaussian Mixtures Test
--written to eliminate a_mix
momentIdealGaussian = method(Options => {K => QQ})
momentIdealGaussian (ZZ, ZZ) := o -> (mix, d)->(
    mn := symbol  mn;
    sd := symbol sd;
    m :=  symbol m;
    t :=  symbol t;
    a := symbol a;
    if mix == 1 then(
	R := o.K[mn_1..mn_mix,sd_1..sd_mix,m_0..m_d][t]/t^(d+1);
	use R;
    	series:= exp(mn_1*t+(1/2)*sd_1^2*t^2);
    	I := ideal for i from 1 to d list i!*coefficient(t^i,series)-m_i;
    	I =  homogenize(eliminate((for i from 1 to mix list mn_i)|(for i from 1 to mix list sd_i),I),m_0);
	return sub(I, o.K[m_0..m_d])
	)
    else( 
	R2 := o.K[mn_1..mn_mix,sd_1..sd_mix,a_1..a_(mix-1),m_0..m_d][t]/t^(d+1);
	use R2;
    	amix := 1 - sum for i from 1 to mix-1 list a_i;

    	series2 := sum for i from 1 to mix-1 list a_i*exp(mn_i*t+(1/2)*sd_i^2*t^2) + amix*exp(mn_mix*t+(1/2)*sd_mix^2*t^2);
    	I2 := ideal for i from 1 to d list i!*coefficient(t^i,series2)-m_i;
    	I = homogenize(eliminate((for i from 1 to mix-1 list a_i)|(for i from 1 to mix list mn_i)|(for i from 1 to mix list sd_i),I2),m_0);
	return sub(I, o.K[m_0..m_d])
	)
)


--computing the moment ideal for the multinomial distribution
--r = #possible outcomes
--n = #trials
--p_1,..,p_r are the probabilities of each outcome so that their sum is 1
--t_1..t_r are the variables of the moment generating function
--d is the truncation order

momentIdealMultinomial = method(Options => {K => QQ})
momentIdealMultinomial (ZZ, ZZ, ZZ) := o -> (k, n, d) -> (
    t := symbol t;
    S := o.K[t_1..t_k];
    exps := flatten apply(toList(0..d), i->flatten entries basis(i,S) / exponents / flatten);
    quotientExps := flatten entries basis(d+1,S) / exponents / flatten;
    Mons := ideal(apply(quotientExps, e->S_e));
    p := symbol p;
    m := symbol m;
    R := o.K[p_1..p_k,apply(exps,i->m_i)][t_1..t_k];
    Mons = sub(Mons,R);
    R = R / Mons;
    use R;
    series := (sum apply(toList(1..k), j-> p_j*exp(t_j)))^n; --moment gen fxn of the multinomial distribution
    I := ideal( apply(exps, e-> (sum e)!*coefficient(sub(S_e,R),series)-m_e) ) + ideal( 1 - sum apply(toList(1..k), i -> p_i));
    T := o.K[apply(exps,i->m_i)];
    homogenize(sub((eliminate(toList(p_1..p_k),I),T)),m_(exps#0))
)

--Mixtures of Multinomial Distributions

-- param es a list with the varaibles k and n
momentIdealMultiMixture = method(Options => {K => QQ, Mixture => 1})
momentIdealMultiMixture (ZZ, ZZ, ZZ) := o -> (r, n, d) -> (

    mix := o.Mixture;
    t := symbol t;
    S := o.K[t_1..t_r];

    exps := flatten apply(toList(0..d), i->flatten entries basis(i,S) / exponents / flatten);
    quotientExps := flatten entries basis(d+1,S) / exponents / flatten;
    Mons := ideal(apply(quotientExps, e->S_e));
    --need different parameters p-1...p_k for each distribution in the mixture
    p := symbol p;
    a := symbol a;
    m := symbol m;
    R := o.K[p_(1,1)..p_(mix,r),a_1..a_mix,apply(exps,i->m_i)][t_1..t_r];
    Mons = sub(Mons,R);
    R = R / Mons;
    use R;
    series := sum apply(toList(1..mix),i->a_i*(sum apply(toList(1..r), j-> p_(i,j)*exp(t_j)))^n); --moment gen fxn of the multinomial distribution
    I := ideal( apply(exps, e-> (sum e)!*coefficient(sub(S_e,R),series)-m_e) ) + 
    	ideal( apply(toList(1..mix),i-> 1 - sum apply(toList(1..r), j -> p_(i,j)))) + 
	ideal(1 - sum apply(toList(1..mix),i -> a_i));
    T := o.K[apply(exps,i->m_i)];
    homogenize(sub((eliminate(toList(p_(1,1)..p_(mix,r))|toList(a_1..a_mix),I),T)),m_(exps#0))
)

--Binomial Distribution
--n trials, truncation order d
momentIdealBinomial = method()
momentIdealBinomial = (n,mix,d) -> momentIdealMultiMixture(2,n,mix,d)


--Moment Ideal from Moment Generating function
--takes as input the number of mixtures, the highest degree of moments appearing, a list with the MGF and the parameters of this function, and a Ring.
--computes the homogeneous moment ideal 

momentIdealFromMGF = method(Options => {K => QQ})
momentIdealFromMGF (ZZ, ZZ, Thing, List) := o -> (mix, d, f, param) -> (
    param := symbol param;
    n := #param - 1;
    m := symbol m;
    a := symbol a;
    paramMix := for i to n list (param_i)_1..(param_i)_mix;
    R := o.K[toSequence paramMix, toSequence param, a_1..a_mix, m_0..m_d];
    S := R[t]/t^(d+1);
    use S;
    paramSubs := flatten for i from 1 to mix list
    	for j to n list R_(param_j) => R_(paramMix_j_(i-1));
    f = sub(f,S);
    series := sum for i from 1 to mix list a_i*sub(f,paramSubs_(i-1));
    I := ideal for i from 1 to d list i!*coefficient(t^i,series)-m_i+ideal(-1+sum for i from 1 to mix list a_i);
    I = homogenize(eliminate((gens R)_{0..(#(gens R)-d-2)},I),m_0);
    sub(I, K[m_0..m_d])
)


--Cumulant ideals
formalLog = (f, d) -> (
    sum for k from 1 to d list (-1)^(k-1)/k * (f-1)^k
)

cumulantIdealGaussian = method()
cumulantIdealGaussian (ZZ,ZZ) := Ideal => (mix,d) -> (
    mn := symbol mn;
    sd := symbol sd;
    t := symbol t;
    k := symbol k;
    a := symbol a;
    R := QQ[mn_1..mn_mix,sd_1..sd_mix,a_1..a_mix,k_0..k_d][t]/t^(d+1);
    use R;
    series:=formalLog(sum for i from 1 to mix list a_i*exp(mn_i*t+(1/2)*sd_i^2*t^2),d);
    I:=ideal for i from 1 to d list i!*coefficient(t^i,series)-k_i +ideal(-1+sum for i from 1 to mix list a_i);
    eliminate(toList(mn_1..mn_mix)|toList(sd_1..sd_mix)|toList(a_1..a_mix),I)
)

--note: l_i's are actually (l_i)^(-1)
cumulantIdealExponential = method()
cumulantIdealExponential = (mix,d) -> (
    l := symbol l;
    a := symbol a;
    k := symbol k;
    R:=QQ[l_1..l_mix,a_1..a_mix,k_0..k_d][t]/t^(d+1);
    use R;
    series := formalLog(sum(apply(toList(1..mix),j->sum(apply(toList(1..d),i->a_j*l_j^i*t^i)))),d);
    --am i missing a factorial? unclear.
    I:=ideal for i from 1 to d list i!*coefficient(t^i,series)-k_i +ideal(-1+sum for i from 1 to mix list a_i);
    eliminate(toList(l_1..l_mix)|toList(a_1..a_mix),I)
)

cumulantIdealPoisson = method()
cumulantIdealPoisson = (mix,d) -> (
    l := symbol l;
    a := symbol a;
    k := symbol k;
    t := symbol t;
    R:=QQ[l_1..l_mix,a_1..a_mix,k_0..k_d][t]/t^(d+1);
    use R;
    series:=formalLog(sum for i from 1 to mix list a_i*exp(l_i*(exp(t)-1)),d);
    I:=ideal for i from 1 to d list i!*coefficient(t^i,series)-k_i+ideal(-1+sum for i from 1 to mix list a_i);
    eliminate((for i from 1 to mix list a_i)|(for i from 1 to mix list l_i),I)
)    
    
cumulantIdealMultinomial = (r,n,mix,d) -> (
    t := symbol t;
    S := QQ[t_1..t_r];
    exps := flatten apply(toList(0..d), i->flatten entries basis(i,S) / exponents / flatten);
    quotientExps := flatten entries basis(d+1,S) / exponents / flatten;
    Mons := ideal(apply(quotientExps, e->S_e));
    p := symbol p;
    a := symbol a;
    k := symbol k;
    R := QQ[p_(1,1)..p_(mix,r),a_1..a_mix,apply(exps,i->k_i)][t_1..t_r];
    Mons = sub(Mons,R);
    R = R / Mons;
    use R;
    series := formalLog(sum apply(toList(1..mix),i->a_i*(sum apply(toList(1..r), j-> p_(i,j)*exp(t_j)))^n),d); --moment gen fxn of the multinomial distribution
    I :=ideal( apply(exps, e-> (sum e)!*coefficient(sub(S_e,R),series)-k_e) ) + 
    	ideal( apply(toList(1..mix),i-> 1 - sum apply(toList(1..r), j -> p_(i,j)))) + 
	ideal(1 - sum apply(toList(1..mix),i -> a_i));
    T := QQ[apply(exps,i->k_i)];
    sub((eliminate(toList(p_(1,1)..p_(mix,r))|toList(a_1..a_mix),I),T))
)



--moment ideal of Laplace with parameters mu and b
momentIdealLaplace = method()
momentIdealLaplace (ZZ, ZZ) := Ideal => (mix,d)->(
    mn := symbol mn;
    b := symbol b;
    m := symbol m;
    t := symbol t;
    a := symbol a;
    R := symbol R;
    I := symbol I;
    if mix == 1 then(
	R=QQ[mn_1..mn_mix,b_1..b_mix,m_0..m_d][t]/t^(d+1);
	use R;
    	series:= exp(mn_1*t)*(1+sum for i from 1 to d list (b_1^2*t^2)^i);
    	I=ideal for i from 1 to d list i!*coefficient(t^i,series)-m_i;
    	return homogenize(eliminate((for i from 1 to mix list mn_i)|(for i from 1 to mix list b_i),I),m_0);
	)
    else( 
	R=QQ[mn_1..mn_mix,b_1..b_mix,a_1..a_(mix-1),m_0..m_d][t]/t^(d+1);
	use R;
    	amix := 1 - sum for i from 1 to mix-1 list a_i;
    	series2:=sum for i from 1 to mix-1 list a_i*exp(mn_i*t)*(1+sum for j from 1 to d list (b_i^2*t^2)^j) + amix*exp(mn_mix*t)*(1+sum for j from 1 to d list (b_mix^2*t^2)^j);
    	I=ideal for i from 1 to d list i!*coefficient(t^i,series2)-m_i;
    	return homogenize(eliminate((for i from 1 to mix-1 list a_i)|(for i from 1 to mix list mn_i)|(for i from 1 to mix list b_i),I),m_0)
	)
)


--cumulant ideal of Laplace with parameters mu and b
cumulantIdealLaplace = method()
cumulantIdealLaplace (ZZ, ZZ) := Ideal => (mix,d)->(
    mn := symbol mn;
    b := symbol b;
    k := symbol k;
    t := symbol t;
    a := symbol a;
    R := symbol R;
    I := symbol I;
    if mix == 1 then(
	R=QQ[mn_1..mn_mix,b_1..b_mix,k_0..k_d][t]/t^(d+1);
	use R;
    	series:= formalLog(exp(mn_1*t)*(1+sum for i from 1 to d list (b_1^2*t^2)^i),d);
    	I=ideal for i from 1 to d list i!*coefficient(t^i,series)-k_i;
    	return eliminate((for i from 1 to mix list mn_i)|(for i from 1 to mix list b_i),I);
	)
    else( 
	R=QQ[mn_1..mn_mix,b_1..b_mix,a_1..a_(mix-1),k_0..k_d][t]/t^(d+1);
	use R;
    	amix := 1 - sum for i from 1 to mix-1 list a_i;
    	series2:=sum for i from 1 to mix-1 list a_i*exp(mn_i*t)*(1+sum for j from 1 to d list (b_i^2*t^2)^j) + amix*exp(mn_mix*t)*(1+sum for j from 1 to d list (b_mix^2*t^2)^j);
    	I=ideal for i from 1 to d list i!*coefficient(t^i,series2)-k_i;
    	eliminate((for i from 1 to mix-1 list a_i)|(for i from 1 to mix list mn_i)|(for i from 1 to mix list b_i),I)
	)
)


---------------------- CONVERSIONS BETWEEN MOMENTS AND CUMULANTS ----------------------------------

-- converts a list of the first 'd' cumulants to a list of the first 'd' moments
--     'cumulants' - a list of values in the ring 'R'
--     'R' - a ring
cumulantsToMoments = (cumulants,R) -> (
    assert(cumulants#0 == 0);
    return powerSeriesTransform(cumulants, exp, R)
    )

-- converts a list of the first 'd' moments to a list of the first 'd' cumulants
--     'moments' - a list of values in the ring 'R'
--     'R' - a ring
momentsToCumulants = (moments,R) -> (
    assert(moments#0 == 1);
    return powerSeriesTransform(moments, x -> formalLog(x,#moments), R)
    )

-- helper function to do a change of variables from moments to cumulants and vice versa
powerSeriesTransform = (l, f, R) -> (
    S := R[t]/t^#l;
    use S;
    p := f(sum for i from 0 to #l-1 list t^i * l#i/i!);
    return for i from 0 to #l-1 list i!*coefficient(t^i, p)
    )

-- the first 'd' coefficients of the formal power series of the exponential function
formalExp = (f, d) -> (
    sum for k from 0 to d list 1/(k!) * f^k
)

-- the first 'd' coefficients of the formal power series of the log function
formalLog = (f, d) -> (
    sum for k from 1 to d list (-1)^(k-1)/k * (f-1)^k
)

-- converts the moment ideal 'I' to the corresponding cumulant ideal, where we only 
-- consider the first 'd' cumulants.
--     'I' - the moment ideal
--     'd' - the number of cumulants to consider
momentIdealToCumulants = (I,d) -> (
    idealTransform(I, d, formalExp, 0)  
    )

-- converts the cumulant ideal 'I' to the corresponding moment ideal, where we only 
-- consider the first 'd' moments.
--     'I' - the cumulant ideal
--     'd' - the number of moments to consider
cumulantToMoments = (I, d) -> (
    idealTransform(I, d, formalLog, 1)
    )

-- helper function that applies a change of variables to the ideal I
idealTransform = (I, truncAbove, f, substVal) -> (
    s := symbol s;
    t := symbol t;
    R2 := QQ[s_0..s_truncAbove];
    use R2;
    R2' := R2[t]/t^(truncAbove+1);
    use R2';
    p := sum for i from 0 to truncAbove list 1/i! * s_i * t^i;
    q := f(p,truncAbove);
    -- NOTE: here, we dehomogenize. May think about giving homog/nonhomog as an option.
    li := for i from 0 to truncAbove list i! * coefficient(t^i, q);
    li = for i from 0 to truncAbove list sub(li_i, s_0 => substVal);
    phi := map(R2, ring I,li);
    phi I
    )

-- helper fct to take the factorial of a multiindex
factorial = multiIndex -> fold(apply(multiIndex,k->k!),times)

-- helper fct to raise an indexed variable to a multiindex
raiseTo = (t, multiIndex) -> (
    li := for i from 1 to #multiIndex list t_i^(multiIndex_(i-1));
    fold(li, times)
)

-- converts the multivariate moment ideal 'I' to the corresponding cumulant ideal, where we only 
-- consider the cumulants described by 'cumulantsTensorFormat'. This function generalizes momentIdealToCumulants.
--     'I' - the moment ideal
--     'cumulantsTensorFormat' - a list describing how many cumulants to consider in each dimension
momentIdealToCumulantsMultivariate = (I, cumulantsTensorFormat) -> (
    idealTransformMultivariate(I, cumulantsTensorFormat, formalExp, 0)
)

-- converts the multivariate cumulant ideal 'I' to the corresponding moment ideal, where we only 
-- consider the moments described by 'momentsTensorFormat'. This function generalizes cumulantIdealToMoments.
--     'I' - the cumulant ideal
--     'momentsTensorFormat' - a list describing how many moments to consider in each dimension
cumulantIdealToMomentsMultivariate = (I, momentsTensorFormat) -> (
    idealTransformMultivariate(I, momentsTensorFormat, formalLog, 1)
)

-- helper function that applies a change of variables to the ideal I
idealTransformMultivariate = (I, tensorFormat, f, substVal) -> (
    s := symbol s;
    t := symbol t;
    zeroes := for i from 0 to #tensorFormat - 1 list 0;
    R2 := QQ[s_zeroes..s_tensorFormat];
    R2' := R2[t_1..t_#tensorFormat,MonomialOrder => Lex]/
           (for i from 1 to #tensorFormat list t_i^(tensorFormat_(i - 1)+1));
    use R2';
    p := sum for I in zeroes..tensorFormat list 1/factorial(I) * s_I * raiseTo(t,I);
    q := f(p,max(tensorFormat));
    li := for I in zeroes..tensorFormat list factorial(I) * coefficient(raiseTo(t,I), q);
    li = for x in li list sub(x, s_zeroes => substVal);
    phi := map (R2, ring I, li);
    phi I
    )

--DOCUMENTATION--

beginDocumentation()

doc ///
  Key
     MomentVariety
  Headline
     A package for computing the moments of distributions 
  Description
   Text
    {\em MomentVariety} is a package that computes the moments of some distributions and finds the underlying ideals.

///    

doc ///
  Key    
    GaussianMoments
    (GaussianMoments, ZZ, Ring)
  Headline
    lists all moments of the univariate Gaussian distribution
  Usage 
    GaussianMoments (n,R) 
  Inputs
    n : ZZ
    R : Ring
  Outputs
    : BasicList
  Description
    Text
      GaussianMoments computes and lists all the moments of the univariate Gaussian distribution.
    Text
      Here we show an example.
    Example
      R = QQ[x_0..x_3]
      d = 2
      GaussianMoments (d,R)

///

doc ///
  Key
    momentIdealExponential
    (momentIdealExponential, ZZ, ZZ, Ring)
  Headline
    compute the ideal corresponding to the Exponential
  Usage
    momentIdealExponential (d,mix)
  Inputs
    d : ZZ
    mix : ZZ
  Outputs
    : Ideal
  Description
    Text
      given the highest degree of the moments and the number of mixtures compute the ideal of the Exponential distribution
    Text
      Here we show an example
    -- Example
    --   MISSING
    
///

doc ///
  Key
    momentIdealGaussian
    (momentIdealGaussian, ZZ, ZZ, Ring)
  Headline
    computes the homogeneous moment ideal of the Gaussian
  Usage
    I = momentIdealGaussian (mix,d, K)
  Inputs
    mix : ZZ
    d : ZZ
  Outputs
    I : Ideal
  Description
    Text
      given the number of mixtures and the highest degree of moments, compute the corresponding homogeneous ideal
    Text
      Here we show an example
    Example
      mix = 2
      d = 1
      momentIdealGaussian (mix,d,K)
      
///

doc ///
  Key
    momentIdealPoisson 
    (momentIdealPoisson, ZZ, ZZ)
  Headline
    compute the homogeneous moment ideal
  Usage
    I = momentIdealPoisson(mix,d,K)
  Inputs
    mix : ZZ
    d : ZZ
    K: Ring
  Outputs
    I : Ideal
  Description
    Text
      given the number of mixtures and the hightest degree of moments, compute the corresponding homogeneous moment ideal
    Text
      Here we show an example
    Example
      MISSING
      
 ///
 
 doc ///
   Key
     momentIdealGaussianTest
     (momentIdealGaussianTest, ZZ, ZZ, Ring)
   Headline
     TO BE GIVEN
   Usage
     I = momentIdealGaussianTest(mix,d,K)
   Inputs
     mix : ZZ
     d : ZZ
   Outputs
     I : Ideal
   Description
     Text
       TO BE GIVEN
     Text
       Here we show an example
     -- Example
     --   TO BE GIVEN
       
 ///
 
 doc ///
   Key
     momentVarietyGaussians
     (momentVarietyGaussians, ZZ, ZZ)
   Headline
     compute the homogeneous ideal of the moment variety
   Usage
     momentVarietyGaussians (n,d)
   Inputs
     n : ZZ
     d : ZZ
   Outputs
     : Ideal
   Description
     Text
       compute the homogeneous ideal of the moment variety
     Text
       Here we show an example
     Example
       n = 1
       d = 4
       momentVarietyGaussians (n,d)
 
 ///
 

 doc ///
   Key
     momentIdealMultinomial
     (momentIdealMultinomial, ZZ, ZZ, ZZ, Ring)
   Headline
     multinomial distribution
   Usage
     I = momentIdealMultinomial (k, n, d, K)
   Inputs
     k : ZZ
     n : ZZ
     d : ZZ
     K : Ring
   Outputs 
     I : Ideal
   Description
     Text
       Given the number of possible outcomes, the number of trials in a statistical experiment as well as the truncation order, compute the moment ideal for the multinomial distribution      
     Text
       Here we show an example
     Example
       k = 2
       n = 3
       d = 2
       momentIdealMultinomial (k,n,d) 
 
 ///
 
end--
uninstallPackage "MomentVariety"
restart
installPackage "MomentVariety"



