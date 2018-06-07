restart
load "MomentVariety.m2"


-- All relations between the first three moments
-- of a univariate Gaussian
I1 = momentIdealGaussianTest(1,4)

-- All relation between the first three moments
-- of a mixture of two Gaussians
I2 = momentVarietyGaussians(2,3)

-- These are homogenized, i.e. the zeroeth moment
-- appears as a variable

-- Poisson, first two moments
I3 = momentIdealPoisson(1,2)

-- Mixture of two Poissons, first four moments
I4 = momentIdealPoisson(2,4)

-- Multinomial distribution with three outcomes,
-- ten trials, two moments
I5 = momentIdealMultinomial(3,10,2)

-- Mixture of two such multinomial distributions
I6 = momentIdealMultiMixture(3,10,2,Mixture => 2)

-- R =  QQ[mu,si,t]/t^5
-- I7 = momentIdealFromMGF(1, 4, {exp(t*mu_1 + 1/2 * si_1^2 * t^2),{mu_1,si_1}}, R)

-- We can also find relations between the cumulants
I7 = cumulantIdealGaussian(1,4)

-- Convert moment coordinates into cumulant coordinates
I9 = momentIdealToCumulants(I1, 4)
transpose gens gb I9

-- Let's try another one
I10 = momentIdealToCumulants(I3,2)

-- These operations are mutually inverse, after
-- dehomogenizing
cumulantToMoments(I10,2)
I3

-- We can also transform a list of moments into
-- cumulants and vice versa.
R = QQ[a,b,c]
li = cumulantsToMoments({0,a,b,c,0},R)
li2 = momentsToCumulants(li, R)

-- And now, conversion between multivariate moments
-- and cumulants
I11 = momentVarietyGaussians(2,3)
I12 = momentIdealToCumulantsMultivariate(I11,{3,3})
transpose gens gb I12
