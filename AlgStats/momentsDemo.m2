restart
load "MomentVariety.m2"
load "MomentsToCumulants.m2"
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
I6 = momentIdealMultiMixture(3,10,2,2)

R =  QQ[mu,si,t]/t^5
I7 = momentIdealFromMGF(1, 4, {exp(t*mu_1 + 1/2 * si_1^2 * t^2),{mu_1,si_1}}, R)

-- We can also find relations between the cumulants
I8 = cumulantIdealGaussian(1,4)

-- Convert moment coordinates into cumulant coordinates
I9 = momentIdealToCumulants(I1, 4)
transpose gens gb I9

I10 = momentIdealToCumulants(I9,2)
