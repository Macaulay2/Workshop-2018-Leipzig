cumulantsToMoments = (cumulants,R) -> (
    assert(cumulants#0 == 0);
    return powerSeriesTransform(cumulants, exp, R)
    )

momentsToCumulants = (moments,R) -> (
    assert(moments#0 == 1);
    return powerSeriesTransform(moments, x -> formalLog(x,#moments), R)
    )

powerSeriesTransform = (l, f, R) -> (
    S := R[t]/t^#l;
    use S;
    p := f(sum for i from 0 to #l-1 list t^i * l#i/i!);
    return for i from 0 to #l-1 list i!*coefficient(t^i, p)
    )

formalExp = (f, d) -> (
    sum for k from 0 to d list 1/(k!) * f^k
)

formalLog = (f, d) -> (
    sum for k from 1 to d list (-1)^(k-1)/k * (f-1)^k
)

momentIdealToCumulants = (I,truncAbove) -> (
    k := symbol k;
    t := symbol t;
    R2 := QQ[k_0..k_truncAbove];
    use R2;
    R2' := R2[t]/t^(truncAbove+1);
    use R2';
    p := sum for i from 0 to truncAbove list 1/i! * k_i * t^i;
    q := formalExp(p,truncAbove+1);
    -- NOTE: here, we dehomogenize. May think about giving homog/nonhomog as an option.
    li := for i from 0 to truncAbove list i! * coefficient(t^i, q);
    li = for i from 0 to truncAbove list sub(li_i, k_0 => 0);
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

-- momentsTensorFormat should contain the indices that are _included_ (generalizes truncAbove) 
momentIdealToCumulantsMultivariate = (I, truncAbove, momentsTensorFormat) -> (
    return 0
)

end


-- TEST --

truncAbove = 4
momentsTensorFormat = {truncAbove, truncAbove, truncAbove}

momentsTensorFormat = {3,3,3}
apply(momentsTensorFormat, x -> x - 1)
t = symbol t
zeroes = for i from 0 to #momentsTensorFormat - 1 list 0
R2 = QQ[k_zeroes..k_momentsTensorFormat]
R2' = R2[t_1..t_#momentsTensorFormat,MonomialOrder => Lex]/
           (for i from 1 to #momentsTensorFormat list t_i^(momentsTensorFormat_(i - 1)+1))
peek R2'
use R2'
p = sum for I in zeroes..momentsTensorFormat list 1/factorial(I) * k_I * raiseTo(t,I)
q = formalExp(p,max(momentsTensorFormat))
coefficients q
