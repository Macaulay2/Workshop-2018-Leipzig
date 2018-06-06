restart
cumulantsToMoments = (cumulants,R) -> (
    assert(cumulants#0 == 0);
    return powerSeriesTransform(cumulants, exp, R)
    )

momentsToCumulants = (moments,R) -> (
    assert(moments#0 == 1);
    return powerSeriesTransform(moments, formalLog, R)
    )

powerSeriesTransform = (l, f, R) -> (
    S := R[t]/t^#l;
    p := f(sum for i from 0 to #l-1 list t^i * l#i/i!);
    return for i from 0 to #l-1 list i!*coefficient(t^i, p)
    )

formalLogOld = f -> (
    assert ((reverse flatten entries last coefficients(f))_0 == 1);
    d := #(coefficients f)_0;
    sum for i from 1 to d list (-1)^(i-1)/i * (f-1)^i
    )


cTM = d -> (
    R:=QQ[k_1..k_d][t]/t^(d+1);
    cumulantsum := sum for i from 1 to d list (k_i*t^i/(i!));
    expcumulantsum :=  1+sum for i from 1 to d list (cumulantsum^i/(i!));
    d!*coefficient(t^d,expcumulantsum)
    )

formalExp = (f, d) -> (
    sum for k from 0 to d-1 list 1/(k!) * f^k
)

formalLog = (f, d) -> (
    sum for k from 1 to d-1 list (-1)^(k-1)/k * (f-1)^k
)

-- TEST --

truncAbove = 6
I = momentIdealGaussian(1,truncAbove)
R2 = QQ[k_0..k_truncAbove]
R2' = R2[t]/t^(truncAbove+1)
p = sum for i from 0 to truncAbove list 1/i! * k_i * t^i
q = formalExp(p,truncAbove+1)
li = for i from 0 to truncAbove list i! * coefficient(t^i, q)
li = for i from 0 to truncAbove list sub(li_i,{k_0 => 0})
phi = map(R2, QQ[m_0..m_truncAbove],li)
I = sub(I, source phi)
J = phi I
transpose gens gb J


