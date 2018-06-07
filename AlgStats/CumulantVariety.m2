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