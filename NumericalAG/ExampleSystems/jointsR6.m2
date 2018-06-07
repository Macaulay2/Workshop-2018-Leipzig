export { "jointsR6" }

jointsR6 = method()
jointsR6 Ring := kk -> (
    -- mixed_volume=64, nreal_solutions=10, 
    -- nsolutions=32 (Julia stats; is nsolutions correct?)
    x := symbol x;
    R := kk[x_1..x_8];
    { 
	x_1^2+x_2^2-1,
        x_3^2+x_4^2-1,
        x_5^2+x_6^2-1,
        x_7^2+x_8^2-1,
        (-2.4915068e-01*x_1*x_3+ 1.6091354e00*x_1*x_4+ 2.7942343e-01 *x_2*x_3+ 1.4348016e00*x_2*x_4) + (4.0026384e-01*x_5*x_8-8.0052768e-01*x_6*x_7+ 7.4052388e-02*x_1-8.3050031e-02*x_2) - (3.8615961e-01*x_3-7.5526603e-01*x_4+ 5.0420168e-01*x_5 -1.0916287e00*x_6+ 4.0026384e-01*x_8) + 4.920729e-02,
        (1.2501635e-01*x_1*x_3-6.8660736e-01*x_1*x_4-1.1922812e-01* x_2*x_3-7.1994047e-01*x_2*x_4) - (4.3241927e-01*x_5*x_7-8.6483855e-01*x_6*x_8-3.715727e-02*x_1+ 3.5436896e-02*x_2)+ 8.5383482e-02*x_3-3.9251967e-02*x_5-4.3241927e-01*x_7+ 1.387301e-02,
        (-6.3555007e-01*x_1*x_3-1.1571992e-01*x_1*x_4-6.6640448e-01 *x_2*x_3) + (1.1036211e-01*x_2*x_4+ 2.9070203e-01*x_5*x_7+ 1.2587767e00*x_5*x_8)- (6.2938836e-01*x_6*x_7+ 5.8140406e-01*x_6*x_8+ 1.9594662e-01*x_1)- (1.2280342e00*x_2-7.9034221e-02*x_4+ 2.6387877e-02*x_5)- 5.713143e-02*x_6-1.1628081e00*x_7+1.2587767e00*x_8+ 2.162575e00,
        (1.4894773e00*x_1*x_3+ 2.3062341e-01*x_1*x_4+ 1.3281073e00*x_2*x_3)-(2.5864503e-01*x_2*x_4+ 1.165172e00*x_5*x_7-2.6908494e-01*x_5*x_8)+ (5.3816987e-01*x_6*x_7+ 5.8258598e-01*x_6*x_8-2.0816985e-01*x_1)+(2.686832e00*x_2-6.9910317e-01*x_3+ 3.5744413e-01*x_4)+ 1.2499117e00*x_5+ 1.467736e00*x_6+ 1.165172e00*x_7+ 1.10763397e00*x_8-6.9686809e-01
    }
    )

doc /// 
    Key
    	jointsR6
	(jointsR6,Ring)
    Headline
    	the six-revolute-joint problem
    Description
    	Text
	  A. Morgan and A. Sommese
	  `Computing all solutions to polynomial systems using homotopy continuation',
 	  Appl. Math. Comput., Vol. 24, pp 115-138, 1987.
	Example
	    F = jointsR6 RR_53
	    sols = solveSystem F;
	    #sols
	    assert(#realPoints sols == 10)
	    needsPackage "MonodromySolver"
	    sols = sparseMonodromySolve polySystem F;
	    assert(#realPoints sols == 10)	    
    ///

