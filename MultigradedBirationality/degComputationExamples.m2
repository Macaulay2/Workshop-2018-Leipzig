restart
loadPackage "Cremona" --- Cremona map in Macaulay2
loadPackage "MultiGradedRationalMap" -- new criterion

installPackage "MultiGradedRationalMap"
------ PURPOSE OF THIS FILE-------
-- In this given a long list of examples, where we certify our computations
--   with the package "Cremona" which is provided in Macaualay2.
--
-- For most of these examples we call the function for computing the degree of a rational 
--   map in the "Cremona" package, then we call the functions provided in our package to 
--   compute the degree and check birationality.

---******************----------------------------------
-- Example 1 (some examples without base points)
R = QQ[x,y,z]
I = ideal(random(2, R), random(2, R), random(2, R))
-- call functions of the package ``Cremona''
degree rationalMap gens I 
-- call our implementations
degreeOfMap I 
degreeOfMap(I, Strategy=>SatSpecialFibStrategy)
degreeOfMapIter(I, 5)
isBiratMap I
jDRank I

I = ideal(random(3, R), random(3, R), random(3, R))
degree rationalMap gens I 
degreeOfMap I
degreeOfMap(I, Strategy=>SatSpecialFibStrategy)
degreeOfMapIter(I, 5)
isBiratMap I 
jDRank I


I = ideal(random(4, R), random(4, R), random(4, R))
degree rationalMap gens I 
degreeOfMap I
degreeOfMap(I, Strategy=>SatSpecialFibStrategy)
degreeOfMapIter(I, 5)
isBiratMap I 
jDRank I


I = ideal(random(5, R), random(5, R), random(5, R))
degree rationalMap gens I 
degreeOfMap I
degreeOfMap(I, Strategy=>SatSpecialFibStrategy)
degreeOfMapIter(I, 5)
isBiratMap I 
jDRank I


I = ideal(random(6, R), random(6, R), random(6, R))
degree rationalMap gens I 
degreeOfMap I
degreeOfMap(I, Strategy=>SatSpecialFibStrategy)
degreeOfMapIter(I, 5)
isBiratMap I 
jDRank I


-----*****************---------------------------
-- Example 2 (some examples by playing with our characterization of Theorem 3.17)
R = QQ[x,y,z]
A = matrix{ {x, x^2 + y^2},
            {-y, y^2 + z*x},
	    {0, x^2}
	   }
I = minors(2, A) -- a birational map
degree rationalMap gens I 
degreeOfMap I 
degreeOfMap(I, Strategy=>SatSpecialFibStrategy)
degreeOfMapIter(I, 5)
isBiratMap I
jDRank I
satSpecialFiber I


A = matrix{ {x, x^5 + y^5},
            {-y, y^5 + z*x^2*y^2},
	    {0, x^5}
	   }
I = minors(2, A) -- a birational map
degree rationalMap gens I 
degreeOfMap I 
degreeOfMap(I, Strategy=>SatSpecialFibStrategy)
degreeOfMapIter(I, 5)
isBiratMap I
jDRank I
satSpecialFiber I


A = matrix{ {x, x^6 + y^6 + z*x^5},
            {-y, y^6 + z*x^3*y^2},
	    {0, x^6 + x*y^4*z}
	   }
I = minors(2, A) -- a birational map
degree rationalMap gens I 
degreeOfMap I 
degreeOfMap(I, Strategy=>SatSpecialFibStrategy)
degreeOfMapIter(I, 5)
isBiratMap I
jDRank I
satSpecialFiber I


A = matrix{ {x^2, x^2 + y^2},
            {-y^2, y^2 + z*x},
	    {0, x^2}
	   }
I = minors(2, A) -- a non birational map
degree rationalMap gens I 
degreeOfMap I 
degreeOfMap(I, Strategy=>SatSpecialFibStrategy)
degreeOfMapIter(I, 5)
isBiratMap I
jDRank I
satSpecialFiber I


A = matrix{ {x^3, x^2 + y^2},
            {-y^3, y^2 + z*x},
	    {0, x^2}
	   }
I = minors(2, A) -- a non birational map
gensSatSpecialFib I
degree rationalMap gens I 
degreeOfMap I 
degreeOfMap(I, Strategy=>SatSpecialFibStrategy)
degreeOfMapIter(I, 5)
isBiratMap I
jDRank I
satSpecialFiber I
gensSatSpecialFib I


A = matrix{ {x^3, x^2 + y^2},
            {-y^3, y^2 + z*x},
	    {0, x^2}
	   }
I = minors(2, A) -- a non birational map
gensSatSpecialFib I
degree rationalMap gens I 
degreeOfMap I 
degreeOfMap(I, Strategy=>SatSpecialFibStrategy)
degreeOfMapIter(I, 5)
isBiratMap I
jDRank I
satSpecialFiber I
gensSatSpecialFib I



A = matrix{ {x^3, x^4},
            {-y^3, y^4},
	    {z^3, x^4}
	   }
I = minors(2, A) -- a non birational map
degree rationalMap gens I 
degreeOfMap I 
isBiratMap I
jDRank I


A = matrix{ {x^3, random(1, R)},
            {-y^3, y},
	    {z^3, random(1, R)}
	   }
I = minors(2, A) -- a non birational map
degree rationalMap gens I 
degreeOfMap I 
isBiratMap I
jDRank I


A = matrix{ {random(1, R), random(1, R)},
            {random(1, R), random(1, R)},
	    {random(1, R), random(1, R)}
	   }
I = minors(2, A) -- a birational map
degree rationalMap gens I 
degreeOfMap I 
isBiratMap I
jDRank I



-------------**************------------
------Example 3 (some examples of a map PP^2 --> PP^3)
R = QQ[x,y,z]

I = ideal(random(2,R), random(2, R), random(2, R), random(2, R)) -- birational
degree rationalMap gens I 
degreeOfMap I 
degreeOfMap(I, Strategy=>SatSpecialFibStrategy)
degreeOfMapIter(I, 5)
isBiratMap I
jDRank I

I = ideal(x^4, y^4, z^4, x*y*z^2) -- non birational
degree rationalMap gens I 
degreeOfMap I 
degreeOfMap(I, Strategy=>SatSpecialFibStrategy)
degreeOfMapIter(I, 5) 
isBiratMap I
jDRank I


I = ideal(7*x^3, 2*y^3 + 11*x*z^2, z^3 - z*x^2, x*z^2) -- non birational
degree rationalMap gens I 
degreeOfMap I 
degreeOfMap(I, Strategy=>SatSpecialFibStrategy)
degreeOfMapIter(I, 5) 
isBiratMap I
jDRank I


I = ideal(7*x^4 + x^2*y^2, 13*y^4 + 11*x*y*z^2, z^4 - z^2*x^2, x*z^3) -- birational
degree rationalMap gens I 
degreeOfMap I 
degreeOfMap(I, Strategy=>SatSpecialFibStrategy)
degreeOfMapIter(I, 5) 
isBiratMap I
jDRank I



-------------****************************-------------
---- Example 4 (non generically finite maps)
R = QQ[x,y,z]

I = ideal(random(4, R), random(4, R))
degree rationalMap gens I 
degreeOfMap I 
isBiratMap I
jDRank I

I = ideal(x^2 + z^2, y^2 + z^2, x^2 - y^2)
degree rationalMap gens I 
degreeOfMap I 
isBiratMap I
jDRank I



-------------**************------------
------Example 5 (some examples of a map PP^2 --> PP^4)
R = QQ[x,y,z]

I = ideal(random(2,R), random(2, R), random(2, R), random(2, R), random(2, R)) -- birational
degree rationalMap gens I 
degreeOfMap I 
degreeOfMapIter(I, 5) 
isBiratMap I
jDRank I


I = ideal(x^7, z^7 , y^7, x^7, x^2*y^2*z^3) -- non birational
degree rationalMap gens I 
degreeOfMap I 
degreeOfMap(I, Strategy=>SatSpecialFibStrategy)
degreeOfMapIter(I, 5) 
isBiratMap I
jDRank I



-------------****************************-------------
---- Example 4 (non generically finite maps)
R = QQ[x,y,z]

I = ideal(random(4, R), random(4, R))
degree rationalMap gens I 
degreeOfMap I 
isBiratMap I
jDRank I

I = ideal(x^2 + z^2, y^2 + z^2, x^2 - y^2)
degree rationalMap gens I 
degreeOfMap I 
isBiratMap I
jDRank I



-------------**************------------
------Example 5 (some examples of a map PP^2 --> PP^4)
R = QQ[x,y,z]

I = ideal(random(2,R), random(2, R), random(2, R), random(2, R), random(2, R)) -- birational
degree rationalMap gens I 
degreeOfMap I 
degreeOfMapIter(I, 5) 
isBiratMap I
jDRank I


I = ideal(x^7, z^7 , y^7, x^7, x^2*y^2*z^3) -- non birational
degree rationalMap gens I 
degreeOfMap I 
degreeOfMap(I, Strategy=>SatSpecialFibStrategy)
degreeOfMapIter(I, 5) 
isBiratMap I
jDRank I


-------------**************------------
------Example 6 (some examples of a map PP^3 --> PP^3)
R = QQ[x,y,z,w]

A = matrix{ {x + y,  x, x},
            {3*z - 4*w, y, z},
	    {w,  z, z + w}, 
	    {y - z,  w, x + y}
	   }
I = minors(3, A) -- birational
degree rationalMap gens I 
degreeOfMap I
satSpecialFiber I
degreeOfMap(I, Strategy=>SatSpecialFibStrategy)
degreeOfMapIter(I, 5) 
isBiratMap I
jDRank I

A = matrix{ {x,  x, x^2},
            {z - w, y, z^2},
	    {w,  z, z*w}, 
	    {y - z,  w, x*y}
	   }
I = minors(3, A) -- non birational
degree rationalMap gens I 
degreeOfMap I 
isBiratMap I
jDRank I


I = ideal(random(2, R), random(2, R), random(2, R), random(2, R))
degree rationalMap gens I 
degreeOfMap I 
isBiratMap I
jDRank I
upperBoundDegreeSingleGraded I


I = ideal(x^5, y^5, z^5, w^5)
degree rationalMap gens I 
degreeOfMap I 
isBiratMap I
jDRank I
upperBoundDegreeSingleGraded I

I = ideal(random(2, R), x^2, y^2, z*w)
degree rationalMap gens I 
degreeOfMap I 
isBiratMap I
jDRank I
upperBoundDegreeSingleGraded I

I = ideal(x^4 + x*y^3, w^4, z^4, x*y*z*w)
degree rationalMap gens I 
degreeOfMap I 
isBiratMap I
jDRank I
upperBoundDegreeSingleGraded I


I = ideal(13*x^4 + x*y^3, w^4, z*y^3, x*y*z*w)
degree rationalMap gens I 
degreeOfMap I 
isBiratMap I
jDRank I
upperBoundDegreeSingleGraded I


I = ideal(13*x^4 + x*y^3, 7*w^2*y^2, 11*z*y^3, 7*x*z^3)
degree rationalMap gens I 
degreeOfMap I 
isBiratMap I
jDRank I


I = ideal(random(4, R), 7*w^2*y^2, 11*z*y^3, 7*x*z^3)
degree rationalMap gens I 
degreeOfMap I 
isBiratMap I
jDRank I
upperBoundDegreeSingleGraded I


I = ideal(13*x^4 + 5*y^4, 7*w^2*y^2, 11*z*y^3, 7*x*z^3)
degree rationalMap gens I 
degreeOfMap I 
isBiratMap I
jDRank I


I = ideal(17*z^5 + 11*w^3*x^2 + y^5, 3*w^5 + x^5 + 7*w^5, 2*x^2*y^3 + 13*w*z^4, x^3*z^2 )
degree rationalMap gens I 
degreeOfMap I 
isBiratMap I
jDRank I
upperBoundDegreeSingleGraded I




-------------**************------------
------Example 7 (some examples of a map PP^4 --> PP^4)
R = QQ[x,y,z,v,w]

I = ideal(random(1, R), random(1, R), random(1, R), random(1, R), random(1, R))
degree rationalMap gens I 
degreeOfMap I 
degreeOfMap(I, Strategy=>SatSpecialFibStrategy)
degreeOfMapIter(I, 5)
isBiratMap I
jDRank I
upperBoundDegreeSingleGraded I

I = ideal(29*x^3 + 55*x*y*z, 7*y^3, 14*z^3, 17*v^3, 12*w^3)
degree rationalMap gens I 
degreeOfMap I 
degreeOfMap(I, Strategy=>SatSpecialFibStrategy)
degreeOfMapIter(I, 5)
isBiratMap I
jDRank I
upperBoundDegreeSingleGraded I

I = ideal(x^5, y^5, z^5, v^5, w^5)
degree rationalMap gens I 
degreeOfMap I 
isBiratMap I
jDRank I
upperBoundDegreeSingleGraded I

I = ideal(x^3*y, y^4 + z^4, x*y*z*w, v^4, y^2*z^2 - v^2*w^2)
degree rationalMap gens I 
degreeOfMap I 
isBiratMap I
jDRank I

I = ideal(14*x*y + v^2, 91*y^2 - 4*z^2, 19*z*w, 17*v^2, 31*y^2 - 23*v*w)
degree rationalMap gens I 
degreeOfMap I 
isBiratMap I
jDRank I

I = ideal(x^2*y*z, y^2*z^2, z^2*v^2, v^2*w^2, w^2*x^2)
degree rationalMap gens I 
degreeOfMap I 
isBiratMap I
jDRank I
 
I = ideal(random(2, R), x^2, y^2, z^2, v^2)
degree rationalMap gens I 
degreeOfMap I 
isBiratMap I
jDRank I
upperBoundDegreeSingleGraded I



-- Parametrization of plane curves
R = QQ[x,y]

A = matrix{ {x, x^5 + y^5},
            {-y, y^5 + x^3*y^2},
	    {0, x^5}
	   }
I = minors(2, A)
degree rationalMap gens I 
degreeOfMap I 
degreeOfMapIter(I, 5) 
isBiratMap I
jDRank I


A = matrix{ {x^3, x^6 + y^6},
            {-y^3, y^6 + x^3*y^3},
	    {0, x^6}
	   }
I = minors(2, A)
degree rationalMap gens I 
degreeOfMap I 
degreeOfMapIter(I, 5) 
isBiratMap I
jDRank I



--- Bigraded rational maps
R = QQ[x,y,u,v, Degrees => {{1,0}, {1,0}, {0,1}, {0,1}}]

I = ideal(x*u, y*u, y*v) -- birational map
degreeOfMapIter(I, 5)
isBiratMap I
jDRank I
partialJDRs I
satSpecialFiber(I, 5)


I = ideal(x*u, y*v, x*v + y*u) -- non birational map
gensSatSpecialFib(I, 5)
degreeOfMapIter(I, 5)
isBiratMap I
jDRank I
partialJDRs I

I = ideal(x*u^2, y*u^2, x*v^2) -- non birational map
degreeOfMapIter(I, 5)
isBiratMap I
jDRank I 
partialJDRs I


A = matrix{ {random({1,0}, R), random({0,1}, R)},
            {random({1,0}, R), random({0,1}, R)},
            {random({1,0}, R), random({0,1}, R)}
          }
I = minors(2, A) -- birational map
degreeOfMapIter(I, 5)
isBiratMap I
jDRank I
partialJDRs I


A = matrix{ {random({0,1}, R), random({1,7}, R)},
            {random({0,1}, R), random({1,7}, R)},
            {random({0,1}, R), random({1,7}, R)}
          }
I = minors(2, A) -- birational map
isBiratMap I
jDRank I
partialJDRs I


A = matrix{ {random({0,2}, R), random({1,4}, R)},
            {random({0,2}, R), random({1,4}, R)},
            {random({0,2}, R), random({1,4}, R)}
          }
I = minors(2, A) -- non birational map
isBiratMap I
jDRank I
partialJDRs I

A = matrix{ {x*u + y*v, y*v},
    	    {y*u - (x+y)*v, y*(u + v)},
	    {12*x*v - 5*y*u, (x-4*y)*u}
    	  }
I = minors(2, A)  -- non birational
isBiratMap I
jDRank I
partialJDRs I


A = matrix{ {x^2*u,  x^2*v^2},
    	    {y^2*v, x^2*u^2},
	    {0,     y^2*v^2}
    	  }
I = minors(2, A)  -- non birational
degreeOfMapIter(I, 5)
isBiratMap I
jDRank I
partialJDRs I

A = matrix{ {x^5*u,  x^2*v^2},
    	    {y^5*v, x^2*u^2},
	    {0,     y^2*v^2}
    	  }
I = minors(2, A)  -- non birational
degreeOfMapIter(I, 5)
isBiratMap I
jDRank I
partialJDRs I



