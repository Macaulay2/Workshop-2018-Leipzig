-- -*- coding: utf-8 -*-
newPackage(
        "MomentIdeals",
        Version => "0.1", 
        Date => "June 8th, 2018",
        Authors => {
	             {Name => "Jane Doe", 
                      Email => "doe@math.uiuc.edu", 
                      HomePage => "https://faculty.math.illinois.edu/~doe/"},
 	             {Name => "Jane Doe", 
                      Email => "doe@math.uiuc.edu", 
                      HomePage => "https://faculty.math.illinois.edu/~doe/"},
		     {Name => "Jane Doe", 
                      Email => "doe@math.uiuc.edu", 
                      HomePage => "https://faculty.math.illinois.edu/~doe/"},
		     {Name => "Jane Doe", 
                      Email => "doe@math.uiuc.edu", 
                      HomePage => "https://faculty.math.illinois.edu/~doe/"},
	 	     {Name => "Jane Doe", 
                      Email => "doe@math.uiuc.edu", 
                      HomePage => "https://faculty.math.illinois.edu/~doe/"},
		     {Name => "Jane Doe", 
                      Email => "doe@math.uiuc.edu", 
                      HomePage => "https://faculty.math.illinois.edu/~doe/"},
		     {Name => "Jane Doe", 
                      Email => "doe@math.uiuc.edu", 
                      HomePage => "https://faculty.math.illinois.edu/~doe/"}
		  
	      },
        Headline => "A Macaulay2 package to compute moment and cumulant ideals",
        DebuggingMode => true
        )

export {"firstFunction"}

firstFunction = method(TypicalValue => String)
firstFunction ZZ := String => n -> if n == 1 then "Hello World!" else "D'oh!"

beginDocumentation()
multidoc ///
 Node
  Key
   FirstPackage
  Headline
     an example Macaulay2 package
  Description
   Text
    {\em FirstPackage} is a basic package to be used as an example.
  Caveat
    Still trying to figure this out.
 Node
  Key
   (firstFunction,ZZ)
   firstFunction
  Headline
   a silly first function
  Usage
   firstFunction n
  Inputs
   n:
  Outputs
   :
    a silly string, depending on the value of {\tt n}
  Description
   Text
    Here we show an example.
   Example
    firstFunction 1
    firstFunction 0
///

TEST ///
    assert ( firstFunction 2 == "D'oh!" )
///

end