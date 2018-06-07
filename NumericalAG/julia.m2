-*
export{"JuliaProcess", "solveJulia", "restartJulia"}
*-

-------------------
---GLOBAL VARIABLES
------------------

juliaBinary="/home/tim/julia-d55cadc350/bin/julia"
JuliaProcess=openInOut("!" |juliaBinary)
noProcessError="JuliaProcess not currently open"


---------------------
--- PACKAGE IMPORTS
---------------------

needsPackage "NumericalAlgebraicGeometry"

--------------
---UNEXPORTED-
--------------

--write exported wrapper for first two
writeSys = method(Options=>{IncludeTemplate=>false})
writeSys (PolySystem, File) := o -> (P,f) -> (
    R:=ring P;
    varString:=apply(gens R,g->(toString g));
    varCommas:=(P.NumberOfVariables-1):", ";
    eqnCommas:=(P.NumberOfPolys-1):", ";    
    f << concatenate(mingle(varString,varCommas))| " = " | "[PolyVar{true}(i) for i in [" | concatenate mingle(apply(varString,g-> "\"" | g | "\""),varCommas)|"]];\n";
    f <<"f = [" | concatenate mingle(
	apply(equations P,e->replace("p[0-9]*e","e",toExternalString e)),eqnCommas)| "];\n";
    )
writeSys (PolySystem, String) :=o -> (P,filename) -> (
    f := openOut filename;
    writeSys(P,f);
    close f;
    )
writeSys PolySystem := o -> P -> (
    if not isOpen JuliaProcess then error("JuliaProcess not open");
    writeSys(P,JuliaProcess)
    )


-------------------------------------
--- EXPORTED FUNCTIONS AND METHODS---
-------------------------------------

restartJulia = () -> (
    if (isOpen JuliaProcess) then close JuliaProcess;
    JuliaProcess=openInOut("!" |juliaBinary);
    )


TEST ///
JuliaProcess
close JuliaProcess
restartJulia()
isOpen JuliaProcess
///

importJulia = importList -> JuliaProcess<<concatenate apply(importList,pkg->pkg | "\n")

-- todo: replace arbitrary variable name "f" for julia system
solveJulia = method(Options=>{})
solveJulia PolySystem := o -> P -> (
    if not isOpen JuliaProcess then error(noProcessError);
--    print "Wait for system to be solved:\n";
    writeSys P;
    JuliaProcess<<flush;
    x:=read JuliaProcess;
    R:=ring P;
    varString:=apply(gens R,g->toExternalString g);
    varCommas=toList((P.NumberOfVariables-1):", ");
    polyCommas=toList((P.NumberOfPolys-1):", ");
    JuliaProcess<<concatenate(mingle(varString,varCommas))| " = " | "[PolyVar{true}(i) for i in [" | concatenate(mingle(apply(varString,g-> "\"" | g | "\""),varCommas)) |"]];\n";
    JuliaProcess<<"f=[" | concatenate mingle(
	apply(equations P,e->replace("p[0-9]+e","e",toExternalString e)),
        polyCommas) | "];\n";
    JuliaProcess<<"sols=solve(f)\n";
    JuliaProcess<<"x=[s.solution for s in sols];\n";
    JuliaProcess<<"show(IOContext(STDOUT, :compact=>false),[s.solution for s in sols])\n";
    finished = "done!!!x";
    JuliaProcess<<"print(\"" | finished |"\")\n";
    out="";
    while not match(finished,out) do (
	JuliaProcess<<flush;
    	out = out | read JuliaProcess;	
	);
    L=select("\\[([0-9]|| |-|\\.|e|i|m|\\+|,)*\\]",out);
    sols=apply(L,s->point{apply(separate(",",replace("\\[|\\]","",replace("im","*ii",s))),x->value x)})
    )    
    

-*
finite(sols)
solutions(sols)
results(sols)
*-



end

restart
needs "julia.m2"
importJulia {"using HomotopyContinuation","import DynamicPolynomials: PolyVar"}

R=CC[x,y]
f= x^2+y
g=y^2-pi
P=polySystem {f,g}
solveJulia P

