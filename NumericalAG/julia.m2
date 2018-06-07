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

parsePolynomial = p -> replace("p[0-9]+","0",toExternalString p)--replace("p[0-9]+(\\*|\\}|,)","*",replace("p[0-9]+e","e",toExternalString p))

--todo: write exported wrapper for first two signatures
writeSys = method(Options=>{IncludeTemplate=>false,WithImports=>false,ImportList=>{"using HomotopyContinuation","import DynamicPolynomials: PolyVar"}})
writeSys (PolySystem, File) := o -> (P,f) -> (
    R:=ring P;
    varString:=apply(gens R,g->(toString g));
    varCommas:=(P.NumberOfVariables-1):", ";
    eqnCommas:=(P.NumberOfPolys-1):", ";
    if o.WithImports then importJulia(o.ImportList,f);
    f << concatenate(mingle(varString,varCommas))| " = " | "[PolyVar{true}(i) for i in [" | concatenate mingle(apply(varString,g-> "\"" | g | "\""),varCommas)|"]];\n";
    f <<"f = [" | concatenate mingle(
	apply(equations P,e->parsePolynomial e),eqnCommas)| "];\n";
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

importJulia = (importList,f) -> f<<concatenate apply(importList,pkg->pkg | "\n")


parseSolutions = out -> (
    L=select("\\[([0-9]|| |-|\\.|e|i|m|\\+|,)*\\]",out);
    sols=apply(L,s->point{apply(separate(",",replace("\\[|\\]","",replace("im","*ii",s))),x->value x)})
    )

-- todo: replace arbitrary variable name "f" for julia system
-- todo: remove bottom two lines from print
solveJulia = method(Options=>{})
solveJulia PolySystem := o -> P -> (
    if not isOpen JuliaProcess then error(noProcessError);
--    print "Wait for system to be solved:\n";
    writeSys P;
    JuliaProcess<<flush;
    x:=read JuliaProcess;
    JuliaProcess<<"sols=solve(f);\n";
--    JuliaProcess<<"x=[s.solution for s in sols];\n";
    JuliaProcess<<"show(IOContext(STDOUT, :compact=>false),[s.solution for s in sols])\n";
    finished = "done!!!x";
    JuliaProcess<<"print(\"" | finished |"\")\n";
    out="";
    while not match(finished,out) do (
	JuliaProcess<<flush;
    	out = out | read JuliaProcess;	
	);
    m:=regex("Paths tracked(.|\n)*x",out);
    n:=substring(first m, out);
    solutionString:=replace(finished, "",substring(last regex("Array.*",n),n));
    outputString:=replace("Array","",substring(first regex("Paths tracked:(.|\n)*Array",n),n));
    print outputString;
    parseSolutions solutionString
    )


    

-*
finite(sols)
solutions(sols)
results(sols)
*-

end

close JuliaProcess
restart
needs "julia.m2"
needs "./ExampleSystems/jointsR6.m2"
importJulia({"using HomotopyContinuation","import DynamicPolynomials: PolyVar"},JuliaProcess)

R=CC[x,y]
f= x^2+y
g=y^2-pi
P=polySystem {f,g}
out=solveJulia P;

#oo

out=replace("-","",replace("\n","xxx",solveJulia P));
match("Paths tracked.*Array",out)
substring(first lastMatch,out)

parseSolutions out

writeSys(P,"test.jl",WithImports=>true)

#sols

Q=polySystem jointsR6(CC_53)
sols=solveJulia Q;

#sols--too many?

writeSys(Q,"joints.jl",WithImports=>true)
