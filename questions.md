# Questions for Q&A session

- The ReactionNetworks Package uses indexed variables.  How do they work?  Consider the following code:

```
needsPackage("ReactionNetworks")
N=oneSiteModificationA()
R=createRing N
k=N.ReactionRates
```

  At this point, `k_0` should be recognized as a variable in the ring
  `R`.  It has the same name after all.  How do I get the variable in
  `R`?  How do I get the index `{0,1}` of the variable `k_{0,1}`.


- Your question here.
4) Suppose I wish to interact with another program using "openInOut" and anticipate that I will transfer more than 4096 bites. What are my options? Consider the following code:

```
n=10000
f=openInOut("!M2")
for j from 0 to 9 do (
    for i from 0 to (n-1) do (
    	f<<"1+1\n";
    	f<<flush;
    	);
    x=read f;
    print(#x);
    );
    
restart

n=10000
f=openInOut("!M2")
for j from 0 to 9 do (
    for i from 0 to (n-1) do (
    	f<<"1+1\n";
    	f<<flush;
	x=read f;
    	);
    print(#x);
    );
```

The loop above "restart" executes successfully on my machine, while the code below hangs (consistent with the warning given in the documentation for openInOut.) Is there a better general way communicate then reading immediately after writing?

- How does M2 relate to other computer algebra systems? The closest competitor seems to be singular? When should I choose one over the other? How do they differ in design philosophy?
