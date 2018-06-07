restart 

loadPackage "EdgeIdeals"


num = 10
RING = QQ[x_1..x_num]

G = graph(RING, {{x_1,x_2},{x_2,x_3},{x_1,x_3},{x_3,x_4},{x_4,x_6},
	{x_3,x_5},{x_5,x_6},{x_6,x_7},{x_7,x_8},{x_8,x_9},{x_6,x_9},{x_10,x_8},{x_9,x_10}
	});

I=edgeIdeal G


monVars = (mon) -> ( 
    R := ring mon;
    X := gens R;
    flatten for i to length X -1 list (
       	while mon % X_i ==0 list 
	i+1
	do 
	mon=mon//X_i
	)
    )


toricEdgeIdeal = (GRAPH) -> (
    I:=edgeIdeal(GRAPH);
    genI:= flatten entries gens I;
    indicesE:= for i to length genI -1 list(monVars(genI_i));
    edgeVars := for i to length indicesE -1 list (e_(indicesE_i));
    R:=QQ[edgeVars];
    toricMap := map(RING,R,genI);
    return trim  ker toricMap;
    )

J=toricEdgeIdeal(G)

f = (a) -> print a