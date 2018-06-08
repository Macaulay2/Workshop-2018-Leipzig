--Khatri–Rao product
khatriRaoProduct = method()
khatriRaoProduct(Matrix,Matrix) := (V,W) -> (
    v:=numColumns V;
    w:=numColumns W;
    if v != w then(
	return "error: different columns size"
	);
    return
         matrix apply(for i from 0 to v-1 list V_i,
             for i from 0 to w-1 list W_i,
               (i,j)->i**j) )
----example 
khatriRaoProduct(M,N)
-----Check the Size Of Matrix
isSameSize=method()
isSameSize(Matrix, Matrix):= (V,W) -> ( 
    if numColumns V==numColumns W and numRows V==numRows W then
     return true
     else false
     )
----example 
isSameSize(M,N)
--------------------  
-- Hadmard Product
hadmardProduct = method()
hadmardProduct(Matrix,Matrix) := (V,W) -> (
    if not isSameSize(V,W) then(  return "error: different size"
	);
    matrix apply(entries M,entries N,(i,j)->( apply(i,j,times   )))
	   )    
----example 
hadmardProduct(M,N)     
