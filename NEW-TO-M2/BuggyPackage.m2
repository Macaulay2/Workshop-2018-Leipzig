-----------------------------------------------------------------
-- INSTRUCTIONS -------------------------------------------------
-- 1. Installing this package produces errors; debug it.
-- 2. Eliminate the warning messages
-- 3. Add a new test.
-- 4. Create a second copy of this package in a different directory 
--    and load each copy.  How do you do know which copy you have
--    installed?
-- 5. Bundle some existing Macaulay2 into your own package.
--
-- PREAMBLE -----------------------------------------------------
-- -*- coding: utf-8 -*-
newPackage(
    "BuggyPackage",
    Version => "0.1", 
    Date => "5 June 2018",
    Authors => {{Name => "Your Name", 
	    Email => "your.name@email.edu", 
	    HomePage => "http://your.website.html"}},
    Headline => "a buggy example of a Macaulay2 package",
    AuxiliaryFiles => false,
    DebuggingMode => true
    )

-- EXPORT LIST --------------------------------------------------
export {
    -- Types
    "Stack",
    -- methods
    "pop",
    "push",
    "toStack",
    }


-- CODE ---------------------------------------------------------

Stack = new Type of BasicList
Stack.synonym = "stack"
net Stack := S -> stack apply(S, x -> net x)

toStack = method(TypicalValue => Stack)
toStack BasicList := L -> new Stack from L

pop = method(TypicalValue => Stack)
pop Stack := S -> (
    if #S <= 0 then error "-- expected a nonempty stack"
    toStack drop(S,1)
    )

push = method(TypicalValue => Stack)
push (Thing,Stack) := (x,S) -> toStack prepend(x, S)

-- DOCUMENTATION ------------------------------------------------

beginDocumentation()
doc ///
  Key
    BuggyPackage
  Headline
     an example Macaulay2 package
  Description
   Text
    {\em BuggyPackage} is a small package creating a Stack datatype from BasicLists.
///    

doc ///
  Key
    (toStack,BasicList)
    toStack
  Headline
    creates a stack
  Usage
    toStack L
  Inputs
    L : BasicList 
  Outputs
    : Stack
      a stack containing element in the given order    
  Description
    Text
      A stack is an abstract data type that serves as a collection of elements,
      with two principal operations: push, which adds an element to the
      collection, and pop, which removes the most recently added element.
    Example
      S1 = toStack {1,2,3,4,5}
      S2 = pop S1
      S3 = push(1,S2)
      S1 === S3
    Text
      Here is a random example.
    Example
      S4 = toStack for i to 5 list random(5)
      S5 = pop S4
      S6 = pop S5
      
///

doc ///
  Key
    (toStack,BasicList)
    toStack
  Headline
    creates a stack
  Usage
    toStack L
  Inputs
    L : BasicList 
  Outputs
    : Stack
      a stack containing element in the given order    
  Description
    Text
      A stack is an abstract data type that serves as a collection of elements,
      with two principal operations: push, which adds an element to the
      collection, and pop, which removes the most recently added element.
    Example
      S1 = toStack {1,2,3,4,5}
      S2 = pop S1
      S3 = push(1,S2)
      S1 === S3
    Text
      Here is a random example.
    Example
      S4 = toStack for i to 5 list random(5)
      S5 = pop S4
      S6 = pop S5
///

doc ///
  Key
    Stack
    (net, Stack)
  Headline
    a linear data structure
  Description
    Text
      A stack is an abstract data type that serves as a collection of elements,
      with two principal operations: push, which adds an element to the
      collection, and pop, which removes the most recently added element.
    Example
      S1 = toStack {1,2,3,4,5}
      S2 = pop S1
      S3 = push(1,S2)
      S1 === S3
    Text
      Here is a random example.
    Example
      S4 = toStack for i to 5 list random(5)
      S5 = pop S4
      S6 = pop S5
///

doc ///
  Key
    pop
    (pop, Stack)
  Headline
    removes the most recently added element
  Usage
    pop S
  Inputs
    S : Stack
  Outputs
    : Stack
      a stack with one less element
  Description
    Text
      A stack is an abstract data type that serves as a collection of elements,
      with two principal operations: push, which adds an element to the
      collection, and pop, which removes the most recently added element that
      was not yet removed.
    Example
      S1 = toStack {1,2,3,4,5}
      S2 = pop S1
      S3 = push(1,S2)
      S1 === S3
    Text
      Here is a random example.
    Example
      S4 = toStack for i to 5 list random(5)
      S5 = pop S4
      S6 = pop S5
  SeeAlso
    push
///

-- TESTS --------------------------------------------------------

TEST ///
    S = toStack {}
    assert try (pop S; false) else true
    S1 = toStack {1,2,3,4,5}
    S2 = pop S1
    S3 = push(1,S2)
    assert(S1 === S3)
///

end--------------------------------------------------------------
-----------------------------------------------------------------

uninstallPackage "BuggyPackage"
restart
installPackage "BuggyPackage"
check "BuggyPackage"


path = prepend(path, "~/yourDirectory")


