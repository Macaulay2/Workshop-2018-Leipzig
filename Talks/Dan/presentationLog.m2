Macaulay2, version 1.11.1
with packages: ConwayPolynomials, Elimination, IntegralClosure, InverseSystems, LLLBases,
               PrimaryDecomposition, ReesAlgebra, TangentCone

i1 : viewHelp

i2 : viewHelp resolution 

i3 : help resolution 

o3 = resolution -- projective resolution
     ***********************************

     Synopsis
     ========

       * Optional inputs:
           * DegreeLimit => ...,  -- compute only up to this degree
           * FastNonminimal => ..., 
           * HardDegreeLimit => ..., 
           * LengthLimit => ...,  -- stop when the resolution reaches this length
           * PairLimit => ...,  -- stop when this number of pairs has been handled
           * SortStrategy => ..., 
           * StopBeforeComputation => ...,  -- whether to stop the computation immediately
           * Strategy => ..., 
           * SyzygyLimit => ...,  -- stop when this number of syzygies are obtained


     Ways to use resolution :
     ========================

       * "resolution(Ideal)" -- compute a projective resolution of (the quotient ring
         corresponding to) an ideal
       * resolution(MonomialIdeal), see "resolution(Ideal)" -- compute a projective
         resolution of (the quotient ring corresponding to) an ideal
       * "resolution(Matrix)" -- given a module map represented by a matrix, produce a
         comparison map between resolutions of its source and target
       * "resolution(Module)" -- compute a free resolution of a module

o3 : DIV

i4 : showStructure 

o4 = Thing : BasicList : Command
                         Constant
                         DocumentTag
                         Eliminate
                         Expression : Adjacent
                                      AssociativeExpression : Equation
                                                              Product
                                                              Sum
                                      BinaryOperation
                                      Divide
                                      FunctionApplication
                                      Holder : OneExpression
                                               ZeroExpression
                                      MatrixExpression
                                      Minus
                                      NonAssociativeProduct
                                      Parenthesize
                                      Power
                                      RowExpression
                                      SparseMonomialVectorExpression
                                      SparseVectorExpression
                                      Subscript
                                      Superscript
                                      Table
                         FilePosition
                         ForestNode
                         Hybrid
                         IndeterminateNumber
                         IndexedVariable
                         InfiniteNumber
                         LowerBound
                         Manipulator
                         MutableList : Bag
                         Option
                         Partition
                         ProductOrder
                         PushforwardComputation
                         RingElement
                         SumOfTwists
                         Time
                         TreeNode
                         URL
                         Vector
                         VisibleList : Array
                                       List : VerticalList : NumberedVerticalList
                                       Sequence
             Boolean
             CompiledFunctionBody
             Database
             Dictionary : GlobalDictionary
                          LocalDictionary
             File
             Function : CompiledFunction
                        CompiledFunctionClosure : MethodFunction
                        FunctionClosure : CacheFunction
                                          MethodFunctionWithOptions
             FunctionBody
             HashTable : CoherentSheaf
                         Ideal : MonomialIdeal
                         ImmutableType : Module
                         ModuleMap : Matrix
                         MonoidElement
                         MutableHashTable : CacheTable
                                            Descent
                                            GradedModule : ChainComplex
                                            GradedModuleMap : ChainComplexMap
                                            GroebnerBasis
                                            IndexedVariableTable
                                            Package
                                            Resolution
                                            ScriptedFunctor
                                            Type : HeaderType
                                                   Monoid : OrderedMonoid : GeneralOrderedMonoid
                                                   Ring : EngineRing : FractionField
                                                                       GaloisField
                                                                       InexactField : ComplexField
                                                                                      RealField
                                                                       LocalRing
                                                                       PolynomialRing
                                                                       QuotientRing
                                                   RingFamily : InexactFieldFamily
                                                   SelfInitializingType
                                                   WrapperType
                                            Variety : AffineVariety
                                                      ProjectiveVariety
                         MutableMatrix
                         OptionTable : GroebnerBasisOptions
                         ProjectiveHilbertPolynomial
                         RingMap
                         SheafOfRings
                         VirtualTally : BettiTally
                                        Tally : Set
             Net : String
             NetFile
             Nothing : InexactNumber  : CC
                                    *     *
                                        RR
                                          *
             Number : InexactNumber : CC
                                      RR
                      QQ
                      ZZ
             Pseudocode
             Symbol : Keyword
             SymbolBody
             Task

o4 : Descent

i5 : needsPackage "NormalToricVarieties"
--loading configuration for package "FourTiTwo" from file /Users/dan/Library/Application Support/Macaulay2/init-FourTiTwo.m2

o5 = NormalToricVarieties

o5 : Package

i6 : showStructure 

o6 = Thing : BasicList : Command
                         Constant
                         DocumentTag
                         Eliminate
                         Expression : Adjacent
                                      AssociativeExpression : Equation
                                                              Product
                                                              Sum
                                      BinaryOperation
                                      Divide
                                      FunctionApplication
                                      Holder : OneExpression
                                               ZeroExpression
                                      MatrixExpression
                                      Minus
                                      NonAssociativeProduct
                                      Parenthesize
                                      Power
                                      RowExpression
                                      SparseMonomialVectorExpression
                                      SparseVectorExpression
                                      Subscript
                                      Superscript
                                      Table
                         FilePosition
                         ForestNode
                         Hybrid
                         IndeterminateNumber
                         IndexedVariable
                         InfiniteNumber
                         LowerBound
                         Manipulator
                         MutableList : Bag
                         Option
                         Partition
                         ProductOrder
                         PushforwardComputation
                         RingElement
                         SumOfTwists
                         Time
                         TreeNode
                         URL
                         Vector
                         VisibleList : Array
                                       List : VerticalList : NumberedVerticalList
                                       Sequence
             Boolean
             CompiledFunctionBody
             Database
             Dictionary : GlobalDictionary
                          LocalDictionary
             File
             Function : CompiledFunction
                        CompiledFunctionClosure : MethodFunction
                        FunctionClosure : CacheFunction
                                          MethodFunctionWithOptions
             FunctionBody
             HashTable : CoherentSheaf
                         Ideal : MonomialIdeal
                         ImmutableType : Module
                         ModuleMap : Matrix
                         MonoidElement
                         MutableHashTable : CacheTable
                                            Descent
                                            GradedModule : ChainComplex
                                            GradedModuleMap : ChainComplexMap
                                            GroebnerBasis
                                            IndexedVariableTable
                                            Package
                                            PolyhedralObject : Cone
                                                               Fan
                                                               PolyhedralComplex
                                                               Polyhedron
                                            Resolution
                                            ScriptedFunctor
                                            Type : HeaderType
                                                   Monoid : OrderedMonoid : GeneralOrderedMonoid
                                                   Ring : EngineRing : FractionField
                                                                       GaloisField
                                                                       InexactField : ComplexField
                                                                                      RealField
                                                                       LocalRing
                                                                       PolynomialRing
                                                                       QuotientRing
                                                   RingFamily : InexactFieldFamily
                                                   SelfInitializingType
                                                   WrapperType
                                            Variety : AffineVariety
                                                      NormalToricVariety
                                                      ProjectiveVariety
                         MutableMatrix
                         OptionTable : GroebnerBasisOptions
                         ProjectiveHilbertPolynomial
                         RingMap
                         SheafOfRings
                         ToricDivisor
                         VirtualTally : BettiTally
                                        Tally : Set
             Net : String
             NetFile
             Nothing : InexactNumber  : CC
                                    *     *
                                        RR
                                          *
             Number : InexactNumber : CC
                                      RR
                      QQ
                      ZZ
             Pseudocode
             Symbol : Keyword
             SymbolBody
             Task

o6 : Descent

i7 : parent GlobalDictionary 

o7 = Dictionary

o7 : Type

i8 : class GlobalDictionary 

o8 = Type

o8 : Type

i9 : showClassStructure 

o9 = Type : AffineVariety
            Array
            AssociativeExpression
            BasicList
            BettiTally
            Boolean : false
                      true
            CacheFunction
            CacheTable
            CC
              *
            ChainComplex
            ChainComplexMap
            CoherentSheaf
            CompiledFunction : abs
                               acos
                               addCancelTask
                               addDependencyTask
                               addStartTask
                               agm
                               alarm
                               ancestor
                               any
                               append
                               apply
                               applyKeys
                               applyPairs
                               applyValues
                               ascii
                               asin
                               atan
                               atan2
                               atEndOfFile
                               BesselJ
                               BesselY
                               cancelTask
                               characters
                               class
                               clearEcho
                               collectGarbage
                               combine
                               commandInterpreter
                               concatenate
                               connectionCount
                               copy
                               cos
                               cosh
                               cot
                               coth
                               cpuTime
                               createTask
                               csc
                               csch
                               currentDirectory
                               currentLineNumber
                               currentTime
                               deepSplice
                               difference
                               disassemble
                               drop
                               dumpdata
                               echoOff
                               echoOn
                               eint
                               erase
                               erf
                               erfc
                               exec
                               expm1
                               fileExecutable
                               fileExists
                               fileLength
                               fileMode
                               fileReadable
                               fileTime
                               fileWritable
                               firstkey
                               flagLookup
                               fork
                               format
                               frames
                               functionBody
                               Gamma
                               GCstats
                               get
                               getc
                               getenv
                               getGlobalSymbol
                               getNetFile
                               groupID
                               hash
                               hashTable
                               horizontalJoin
                               identity
                               imaginaryPart
                               installMethod
                               instance
                               isANumber
                               isCanceled
                               isDirectory
                               isFinite
                               isGlobalSymbol
                               isInfinite
                               isInputFile
                               isListener
                               isOpen
                               isOutputFile
                               isReady
                               isRegularFile
                               join
                               keys
                               kill
                               limitFiles
                               limitProcesses
                               linkFile
                               loaddata
                               localDictionaries
                               locate
                               log
                               log1p
                               lookup
                               lookupCount
                               merge
                               mergePairs
                               mingle
                               minimizeFilename
                               minus
                               mkdir
                               mutable
                               newClass
                               newNetFile
                               nextkey
                               openDatabase
                               openDatabaseOut
                               openFiles
                               openIn
                               openInOut
                               openListener
                               openOut
                               openOutAppend
                               override
                               pack
                               pairs
                               parent
                               plus
                               power
                               powermod
                               prepend
                               printString
                               processID
                               protect
                               pseudocode
                               read
                               readDirectory
                               readlink
                               realPart
                               realpath
                               recursionDepth
                               regex
                               registerFinalizer
                               relativizeFilename
                               remove
                               removeDirectory
                               removeFile
                               reorganize
                               reverse
                               run
                               scan
                               scanPairs
                               schedule
                               sec
                               sech
                               select
                               separate
                               sequence
                               serialNumber
                               set
                               setEcho
                               setGroupID
                               setIOExclusive
                               setIOSynchronized
                               setIOUnSynchronized
                               sin
                               sinh
                               size2
                               sleep
                               splice
                               sqrt
                               stack
                               substring
                               symbolBody
                               symlinkFile
                               take
                               tally
                               tan
                               tanh
                               taskResult
                               times
                               toCC
                               toList
                               toRR
                               toSequence
                               uncurry
                               unsequence
                               unstack
                               utf8
                               utf8check
                               values
                               wait
                               wrap
                               xor
                               youngest
                               zeta
            CompiledFunctionBody
            CompiledFunctionClosure : code
                                      commonest
                                      directSum
                                      EXAMPLE
                                      examples
                                      export
                                      exportMutable
                                      expression
                                      flatten
                                      gcd
                                      gradedModule
                                      gradedModuleMap
                                      hold
                                      html
                                      hypertext
                                      ideal
                                      info
                                      intersect
                                      isSorted
                                      lcm
                                      length
                                      makePackageIndex
                                      mathML
                                      max
                                      maxPosition
                                      methods
                                      min
                                      minPosition
                                      monomialIdeal
                                      net
                                      options
                                      package
                                      pretty
                                      runLengthEncode
                                      tex
                                      texMath
                                      toExternalString
                                      toString
                                      transpose
                                      undocumented
                                      unique
                                      vars
            ComplexField
            Cone
            Constant : EulerConstant
                       ii
                       pi
            Database
            Descent
            Dictionary
            DocumentTag
            EngineRing
            Expression
            Fan
            File : stderr
                   stdio
            FilePosition
            ForestNode
            FractionField
            Function
            FunctionBody
            FunctionClosure : addEndFunction
                              addStartFunction
                              ancestors
                              applicationDirectory
                              applicationDirectorySuffix
                              applyTable
                              assert
                              baseFilename
                              beginDocumentation
                              benchmark
                              cacheValue
                              centerString
                              chainComplex
                              columnate
                              delete
                              demark
                              End
                              error
                              even
                              first
                              getNonUnit
                              getSymbol
                              globalAssign
                              globalAssignFunction
                              globalAssignment
                              globalReleaseFunction
                              infoHelp
                              input
                              installedPackages
                              integrate
                              inversePermutation
                              isFinitePrimeField
                              isPrimitive
                              isTable
                              last
                              lines
                              load
                              makeDocumentTag
                              method
                              mod
                              monoid
                              monomialCurveIdeal
                              needs
                              notImplemented
                              number
                              odd
                              on
                              pager
                              peek
                              print
                              same
                              seeParsing
                              showHtml
                              stashValue
                              subtable
                              synonym
                              SYNOPSIS
                              syzygyScheme
                              table
                              temporaryFileName
                              testPackage
                              toAbsolutePath
                              toLower
                              toUpper
                              tutorial
                              uniform
                              uninstallAllPackages
                              userSymbols
                              zero
            GaloisField
            GeneralOrderedMonoid
            GlobalDictionary : OutputDictionary
                               PackageDictionary
            GradedModule
            GradedModuleMap
            GroebnerBasis
            GroebnerBasisOptions
            HashTable
            HeaderType : Adjacent
                         BinaryOperation
                         Divide
                         Equation
                         FunctionApplication
                         MatrixExpression
                         Power
                         RowExpression
                         SparseMonomialVectorExpression
                         SparseVectorExpression
                         Subscript
                         Superscript
                         Table
            Ideal
            ImmutableType
            IndeterminateNumber : indeterminate
            IndexedVariable
            IndexedVariableTable
            InexactField
            InexactFieldFamily : CC
                                 RR
            InexactNumber
            InexactNumber
                         *
            InfiniteNumber : infinity
            Keyword
            List
            LocalDictionary
            LocalRing
            Manipulator : close
                          closeIn
                          closeOut
                          endl
                          flush
            Matrix
            MethodFunction : accumulate
                             acosh
                             acot
                             addCone
                             addHook
                             addPolyhedron
                             adjoint
                             adjoint'
                             affineHull
                             affineImage
                             affinePreimage
                             all
                             ambDim
                             ambient
                             antipode
                             apropos
                             areCompatible
                             asinh
                             autoload
                             baseName
                             between
                             binomial
                             bipyramid
                             blowup
                             borel
                             capture
                             cartierDivisorGroup
                             ccRefinement
                             ceiling
                             cellDecompose
                             char
                             chi
                             classGroup
                             clean
                             coefficient
                             coefficientRing
                             coimage
                             cokernel
                             columnAdd
                             columnMult
                             columnPermute
                             columnRankProfile
                             columnSwap
                             commonFace
                             commonRing
                             comodule
                             complement
                             complete
                             components
                             compose
                             compositions
                             compress
                             conductor
                             cone
                             coneFromHData
                             coneFromVData
                             cones
                             conjugate
                             contains
                             content
                             contract
                             contract'
                             convexHull
                             conwayPolynomial
                             cover
                             coverMap
                             crossPolytope
                             cyclicPolytope
                             debug
                             decompose
                             default
                             degree
                             degreeLength
                             degrees
                             degreesMonoid
                             degreesRing
                             denominator
                             depth
                             describe
                             diagonalMatrix
                             dictionary
                             diff
                             diff'
                             dim
                             directProduct
                             discriminant
                             dismiss
                             divideByVariable
                             dropAlternatives
                             dualCone
                             dualFaceRepresentationMap
                             eagonNorthcott
                             ehrhart
                             elements
                             eliminate
                             emptyPolyhedron
                             endPackage
                             entries
                             euler
                             eulers
                             exp
                             expectedReesIdeal
                             exponents
                             exportFrom
                             faceFan
                             faces
                             facesAsCones
                             facesAsPolyhedra
                             facets
                             fan
                             fanFromGfan
                             Fano
                             findSynonyms
                             fittingIdeal
                             flip
                             floor
                             fold
                             frac
                             fraction
                             fromCDivToPic
                             fromCDivToWDiv
                             fromDividedPowers
                             fromPicToCl
                             fromWDivToCl
                             fVector
                             gbRemove
                             gbSnapshot
                             gcdCoefficients
                             genera
                             generateAssertions
                             generator
                             genericMatrix
                             genericSkewMatrix
                             genericSymmetricMatrix
                             genus
                             getChangeMatrix
                             getMatrix
                             getWWW
                             gramm
                             halfspaces
                             heft
                             height
                             hilbertFunction
                             hirzebruch
                             Hom
                             homogenize
                             homomorphism
                             homomorphism'
                             httpHeaders
                             hypercube
                             hyperplanes
                             icFractions
                             icMap
                             icPIdeal
                             image
                             imageFan
                             incompCones
                             incompPolyhedra
                             index
                             indices
                             inducesWellDefinedMap
                             inInterior
                             insert
                             installAssignmentMethod
                             installHilbertFunction
                             instances
                             interiorLatticePoints
                             interiorPoint
                             interiorVector
                             intersection
                             inverse
                             irreducibleCharacteristicSeries
                             irreducibleDecomposition
                             isAffineRing
                             isAmple
                             isBorel
                             isCartier
                             isCommutative
                             isCompact
                             isComplete
                             isConstant
                             isDegenerate
                             isDirectSum
                             isEffective
                             isEmpty
                             isFace
                             isFano
                             isField
                             isFreeModule
                             isFullDimensional
                             isHomogeneous
                             isIdeal
                             isInjective
                             isIsomorphism
                             isLatticePolytope
                             isModule
                             isMonomialIdeal
                             isNef
                             isNormal
                             isPointed
                             isPolynomialRing
                             isPolytopal
                             isPrimary
                             isPrime
                             isProjective
                             isPseudoprime
                             isPure
                             isQQCartier
                             isQuotientModule
                             isQuotientOf
                             isQuotientRing
                             isReal
                             isReflexive
                             isRing
                             isSimplicial
                             isSkewCommutative
                             isSmooth
                             isSquareFree
                             isStandardGradedPolynomialRing
                             isSubmodule
                             isSubquotient
                             isSubset
                             isSurjective
                             isUnit
                             isVeryAmple
                             isWellDefined
                             isWeylAlgebra
                             jacobian
                             koszul
                             latticePoints
                             latticeVolume
                             leadCoefficient
                             leadComponent
                             leadMonomial
                             leadTerm
                             liftable
                             linealitySpace
                             linearTransform
                             linSpace
                             listForm
                             listSymbols
                             lngamma
                             loadAlternative
                             localRing
                             LUdecomposition
                             makeDirectory
                             match
                             maxCones
                             maxFace
                             maxPolyhedra
                             member
                             memoize
                             methodOptions
                             minFace
                             minimalPrimes
                             minkowskiSum
                             minkSummandCone
                             mixedVolume
                             module
                             monomialSubideal
                             multidegree
                             nefGenerators
                             newCoordinateSystem
                             newtonPolytope
                             nextPrime
                             norm
                             normalFan
                             nullhomotopy
                             nullSpace
                             numColumns
                             numerator
                             numeric
                             numgens
                             numRows
                             nVertices
                             objectiveVector
                             ofClass
                             orbits
                             pad
                             part
                             partition
                             partitions
                             parts
                             pdim
                             peek'
                             permanents
                             permutations
                             pfaffians
                             picardGroup
                             pivots
                             poincare
                             poincareN
                             polar
                             polarFace
                             poly
                             polyhedra
                             polyhedralComplex
                             polyhedron
                             polyhedronFromHData
                             polytope
                             positions
                             posOrthant
                             precision
                             preimage
                             presentation
                             product
                             profile
                             Proj
                             projectiveHilbertPolynomial
                             promote
                             proximum
                             pseudoRemainder
                             putMatrix
                             pyramid
                             QRDecomposition
                             quotient'
                             quotientRemainder
                             quotientRemainder'
                             randomKRationalPoint
                             rank
                             rays
                             reduceHilbert
                             reductionNumber
                             relations
                             remainder
                             remainder'
                             removeHook
                             removeLowestDimension
                             replace
                             reshape
                             resultant
                             ring
                             rotate
                             round
                             rowAdd
                             rowMult
                             rowPermute
                             rowRankProfile
                             rowSwap
                             runHooks
                             saveSession
                             scanKeys
                             scanLines
                             scanValues
                             schreyerOrder
                             searchPath
                             secondaryPolytope
                             selectInSubring
                             selectVariables
                             separateRegexp
                             setRandomSeed
                             setup
                             setupEmacs
                             sheaf
                             sheafHom
                             show
                             singularLocus
                             size
                             skeleton
                             smallestFace
                             smoothSubfan
                             someTerms
                             source
                             Spec
                             splitWWW
                             standardForm
                             standardPairs
                             statePolytope
                             stdSimplex
                             stellarSubdivision
                             sublatticeBasis
                             sublists
                             submatrix
                             submatrix'
                             submatrixByDegrees
                             subquotient
                             subsets
                             substitute
                             sum
                             super
                             support
                             switch
                             sylvesterMatrix
                             symmetricPower
                             tailCone
                             target
                             tensorAssociativity
                             terms
                             TEST
                             toBinomial
                             toDividedPowers
                             toField
                             topCoefficients
                             topComponents
                             toricCircuits
                             toricGraver
                             toricGraverDegrees
                             toSublattice
                             trace
                             triangulate
                             truncate
                             truncateOutput
                             ultimate
                             unbag
                             use
                             value
                             variety
                             vector
                             versalEmbedding
                             vertexEdgeMatrix
                             vertexFacetMatrix
                             vertices
                             volume
                             wedgeProduct
                             weightRange
                             weilDivisorGroup
                             whichGm
                             width
                             Wikipedia
            MethodFunctionWithOptions : about
                                        affineSpace
                                        analyticSpread
                                        annihilator
                                        associatedPrimes
                                        basis
                                        betti
                                        check
                                        codim
                                        coefficients
                                        cohomology
                                        copyDirectory
                                        copyFile
                                        cotangentSheaf
                                        determinant
                                        distinguished
                                        document
                                        dual
                                        eigenvalues
                                        eigenvectors
                                        extend
                                        exteriorPower
                                        factor
                                        fillMatrix
                                        findFiles
                                        findHeft
                                        flattenRing
                                        forceGB
                                        fromDual
                                        gb
                                        gcdLLL
                                        generators
                                        getPrimeWithRootOfUnity
                                        GF
                                        graphIdeal
                                        graphRing
                                        Grassmannian
                                        groebnerBasis
                                        hermite
                                        hilbertBasis
                                        hilbertPolynomial
                                        hilbertSeries
                                        hirzebruchSurface
                                        homology
                                        icFracP
                                        idealizer
                                        independentSets
                                        inducedMap
                                        installPackage
                                        integralClosure
                                        integralClosures
                                        intersectInP
                                        inverseSystem
                                        isLinearType
                                        isLLL
                                        isReduction
                                        jacobianDual
                                        kernel
                                        kernelLLL
                                        kleinschmidt
                                        lift
                                        LLL
                                        loadPackage
                                        localize
                                        makeS2
                                        makeSimplicial
                                        makeSmooth
                                        map
                                        markedGB
                                        matrix
                                        mingens
                                        minimalBetti
                                        minimalPresentation
                                        minimalReduction
                                        minors
                                        modulo
                                        monomials
                                        moveFile
                                        multiplicity
                                        mutableIdentity
                                        mutableMatrix
                                        needsPackage
                                        netList
                                        newPackage
                                        newRing
                                        normalCone
                                        normalToricVariety
                                        position
                                        primaryComponent
                                        primaryDecomposition
                                        projectiveSpace
                                        prune
                                        pushForward
                                        quotient
                                        radical
                                        random
                                        randomMutableMatrix
                                        reesAlgebra
                                        reesIdeal
                                        regularity
                                        resolution
                                        ringFromFractions
                                        roots
                                        rsort
                                        saturate
                                        Schubert
                                        showTex
                                        smallAmpleToricDivisor
                                        smithNormalForm
                                        smoothFanoToricVariety
                                        solve
                                        sort
                                        sortColumns
                                        specialFiber
                                        specialFiberIdeal
                                        status
                                        SVD
                                        symlinkDirectory
                                        symmetricAlgebra
                                        symmetricAlgebraIdeal
                                        symmetricKernel
                                        syz
                                        tangentCone
                                        tangentSheaf
                                        tensor
                                        toDual
                                        toricDivisor
                                        toricGroebner
                                        toricMarkov
                                        trim
                                        uninstallPackage
                                        weightedProjectiveSpace
            Module
            ModuleMap
            Monoid
            MonoidElement
            MonomialIdeal
            MutableHashTable
            MutableList
            MutableMatrix
            Net
            NetFile
            NormalToricVariety
            Nothing
            Number
            OneExpression
            Option
            OptionTable
            OrderedMonoid
            Package : Classic
                      ConwayPolynomials
                      Core
                      Elimination
                      FourierMotzkin
                      FourTiTwo
                      IntegralClosure
                      InverseSystems
                      LLLBases
                      Macaulay2Doc
                      Normaliz
                      NormalToricVarieties
                      Parsing
                      Polyhedra
                      PrimaryDecomposition
                      ReesAlgebra
                      TangentCone
                      Text
                      User
            Partition
            PolyhedralComplex
            PolyhedralObject
            Polyhedron
            PolynomialRing
            ProjectiveHilbertPolynomial
            ProjectiveVariety
            Pseudocode
            QuotientRing
            RealField
            Resolution
            Ring : QQ
                   ZZ
            RingElement
            RingFamily
            RingMap
            RR
              *
            ScriptedFunctor : Ext
                              HH
                              hh
                              id
                              OO
                              sheafExt
                              Tor
            SelfInitializingType : Bag
                                   Command : clearAll
                                             clearOutput
                                             edit
                                             exit
                                             help
                                             listLocalSymbols
                                             listUserSymbols
                                             profileSummary
                                             quit
                                             restart
                                             showClassStructure
                                             showStructure
                                             showUserStructure
                                             viewHelp
                                   Eliminate
                                   Hybrid
                                   LowerBound
                                   NumberedVerticalList
                                   ProductOrder
                                   PushforwardComputation
                                   URL
                                   VerticalList
            Sequence
            Set
            SheafOfRings
            String
            SumOfTwists
            Symbol
            SymbolBody
            Tally
            Task
            Thing
            Time
            ToricDivisor
            TreeNode
            Variety
            Vector
            VirtualTally
            VisibleList
            WrapperType : Holder
                          Minus
                          NonAssociativeProduct
                          Parenthesize
                          Product
                          Sum
            ZeroExpression

o9 : Descent

i10 : NormalToricVarieties 

o10 = NormalToricVarieties

o10 : Package

i11 : 3.4/1

o11 = 3.4

o11 : RR (of precision 53)

i12 : ancestors List

o12 = {List, VisibleList, BasicList, Thing}

o12 : List

i13 : {2,3,4}

o13 = {2, 3, 4}

o13 : List

i14 : methods List

o14 = {(%, List, Number)                                   }
      {(%, List, RingElement)                              }
      {(*, Thing, List)                                    }
      {(++, OptionTable, List)                             }
      {(+, List, List)                                     }
      {(-*Function*-, List)                                }
      {(-*Function*-, List)                                }
      {(-*Function*-, List)                                }
      {(-*Function*-, List)                                }
      {(-*Function*-, List, Function)                      }
      {(-*Function*-, Matrix, Matrix, List)                }
      {(-*Function*-, Ring, List)                          }
      {(-, List)                                           }
      {(-, List, List)                                     }
      {(-, List, Set)                                      }
      {(-, Set, List)                                      }
      {(.., List, List)                                    }
      {(..<, List, List)                                   }
      {(/, List, Command)                                  }
      {(/, List, Function)                                 }
      {(/, List, RingMap)                                  }
      {(/, List, SelfInitializingType)                     }
      {(/, List, Thing)                                    }
      {(/, Module, List)                                   }
      {(/, Ring, List)                                     }
      {(//, List, Number)                                  }
      {(//, List, RingElement)                             }
      {(<<, List, Thing)                                   }
      {(>>, List, Function)                                }
      {(?, List, List)                                     }
      {(\, RingMap, List)                                  }
      {(^, Matrix, List)                                   }
      {(^, Module, List)                                   }
      {(^, MutableMatrix, List)                            }
      {(^, Ring, List)                                     }
      {(^, SheafOfRings, List)                             }
      {(_, Ideal, List)                                    }
      {(_, Matrix, List)                                   }
      {(_, Module, List)                                   }
      {(_, Monoid, List)                                   }
      {(_, MutableMatrix, List)                            }
      {(_, PolynomialRing, List)                           }
      {(_, Ring, List)                                     }
      {(_, VisibleList, List)                              }
      {(|, List, List)                                     }
      {(addCone, List, Fan)                                }
      {(addPolyhedron, List, PolyhedralComplex)            }
      {(ascii, List)                                       }
      {(basis, InfiniteNumber, List, Ideal)                }
      {(basis, InfiniteNumber, List, Matrix)               }
      {(basis, InfiniteNumber, List, Module)               }
      {(basis, InfiniteNumber, List, Ring)                 }
      {(basis, List, Ideal)                                }
      {(basis, List, InfiniteNumber, Ideal)                }
      {(basis, List, InfiniteNumber, Matrix)               }
      {(basis, List, InfiniteNumber, Module)               }
      {(basis, List, InfiniteNumber, Ring)                 }
      {(basis, List, List, Ideal)                          }
      {(basis, List, List, Matrix)                         }
      {(basis, List, List, Module)                         }
      {(basis, List, List, Ring)                           }
      {(basis, List, Matrix)                               }
      {(basis, List, Module)                               }
      {(basis, List, Ring)                                 }
      {(basis, List, ZZ, Ideal)                            }
      {(basis, List, ZZ, Ring)                             }
      {(basis, ZZ, List, Ideal)                            }
      {(basis, ZZ, List, Ring)                             }
      {(blowup, List, NormalToricVariety)                  }
      {(blowup, List, NormalToricVariety, List)            }
      {(chainComplex, List)                                }
      {(checkConesForMaximality, List)                     }
      {(code, List)                                        }
      {(columnPermute, MutableMatrix, ZZ, List)            }
      {(commonFace, List)                                  }
      {(commonRing, List)                                  }
      {(coneFromVData, List)                               }
      {(contains, List, Cone)                              }
      {(contains, List, Polyhedron)                        }
      {(convexHull, List)                                  }
      {(createMonomialSubalgebra, List)                    }
      {(degreesMonoid, List)                               }
      {(degreesRing, List)                                 }
      {(diagonalMatrix, List)                              }
      {(diagonalMatrix, Ring, List)                        }
      {(diagonalMatrix, Ring, ZZ, ZZ, List)                }
      {(diagonalMatrix, RingFamily, List)                  }
      {(diagonalMatrix, RingFamily, ZZ, ZZ, List)          }
      {(diagonalMatrix, ZZ, ZZ, List)                      }
      {(directSum, List)                                   }
      {(document, List)                                    }
      {(doWriteNmzData, List)                              }
      {(drop, BasicList, List)                             }
      {(dual, MonomialIdeal, List)                         }
      {(ehrhartRing, List)                                 }
      {(ehrhartRing, List, RingElement)                    }
      {(eliminate, Ideal, List)                            }
      {(eliminate, List, Ideal)                            }
      {(export, List)                                      }
      {(exportFrom, Package, List)                         }
      {(exportMutable, List)                               }
      {(facesOfCone, List, ZZ)                             }
      {(facesOfCone, Matrix, List)                         }
      {(fan, List)                                         }
      {(fan, Matrix, List)                                 }
      {(fan, Matrix, Matrix, List)                         }
      {(fanFromGfan, List)                                 }
      {(findFiles, List)                                   }
      {(findHeft, List)                                    }
      {(fourierMotzkinElimination, List, List, ZZ)         }
      {(gcd, List)                                         }
      {(gcdLLL, List)                                      }
      {(generateAssertions, List)                          }
      {(getSaneOutputPart, List)                           }
      {(gradedModule, List)                                }
      {(gradedModuleMap, List)                             }
      {(hashTable, List)                                   }
      {(hasProperties, PolyhedralObject, List)             }
      {(help, List)                                        }
      {(hilbertFunction, List, CoherentSheaf)              }
      {(hilbertFunction, List, Ideal)                      }
      {(hilbertFunction, List, Module)                     }
      {(hilbertFunction, List, ProjectiveVariety)          }
      {(hilbertFunction, List, Ring)                       }
      {(homogenize, Matrix, RingElement, List)             }
      {(homogenize, Module, RingElement, List)             }
      {(homogenize, RingElement, RingElement, List)        }
      {(homogenize, Vector, RingElement, List)             }
      {(hypertext, List)                                   }
      {(ideal, List)                                       }
      {(incompCones, List)                                 }
      {(incompPolyhedra, List)                             }
      {(intclToricRing, List)                              }
      {(intersect, List)                                   }
      {(intersection, List)                                }
      {(isRedundant, List, Set)                            }
      {(kleinschmidt, ZZ, List)                            }
      {(lcm, List)                                         }
      {(listSymbols, List)                                 }
      {(makePackageIndex, List)                            }
      {(makePrimitive, List)                               }
      {(map, Module, Module, List)                         }
      {(map, Module, Module, RingMap, List)                }
      {(map, Module, Nothing, List)                        }
      {(map, Module, Nothing, RingMap, List)               }
      {(map, Module, ZZ, List)                             }
      {(map, Ring, Ring, List)                             }
      {(mathML, List)                                      }
      {(matrix, List)                                      }
      {(matrix, Ring, List)                                }
      {(matrix, RingFamily, List)                          }
      {(matrixFromVectorList, List, ZZ, Ring)              }
      {(memoize, Function, List)                           }
      {(mixedVolume, List)                                 }
      {(monoid, List)                                      }
      {(monomialIdeal, List)                               }
      {(mutableMatrix, List)                               }
      {(net, List)                                         }
      {(NewFromMethod, DocumentTag, List)                  }
      {(NewFromMethod, HashTable, List)                    }
      {(NewFromMethod, Module, List)                       }
      {(NewFromMethod, Set, List)                          }
      {(norm, List)                                        }
      {(normaliz, List)                                    }
      {(normalToricRing, List)                             }
      {(normalToricVariety, List, List)                    }
      {(ofClass, List)                                     }
      {(part, List, RingElement)                           }
      {(peek', ZZ, List)                                   }
      {(polyhedralComplex, List)                           }
      {(polyhedralComplex, Matrix, List)                   }
      {(polyhedralComplex, Matrix, Matrix, Matrix, List)   }
      {(precedence, List)                                  }
      {(preMatrixFromStringList, List, List)               }
      {(pretty2, List)                                     }
      {(primitive, List)                                   }
      {(primitive, List)                                   }
      {(product, List)                                     }
      {(promote, List, CC , CC )                           }
      {                  *    *                            }
      {(promote, List, QQ, CC )                            }
      {                      *                             }
      {(promote, List, QQ, QQ)                             }
      {(promote, List, QQ, RR )                            }
      {                      *                             }
      {(promote, List, RR , CC )                           }
      {                  *    *                            }
      {(promote, List, RR , RR )                           }
      {                  *    *                            }
      {(promote, List, ZZ, CC )                            }
      {                      *                             }
      {(promote, List, ZZ, QQ)                             }
      {(promote, List, ZZ, RR )                            }
      {                      *                             }
      {(promote, List, ZZ, ZZ)                             }
      {(random, List)                                      }
      {(random, List, Ring)                                }
      {(regularSubdivision, NormalToricVariety, List, List)}
      {(rowPermute, MutableMatrix, ZZ, List)               }
      {(rsort, List)                                       }
      {(runFMAlternativeOnInput, String, List)             }
      {(runNormaliz, List)                                 }
      {(scanLines, Function, List)                         }
      {(searchPath, List, String)                          }
      {(selectVariables, List, PolynomialRing)             }
      {(slice, Vector, List)                               }
      {(sort, List)                                        }
      {(SPACE, HeaderType, List)                           }
      {(SPACE, Ring, List)                                 }
      {(SPACE, WrapperType, List)                          }
      {(standardPairs, MonomialIdeal, List)                }
      {(submatrixByDegrees, Matrix, List, List)            }
      {(subsets, List)                                     }
      {(subsets, List, ZZ)                                 }
      {(substitute, Ideal, List)                           }
      {(substitute, Matrix, List)                          }
      {(substitute, Module, List)                          }
      {(substitute, RingElement, List)                     }
      {(substitute, Vector, List)                          }
      {(sum, List)                                         }
      {(SYNOPSIS, List)                                    }
      {(take, BasicList, List)                             }
      {(TEST, List)                                        }
      {(texMath, List)                                     }
      {(toricDivisor, List, NormalToricVariety)            }
      {(toZZ, List)                                        }
      {(toZZ, List)                                        }
      {(transpose, List)                                   }
      {(truncate, List, Ideal)                             }
      {(truncate, List, Module)                            }
      {(undocumented, List)                                }
      {(unique, List)                                      }
      {(vars, List)                                        }
      {(vector, List)                                      }
      {(weightedProjectiveSpace, List)                     }
      {(weightRange, List, RingElement)                    }
      {(writeNmzData, List)                                }

o14 : VerticalList

i15 : {2,3}+{3,4}

o15 = {5, 7}

o15 : List

i16 : x = new BasicList from {2,3}

o16 = BasicList{2, 3}

o16 : BasicList

i17 : x+x
stdio:17:2:(3): error: no method for binary operator + applied to objects:
--            BasicList{2, 3} (of class BasicList)
--      +     BasicList{2, 3} (of class BasicList)

i18 : R = QQ[x]
stdio:18:5:(3): error: encountered object not usable as variable at position 0 in list:
        BasicList{2, 3} (of class BasicList)

i19 : R = QQ[a]

o19 = R

o19 : PolynomialRing

i20 : a

o20 = a

o20 : R

i21 : class a

o21 = R

o21 : PolynomialRing

i22 : parent oo

o22 = RingElement

o22 : Type

i23 : parent oo

o23 = BasicList

o23 : Type

i24 : peek a

o24 = R{a}

i25 : a#0

o25 = a

o25 : RawRingElement

i26 : NormalToricVariety

o26 = NormalToricVariety

o26 : Type

i27 : NormalToricVarieties 

o27 = NormalToricVarieties

o27 : Package

i28 : peek oo

o28 = Package{auxiliary files => /Users/dan/src/M2/M2-Macaulay2/M2/Macaulay2/packages/NormalToricVarieties/                                                                                                                                                                                                                                                                                                                                                                                                                                                 }
              configuration file name => /Users/dan/Library/Application Support/Macaulay2/init-NormalToricVarieties.m2
              Dictionary => NormalToricVarieties.Dictionary
              documentation not loaded => true
              example data files => MutableHashTable{}
              example inputs => MutableHashTable{}
              example results => MutableHashTable{}
              exported mutable symbols => {}
              exported symbols => {makeSimplicial, cartierDivisorGroup, kleinschmidt, isEffective, normalToricVariety, fromPicToCl, affineSpace, projectiveSpace, classGroup, orbits, isProjective, WeilToClass, isCartier, toricDivisor, nefGenerators, emsBound, isAmple, picardGroup, fromCDivToPic, makeSmooth, NormalToricVariety, hirzebruchSurface, isDegenerate, smoothFanoToricVariety, smallAmpleToricDivisor, isNef, isFano, rawHHOO, fromWDivToCl, isQQCartier, blowup, fromCDivToWDiv, ToricDivisor, weightedProjectiveSpace, weilDivisorGroup}
              index.html => /Users/dan/src/M2/M2-Macaulay2/M2/BUILD/dan/builds.tmp/einsteinium-release-1.11.1/usr-dist/common/share/doc/Macaulay2/NormalToricVarieties/html/index.html
              loadDepth => 3
              Options => OptionTable{Authors => {{Name => Gregory G. Smith, Email => ggsmith@mast.queensu.ca, HomePage => http://www.mast.queensu.ca/~ggsmith}}}
                                     AuxiliaryFiles => true
                                     CacheExampleOutput => null
                                     Certification => null
                                     Configuration => OptionTable{}
                                     Date => 31 May 2017
                                     DebuggingMode => false
                                     Headline => normal toric varieties
                                     HomePage => null
                                     InfoDirSection => Macaulay2 and its packages
                                     OptionalComponentsPresent => true
                                     PackageExports => {Polyhedra}
                                     PackageImports => {FourierMotzkin, Normaliz}
                                     Reload => false
                                     UseCachedExampleOutput => false
                                     Version => 1.5
              package prefix => /Users/dan/src/M2/M2-Macaulay2/M2/BUILD/dan/builds.tmp/einsteinium-release-1.11.1/usr-dist/
              private dictionary => NormalToricVarieties#"private dictionary"
              processed documentation => MutableHashTable{}
              raw documentation => MutableHashTable{}
              raw documentation database => /Users/dan/src/M2/M2-Macaulay2/M2/BUILD/dan/builds.tmp/einsteinium-release-1.11.1/usr-dist/x86_64-Darwin-MacOS-10.13.4/lib/Macaulay2/x86_64-Darwin-MacOS-10.13.4/NormalToricVarieties/cache/rawdocumentation-dcba-8.db
              source directory => /Users/dan/src/M2/M2-Macaulay2/M2/Macaulay2/packages/
              source file => /Users/dan/src/M2/M2-Macaulay2/M2/Macaulay2/packages/NormalToricVarieties.m2
              test inputs => MutableHashTable{}
              test number => 0
              title => NormalToricVarieties
              undocumented keys => MutableHashTable{}

i29 : NormalToricVarieties#"source file"

o29 = /Users/dan/src/M2/M2-Macaulay2/M2/Macaulay2/packages/NormalToricVarieties.m2

i30 : x = new HashTable from { "a" => 3 , NormalToricVarieties => "hi" }

o30 = HashTable{a => 3                    }
                NormalToricVarieties => hi

o30 : HashTable

i31 : x#a
stdio:31:2:(3): error: key not found in hash table

i32 : a

o32 = a

o32 : R

i33 : x#"a"

o33 = 3

i34 : x#NormalToricVarieties

o34 = hi

i35 : showStructure 

o35 = Thing : BasicList : Command
                          Constant
                          DocumentTag
                          Eliminate
                          Expression : Adjacent
                                       AssociativeExpression : Equation
                                                               Product
                                                               Sum
                                       BinaryOperation
                                       Divide
                                       FunctionApplication
                                       Holder : OneExpression
                                                ZeroExpression
                                       MatrixExpression
                                       Minus
                                       NonAssociativeProduct
                                       Parenthesize
                                       Power
                                       RowExpression
                                       SparseMonomialVectorExpression
                                       SparseVectorExpression
                                       Subscript
                                       Superscript
                                       Table
                          FilePosition
                          ForestNode
                          Hybrid
                          IndeterminateNumber
                          IndexedVariable
                          InfiniteNumber
                          LowerBound
                          Manipulator
                          MutableList : Bag
                          Option
                          Partition
                          ProductOrder
                          PushforwardComputation
                          RingElement : R
                          SumOfTwists
                          Time
                          TreeNode
                          URL
                          Vector
                          VisibleList : Array
                                        List : VerticalList : NumberedVerticalList
                                        Sequence
              Boolean
              CompiledFunctionBody
              Database
              Dictionary : GlobalDictionary
                           LocalDictionary
              File
              Function : CompiledFunction
                         CompiledFunctionClosure : MethodFunction
                         FunctionClosure : CacheFunction
                                           MethodFunctionWithOptions
              FunctionBody
              HashTable : CoherentSheaf
                          Ideal : MonomialIdeal
                          ImmutableType : Module
                          ModuleMap : Matrix
                          MonoidElement
                          MutableHashTable : CacheTable
                                             Descent
                                             GradedModule : ChainComplex
                                             GradedModuleMap : ChainComplexMap
                                             GroebnerBasis
                                             IndexedVariableTable
                                             Package
                                             PolyhedralObject : Cone
                                                                Fan
                                                                PolyhedralComplex
                                                                Polyhedron
                                             Resolution
                                             ScriptedFunctor
                                             Type : HeaderType
                                                    Monoid : OrderedMonoid : GeneralOrderedMonoid
                                                    Ring : EngineRing : FractionField
                                                                        GaloisField
                                                                        InexactField : ComplexField
                                                                                       RealField
                                                                        LocalRing
                                                                        PolynomialRing
                                                                        QuotientRing
                                                    RingFamily : InexactFieldFamily
                                                    SelfInitializingType
                                                    WrapperType
                                             Variety : AffineVariety
                                                       NormalToricVariety
                                                       ProjectiveVariety
                          MutableMatrix
                          OptionTable : GroebnerBasisOptions
                          ProjectiveHilbertPolynomial
                          RingMap
                          SheafOfRings
                          ToricDivisor
                          VirtualTally : BettiTally
                                         Tally : Set
              Net : String
              NetFile
              Nothing : InexactNumber  : CC
                                     *     *
                                         RR
                                           *
              Number : InexactNumber : CC
                                       RR
                       QQ
                       ZZ
              Pseudocode
              Symbol : Keyword
              SymbolBody
              Task

o35 : Descent

i36 : "asdfasd asdf"

o36 = asdfasd asdf

i37 : oo || oo

o37 = asdfasd asdf
      asdfasd asdf

i38 : oo^1

      asdfasd asdf
o38 = asdfasd asdf

i39 : height o37

o39 = 1

i40 : height o38

o40 = 2

i41 : depth o37

o41 = 1

i42 : depth o38

o42 = 0

i43 : width o38

o43 = 12

i44 : x

o44 = HashTable{a => 3                    }
                NormalToricVarieties => hi

o44 : HashTable

i45 : x=symbol x

o45 = x

o45 : Symbol

i46 : QQ[x_1 .. x_4]

o46 = QQ[x , x , x , x ]
          1   2   3   4

o46 : PolynomialRing

i47 : (x_1 + x_2)^3

       3     2         2    3
o47 = x  + 3x x  + 3x x  + x
       1     1 2     1 2    2

o47 : QQ[x , x , x , x ]
          1   2   3   4

i48 : net o47

       3     2         2    3
o48 = x  + 3x x  + 3x x  + x
       1     1 2     1 2    2

i49 : height oo

o49 = 2

i50 : depth ooo

o50 = 1

i51 : o37

o51 = asdfasd asdf
      asdfasd asdf

i52 : o37|o37

o52 = asdfasd asdfasdfasd asdf
      asdfasd asdfasdfasd asdf

i53 : o37||o37

o53 = asdfasd asdf
      asdfasd asdf
      asdfasd asdf
      asdfasd asdf

i54 : stack (o37,o37,o37,o37)

o54 = asdfasd asdf
      asdfasd asdf
      asdfasd asdf
      asdfasd asdf
      asdfasd asdf
      asdfasd asdf
      asdfasd asdf
      asdfasd asdf

i55 : horizontalJoin  (o37,o37,o37,o37)

o55 = asdfasd asdfasdfasd asdfasdfasd asdfasdfasd asdf
      asdfasd asdfasdfasd asdfasdfasd asdfasdfasd asdf

i56 : x

o56 = x

o56 : IndexedVariableTable

i57 : x_1

o57 = x
       1

o57 : QQ[x , x , x , x ]
          1   2   3   4

i58 : I = ideal (x_1,x_2)

o58 = ideal (x , x )
              1   2

o58 : Ideal of QQ[x , x , x , x ]
                   1   2   3   4

i59 : J = ideal (x_1,x_2+x_1)

o59 = ideal (x , x  + x )
              1   1    2

o59 : Ideal of QQ[x , x , x , x ]
                   1   2   3   4

i60 : I == J

o60 = true

i61 : I === J

o61 = false

i62 : peek I

o62 = Ideal{cache => CacheTable{...1...}}
            generators => | x_1 x_2 |
            ring => QQ[x , x , x , x ]
                        1   2   3   4

i63 : peek J

o63 = Ideal{cache => CacheTable{...1...} }
            generators => | x_1 x_1+x_2 |
            ring => QQ[x , x , x , x ]
                        1   2   3   4

i64 : S = set {I,J}

o64 = set {ideal (x , x  + x ), ideal (x , x )}
                   1   1    2           1   2

o64 : Set

i65 : S#I

o65 = 1

i66 : peek S

o66 = Set{ideal (x , x  + x ) => 1}
                  1   1    2
          ideal (x , x ) => 1
                  1   2

i67 : tally {3,4,5,5,5}

o67 = Tally{3 => 1}
            4 => 1
            5 => 3

o67 : Tally

i68 : S == S
stdio:68:3:(3): error: no method for binary operator == applied to objects:
--            set {ideal (x , x  + x ), ideal (x , x )} (of class Set)
--                         1   1    2           1   2
--     ==     set {ideal (x , x  + x ), ideal (x , x )} (of class Set)
--                         1   1    2           1   2

i69 : S

o69 = set {ideal (x , x  + x ), ideal (x , x )}
                   1   1    2           1   2

o69 : Set

i70 : {S,S,S}

o70 = {set {ideal (x , x  + x ), ideal (x , x )}, set {ideal (x , x  + x ), ideal (x , x )},
                    1   1    2           1   2                 1   1    2           1   2   
      --------------------------------------------------------------------------------------
      set {ideal (x , x  + x ), ideal (x , x )}}
                   1   1    2           1   2

o70 : List

i71 : symbol S

o71 = S

o71 : Symbol

i72 : value oo

o72 = set {ideal (x , x  + x ), ideal (x , x )}
                   1   1    2           1   2

o72 : Set

i73 : S=4

o73 = 4

i74 : o72

o74 = set {ideal (x , x  + x ), ideal (x , x )}
                   1   1    2           1   2

o74 : Set

i75 : mutable S

o75 = false

i76 : R

o76 = R

o76 : PolynomialRing

i77 : class R

o77 = PolynomialRing

o77 : Type

i78 : ancestors oo

o78 = {PolynomialRing, EngineRing, Ring, Type, MutableHashTable, HashTable, Thing}

o78 : List

i79 : R#"hi" = "there"

o79 = there

i80 : R#"hi"

o80 = there

i81 : R

o81 = R

o81 : PolynomialRing

i82 : x_1

o82 = x
       1

o82 : QQ[x , x , x , x ]
          1   2   3   4

i83 : x_1 ** x_2

o83 = | x_1x_2 |

                                 1                          1
o83 : Matrix (QQ[x , x , x , x ])  <--- (QQ[x , x , x , x ])
                  1   2   3   4              1   2   3   4

i84 : x_1 ++ x_2

o84 = | x_1 0   |
      | 0   x_2 |

                                 2                          2
o84 : Matrix (QQ[x , x , x , x ])  <--- (QQ[x , x , x , x ])
                  1   2   3   4              1   2   3   4

i85 : x_1 \ x_2
stdio:85:5:(3): error: no method for binary operator \ applied to objects:
--            x  (of class QQ[x , x , x , x ])
--             1               1   2   3   4
--      \     x  (of class QQ[x , x , x , x ])
--             2               1   2   3   4

i86 : R\R := (r,s) -> "hi there"

o86 = -*Function[stdio:86:14-86:14]*-

o86 : FunctionClosure

i87 : x_1 \ x_2
stdio:87:5:(3): error: no method for binary operator \ applied to objects:
--            x  (of class QQ[x , x , x , x ])
--             1               1   2   3   4
--      \     x  (of class QQ[x , x , x , x ])
--             2               1   2   3   4

i88 : gens R

o88 = {a}

o88 : List

i89 : a\a

o89 = hi there

i90 : methods symbol \

o90 = {((\, =), Thing, Thing)                }
      {((\, =), Type, Type)                  }
      {(\, Command, Set)                     }
      {(\, Command, VirtualTally)            }
      {(\, Command, VisibleList)             }
      {(\, Function, Ideal)                  }
      {(\, Function, Set)                    }
      {(\, Function, VirtualTally)           }
      {(\, Function, VisibleList)            }
      {(\, R, R)                             }
      {(\, RingMap, List)                    }
      {(\, SelfInitializingType, VisibleList)}
      {(\, Thing, Thing)                     }

o90 : VerticalList

i91 : peek R

o91 = PolynomialRing of RingElement{(?, R, R) => -*Function[../../../../Macaulay2/m2/enginering.m2:170:21-170:36]*-                                                                                                     }
                                    (\, R, R) => -*Function[stdio:86:14-86:14]*-
                                    (_, R, monoid[a, Degrees => {1}, Heft => {1}, MonomialOrder => {MonomialSize => 32}, DegreeRank => 1]) => -*Function[../../../../Macaulay2/m2/orderedmonoidrings.m2:240:27-240:78]*-
                                                                                                   {GRevLex => {1}    }
                                                                                                   {Position => Up    }
                                    (lift, List, R, QQ) => -*Function[../../../../Macaulay2/m2/enginering.m2:195:44-195:65]*-
                                    (lift, List, R, R) => -*Function[../../../../Macaulay2/m2/enginering.m2:159:29-159:43]*-
                                    (lift, List, R, ZZ) => -*Function[../../../../Macaulay2/m2/enginering.m2:216:50-216:75]*-
                                    (lift, Matrix, R, QQ) => -*Function[../../../../Macaulay2/m2/enginering.m2:191:46-191:87]*-
                                    (lift, Matrix, R, R) => -*Function[../../../../Macaulay2/m2/enginering.m2:160:31-160:45]*-
                                    (lift, Matrix, R, ZZ) => -*Function[../../../../Macaulay2/m2/enginering.m2:212:49-212:91]*-
                                    (lift, Module, R, QQ) => -*Function[../../../../Macaulay2/m2/enginering.m2:189:46-189:80]*-
                                    (lift, Module, R, R) => -*Function[../../../../Macaulay2/m2/enginering.m2:162:31-162:45]*-
                                    (lift, Module, R, ZZ) => -*Function[../../../../Macaulay2/m2/enginering.m2:210:49-210:84]*-
                                    (lift, MutableMatrix, R, QQ) => -*Function[../../../../Macaulay2/m2/enginering.m2:193:53-193:99]*-
                                    (lift, MutableMatrix, R, R) => -*Function[../../../../Macaulay2/m2/enginering.m2:163:38-163:52]*-
                                    (lift, MutableMatrix, R, ZZ) => -*Function[../../../../Macaulay2/m2/enginering.m2:214:56-214:105]*-
                                    (lift, R, QQ) => -*Function[../../../../Macaulay2/m2/enginering.m2:187:39-187:75]*-
                                    (lift, R, R) => -*Function[../../../../Macaulay2/m2/enginering.m2:155:24-155:36]*-
                                    (lift, R, ZZ) => -*Function[../../../../Macaulay2/m2/enginering.m2:208:39-208:74]*-
                                    (promote, List, QQ, R) => -*Function[../../../../Macaulay2/m2/enginering.m2:194:50-194:62]*-
                                    (promote, List, R, R) => -*Function[../../../../Macaulay2/m2/enginering.m2:157:35-157:38]*-
                                    (promote, List, ZZ, R) => -*Function[../../../../Macaulay2/m2/enginering.m2:215:50-215:78]*-
                                    (promote, Matrix, QQ, R) => -*Function[../../../../Macaulay2/m2/enginering.m2:190:52-190:78]*-
                                    (promote, Matrix, R, R) => -*Function[../../../../Macaulay2/m2/enginering.m2:158:37-158:40]*-
                                    (promote, Matrix, ZZ, R) => -*Function[../../../../Macaulay2/m2/enginering.m2:211:52-211:79]*-
                                    (promote, Module, QQ, R) => -*Function[../../../../Macaulay2/m2/enginering.m2:188:52-188:78]*-
                                    (promote, Module, R, R) => -*Function[../../../../Macaulay2/m2/enginering.m2:161:37-161:40]*-
                                    (promote, Module, ZZ, R) => -*Function[../../../../Macaulay2/m2/enginering.m2:209:52-209:79]*-
                                    (promote, MutableMatrix, QQ, R) => -*Function[../../../../Macaulay2/m2/enginering.m2:192:59-192:90]*-
                                    (promote, MutableMatrix, R, R) => -*Function[../../../../Macaulay2/m2/enginering.m2:164:44-164:47]*-
                                    (promote, MutableMatrix, ZZ, R) => -*Function[../../../../Macaulay2/m2/enginering.m2:213:59-213:93]*-
                                    (promote, QQ, R) => -*Function[../../../../Macaulay2/m2/enginering.m2:185:37-185:71]*-
                                    (promote, R, R) => -*Function[../../../../Macaulay2/m2/enginering.m2:156:28-156:31]*-
                                    (promote, ZZ, R) => -*Function[../../../../Macaulay2/m2/enginering.m2:206:37-206:71]*-
                                    0 => 0
                                    1 => 1
                                    basering => QQ
                                    baseRings => {ZZ, QQ}
                                    char => 0
                                    degreesMonoid => monoid[T, Degrees => {1}, MonomialOrder => {MonomialSize => 32}, DegreeRank => 1, Inverses => true, Global => false]
                                                                                                {Weights => {-1}   }
                                                                                                {GroupLex => 1     }
                                                                                                {Position => Up    }
                                    degreesRing => ZZ[T]
                                    Engine => true
                                    expression => -*Function[../../../../Macaulay2/m2/orderedmonoidrings.m2:246:35-253:43]*-
                                    factor => -*Function[../../../../Macaulay2/m2/orderedmonoidrings.m2:267:29-297:76]*-
                                    FlatMonoid => monoid[a, Degrees => {1}, Heft => {1}, MonomialOrder => {MonomialSize => 32}, DegreeRank => 1]
                                                                                                          {GRevLex => {1}    }
                                                                                                          {Position => Up    }
                                    generatorExpressions => {a}
                                    generators => {a}
                                    generatorSymbols => {a}
                                    hi => there
                                    indexStrings => HashTable{a => a}
                                    indexSymbols => HashTable{a => a}
                                    isCommutative => true
                                    isPrime => -*Function[../../../../Macaulay2/m2/orderedmonoidrings.m2:298:27-302:22]*-
                                    liftDegree => -*Function[../../../../Macaulay2/m2/enginering.m2:59:20-61:38]*-
                                    monoid => monoid[a, Degrees => {1}, Heft => {1}, MonomialOrder => {MonomialSize => 32}, DegreeRank => 1]
                                                                                                      {GRevLex => {1}    }
                                                                                                      {Position => Up    }
                                    numallvars => 1
                                    promoteDegree => -*Function[../../../../Macaulay2/m2/enginering.m2:59:20-61:38]*-
                                    raw creation log => -*a bagged function application expression*-
                                    RawRing => QQGMP[a,
                                                 DegreeLength => 1,
                                                 Degrees => {1},
                                                 Heft => {1},
                                                 MonomialOrder => {
                                                   GRevLex => {1},
                                                   Position => Up
                                                   }
                                                 ]

i92 : {a,b,c}

o92 = {a, b, c}

o92 : List

i93 : X = new Type of BasicList

o93 = X

o93 : Type

i94 : new X from {a,b,c}

o94 = X{a, b, c}

o94 : X

i95 : parent X

o95 = BasicList

o95 : Type

i96 : class X

o96 = Type

o96 : Type

i97 : o94 # 1

o97 = b

o97 : Symbol

i98 : o94 _ 1
stdio:98:5:(3): error: no method for binary operator _ applied to objects:
--            X{a, b, c} (of class X)
--      _     1 (of class ZZ)

i99 : # o94

o99 = 3

i100 : x = new X from {a,b,c}

o100 = X{a, b, c}

o100 : X

i101 : x

o101 = X{a, b, c}

o101 : X

i102 : net X := t -> "hi there"

o102 = -*Function[stdio:102:12-102:12]*-

o102 : FunctionClosure

i103 : x

o103 = hi there

o103 : X

i104 : new X from {a,b,c}

o104 = hi there

o104 : X

i105 : net X := t -> apply(t,net)

o105 = -*Function[stdio:105:12-105:23]*-

o105 : FunctionClosure

i106 : x

o106 = stdio:106:1:(3): error: at print: expected argument 2 to be a string, or list or sequence of strings and integers, or null

-- comment: what happened just above is that our function above did not return a net, so it confused the printing routines

o106 : X

i115 : net X := t -> horizontalJoin apply(toList t,net)

o115 = -*Function[stdio:115:12-115:45]*-

o115 : FunctionClosure

i116 : x

o116 = abc

o116 : X

i117 : net X := t -> horizontalJoin between_"," apply(toList t,net)

o117 = -*Function[stdio:117:12-117:57]*-

o117 : FunctionClosure

i118 : x

o118 = a,b,c

o118 : X

i119 : T = QQ[r,s,t]

o119 = T

o119 : PolynomialRing

i120 : C = res coker vars T

        1      3      3      1
o120 = T  <-- T  <-- T  <-- T  <-- 0
                                    
       0      1      2      3      4

o120 : ChainComplex

i121 : C.dd

            1                 3
o121 = 0 : T  <------------- T  : 1
                 | r s t |

            3                        3
       1 : T  <-------------------- T  : 2
                 {1} | -s -t 0  |
                 {1} | r  0  -t |
                 {1} | 0  r  s  |

            3                  1
       2 : T  <-------------- T  : 3
                 {2} | t  |
                 {2} | -s |
                 {2} | r  |

            1
       3 : T  <----- 0 : 4
                 0

o121 : ChainComplexMap

i122 : 

Process M2 finished
