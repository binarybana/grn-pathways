{-# LANGUAGE BangPatterns #-}
-- |
-- Module    : GRN.StateTransition
-- Copyright : (c) 2011 Jason Knight
-- License   : BSD3
--
-- Maintainer  : jason@jasonknight.us
-- Stability   : experimental
-- Portability : portable
-- 
-- Handles the creation of a Karnaugh map or species pathway diagram from
-- ParseData. Also handles the simulation of the long run probabilities and
-- calcuation of the SSA Transform
--

module GRN.StateTransition where

import Data.Graph.Inductive
import GRN.Parse
import GRN.Types

import qualified Data.Vector as V
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as U

import System.Cmd
import Control.Monad
import Data.List
import Data.Ord (comparing)
import Data.Maybe
import qualified Data.Map as M
import Text.Printf
import System.Environment
import System.Console.ParseArgs
import Data.Graph.Analysis
import System.Random
import System.IO.Unsafe

dec2bin :: Int -> Int -> [Int]
dec2bin len y = (replicate (len-(length res)) 0) ++ res -- need to add padding
    where   res = reverse $ dec2bin' y
            dec2bin' 0 = []
            dec2bin' y = let (a,b) = quotRem y 2
                         in b : dec2bin' a
bin2dec :: [Int] -> Int
bin2dec st = foldl' (\a x->2*a+x) 0 st

fullIteration = {-# SCC "fullIteration" #-} iterateProbs2.iterateProbs 

iterateProbs2 :: ColoredStateGraph -> ColoredStateGraph
iterateProbs2 !sgraph = {-# SCC "iterate2" #-} gmap flowProbs sgraph
    where 
        flowProbs (adj1, node, !(NodeInfo a oldprob), adj2)=
            {-# SCC "flowProbs2" #-}(newIns, node, NodeInfo a newProb, newOuts)
            where   
                newProb = {-# SCC "newProb" #-} 
                    sum $ map (\(_,_,(EdgeInfo a _))->a) (inn sgraph node)
                newOuts = {-# SCC "newOuts" #-} map (\((EdgeInfo _ w),e)->((EdgeInfo 0.0 w),e)) adj2
                newIns = {-# SCC "newIns" #-} map (\((EdgeInfo _ w),e)->((EdgeInfo 0.0 w),e)) adj1

iterateProbs :: ColoredStateGraph -> ColoredStateGraph
iterateProbs !sgraph = {-# SCC "iterate1" #-} gmap flowProbs sgraph
    where   
        flowProbs  (adj1, node, !(NodeInfo a oldprob), adj2) =
            {-# SCC "flowProbs1" #-} (newIns, node, NodeInfo a 0, newOuts)  
            where   
                newOuts = 
                    {-# SCC "newOuts1" #-} map (\((EdgeInfo _ w),n) ->((EdgeInfo (oldprob*w) w),n)) adj2
                newIns = {-# SCC "newIns1" #-} map updateIns adj1
                updateIns ((EdgeInfo _ w),n) = {-# SCC "updateIns1" #-}
                    let (_,_,(NodeInfo _ p),_) = context sgraph n 
                        in ((EdgeInfo (p*w) w),n)

-- Take probabilities from oldg, and apply them to newg to allow for a change
-- in external conditions or something like that.
convertProbs :: ColoredStateGraph -> ColoredStateGraph -> ColoredStateGraph
convertProbs oldg blankg = gmap copyProb (emap (\(EdgeInfo _ w)->(EdgeInfo 0.0 w)) blankg)
    where
        copyProb (adj1, node, NodeInfo a y, adj2) = (adj1, node, NodeInfo a p, adj2)
            where p = case (lab oldg node) of
                    Just (NodeInfo _ p) -> p 
                    Nothing -> 0.0

buildKmaps :: ParseData -> KmapSet
buildKmaps = M.map buildKmap

buildKmap :: GeneInfo -> Kmap
buildKmap (GeneInfo n (Just val) _ _) = Kmap n [] (M.singleton [] (V 0 (valconv val)))
buildKmap (GeneInfo n _ [] []) = Kmap n [] (M.singleton [] Const)
buildKmap gi = Kmap (name gi) predictors (foldl updateMap initMap spways)
    where
        -- Sorted Pathways by ascending number of predictors
        spways = sortBy (comparing (\(Pathway x _ _ _)->length x)) $ pathways gi
        -- Sorted list of predictors
        predictors = sort.nub.concatMap (\(Pathway x _ _ _)->x) $ pathways gi
        nump = length predictors
        positions :: [[Int]]
        positions = [[0],[1]] >>= foldr (<=<) return (replicate (nump-1) (\x->[x++[0],x++[1]]))
        initMap = M.fromList $ zip positions (repeat X)
        updateMap :: M.Map [Int] Kentry -> Pathway -> M.Map [Int] Kentry 
        updateMap k (Pathway up pre post _) = foldl (\k pos-> updateLoc k pos (valconv post) (length up)) k upposses 
            where upposses = foldl (filterPoses pre) positions up
        updateLoc :: M.Map [Int] Kentry  -> [Int] -> Double -> Int -> M.Map [Int] Kentry 
        updateLoc k loc val num = M.insert loc (updateRule (k M.! loc) (V num val) num) k
        updateRule X val num = val
        updateRule (C i) val num 
            | i < num = val
            | i == num = C i
            | otherwise = error "Should never reach another pathway with less than the conflict"
        updateRule (V i d) val@(V ni nd) num
            | i < num = val
            | i == num = if d == nd then V i d else C i
            | otherwise = error "Should never reach this pathway 2"
        filterPoses :: Bool -> [[Int]] -> Gene -> [[Int]]
        filterPoses val xs gene
            | head gene == '!' = filterPoses (not val) xs gene
            | otherwise = filter ((==(if val then 1 else 0)).(!!ind)) xs
                where Just ind = elemIndex gene predictors

-- Boolean to Double
valconv :: Bool -> Double
valconv val = if val then 1 else 0

reDrawEdges :: ParseData -> KmapSet-> (Int, ColoredStateGraph) -> (Int, ColoredStateGraph)
reDrawEdges pdata ks (x,orig) = (newX, mkGraph permedNodeStates edges)
    where
        newX = if x>0 then x+1 else x -- 0 is the deterministic seed
        permedStates = map ((dec2bin nGenes).fst) permedNodeStates
        permedNodeStates = labNodes orig
        nGenes = length genes
        genes = M.keys pdata
        edges = concatMap genEdges permedStates 
        genEdges :: [Int] -> [(Int,Int,EdgeInfo)]
        genEdges st = map (\(val,to)-> (from,bin2dec to,EdgeInfo 0.0 val)) tos
            where   
                from = bin2dec st
                tos = [(1.0,[])] >>= foldr (>=>) return (stateKmapLus x st genes ks)


kmapToStateGraph :: ParseData -> KmapSet -> ColoredStateGraph
kmapToStateGraph pdata kset = mkGraph permedNodeStates edges 
    where
        permedStates = map (dec2bin nGenes) intStates 
        permedNodeStates = 
            zip intStates (map (\x->NodeInfo x (1.0/(fromIntegral nStates))) stringStates)
        stringStates = map (concatMap show) permedStates
        nStates = 2^nGenes 
        intStates = [0..nStates-1]
        nGenes = length genes
        genes = M.keys pdata
        edges = concatMap genEdges permedStates 
        genEdges :: [Int] -> [(Int,Int,EdgeInfo)]
        genEdges st = map (\(val,to)-> (from,bin2dec to,EdgeInfo 0.0 val)) tos
            where   
                from = bin2dec st
                tos = [(1.0,[])] >>= foldr (>=>) return (stateKmapLus 0 st genes kset)
        
stateKmapLus :: Int -> [Int] -> [Gene] -> KmapSet -> [(Double,[Int]) -> [(Double,[Int])]]
stateKmapLus x st genes ks = map (evaluator.stateKmapLu x st genes) $ M.elems ks 

evaluator :: Double -> ((Double,[Int]) -> [(Double,[Int])])
evaluator val 
    | val == 1.0 = det1
    | val == 0.0 = det0
    | otherwise = rand
    where rand (d,xs) = [(d*ival,xs++[0]),(d*val,xs++[1])]
          det0 (d,xs) = [(d,xs++[0])]
          det1 (d,xs) = [(d,xs++[1])]
          ival = 1.0 - val

randUnsafe x = unsafePerformIO $ do
            setStdGen (mkStdGen x)
            y <- randomRIO (0,1)
            if x==0 then return 0.5 else return y

stateKmapLu :: Int -> [Int] -> [Gene] -> Kmap -> Double
stateKmapLu x st genes (Kmap g gs k) = case M.lookup loc k of
    Just X -> randUnsafe x
    Just (C _) -> randUnsafe x
    Just (V _ d) -> d
    Just Const -> fromIntegral $ st !! thisInd
    Nothing -> error ("Something wrong here " ++ show loc ++ show k ++ show inds ++ show gs ++ show genes)
    where
        loc = map (st!!) inds
        inds = concatMap (flip elemIndices genes) gs
        Just thisInd = elemIndex g genes 

test = do
    con <- readFile "pws/nfkb-super.pw"
    let pdata = parsePW con
        ks = buildKmaps pdata
    return ks

buildBiGeneGraph :: ParseData -> DirectedGraph
buildBiGeneGraph pdata = 
    mkGraph (map (\(a,b)->(b,a)) name2NodeMap) (edges ++ redges)
    where   
        edges = concatMap splitGis gis 
            where   
                gis = (M.elems pdata)
                splitGis gi = concatMap (buildEdges (name gi)) (pathways gi)
        buildEdges down (Pathway ups b1 b2 _) = map buildEdge ups
                where buildEdge up = case (b1,b2) of
                        (True,True)   -> (n2N up,n2N down,"")
                        (True,False)  -> (n2N up,negn2N down,"")
                        (False,True)  -> (negn2N up,n2N down,"")
                        (False,False) -> (negn2N up,negn2N down,"")
        negn2N name = n2N ('n':name)
        filt name = if head name == '!' then tail name else name
        n2N name = case lookup (filt name) name2NodeMap of
            Just x -> x
            Nothing -> error "no node mapping for this gene"
        name2NodeMap = zip ((M.keys pdata) ++ ngenes) [1..]
        redges = zip3 (map n2N ngenes) (map (n2N.(\x->tail x))ngenes) 
            (repeat "red")
        ngenes = map ('n':) (M.keys pdata)

buildGeneGraph :: ParseData -> DirectedGraph
buildGeneGraph pdata = mkGraph lnodes edges
    where   nodes = zipWith (\a (b,c)->(a,b,c)) [1..] (M.assocs pdata)
            lnodes = map (\(a,b,c)->(a,b)) nodes
            edges = concatMap edges1 nodes 
            edges1 (n,g,gi) = 
                zip3 (map n2N (allUpStream gi)) (repeat n) (repeat "")
            filt name = if head name == '!' then tail name else name
            n2N name = case lookup (filt name) name2NodeMap of
                Just x -> x
                Nothing -> error (name ++ " gene: no node mapping for this gene")
            name2NodeMap = zip (M.keys pdata) [1..] 
            allUpStream = concatMap (\(Pathway xs _ _ _)->xs) . pathways

simulate :: Args String -> ColoredStateGraph -> ColoredStateGraph
simulate args 
    | gotArg args "reduce" = reduceSim
    | otherwise = regSim
      where 
            n1 = getRequiredArg args "n1"
            n2 = getRequiredArg args "n2"
            reduceSim = (pass n2 fullIteration).
                    (pass n1 (stripTransNodes.fullIteration)) 
            regSim = pass (n1+n2) fullIteration
            pass n f = last.take n.iterate f 

stripTransNodes :: Graph gr => gr NodeInfo b -> gr NodeInfo b
stripTransNodes orig = foldr strip orig (nodes orig)
    where   strip node gr
                | (null $ inn orig node) && prob (lab gr node) /= 0.0 = error ("huh? non zero probability? "++show node)
                | null $ inn orig node = delNode node gr
                | otherwise = gr
            prob i = let Just (NodeInfo _ p) = i in p 
    

genSSA :: ColoredStateGraph -> SSA
genSSA g = G.map (/denom) $ foldl1 (G.zipWith (+)) ssaList
    where
        denom = fromIntegral n
        n = 8
        ssaList = map genSSA' glist
        glist = take n $ iterate fullIteration g

sumMass :: ColoredStateGraph -> Double
sumMass g = sum $ map ((\(NodeInfo _ pr)->pr).snd) (labNodes g)

normalizeGraph :: ColoredStateGraph -> (Double, ColoredStateGraph)
normalizeGraph g = (tot, nmap (\(NodeInfo a pr)->(NodeInfo a (pr/tot))) g)
    where tot = sumMass g

genSSA' :: ColoredStateGraph -> SSA
genSSA' g = G.foldr (G.zipWith (+)) (U.replicate len 0) filt2
    where   len = length $ name $ G.head filt1
            filt1 = V.fromList $ map snd (labNodes g)
            filt2 = G.map weight filt1 
            weight (NodeInfo a pr) = U.fromList $ map ((*pr).conv) a
            conv = (\x->if x=='1' then 1 else 0) 
            name (NodeInfo a pr) = a

printSSA :: SSA -> IO()
printSSA vec = do
    G.mapM_ (printf "%7.3f") vec
    putStrLn ""

testFun = do
    con <- readFile "./pws/nfkb-super.pw"
    return $ testPure con
    
testPure con = finalGraph
    where
        pdata = parsePW $ con ++ (unlines.words $ "tnf=1 ltbr=0 lps=0") 
        ks = buildKmaps pdata
        sg = kmapToStateGraph pdata ks
        finalGraph = (pass 100 fullIteration).
                    (pass 20 (stripTransNodes.fullIteration)) $ sg
        pass n f = last.take n.iterate f 
