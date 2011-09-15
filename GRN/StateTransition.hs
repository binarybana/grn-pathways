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

import GRN.Parse
import GRN.Types
import GRN.Utils

import qualified Data.Vector as V
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as U

import Data.Vector.Strategies
import Control.Parallel.Strategies

import System.Cmd
import Control.Monad
import Data.List
import Data.Bits
import Data.Ord (comparing)
import Data.Maybe
import qualified Data.Map as M
import Text.Printf
import System.Environment
import System.Console.ParseArgs
import Data.Graph.Analysis 
import System.Random
import System.IO.Unsafe


fullIteration = {-# SCC "fullIteration" #-} iterateProbs2.iterateProbs 

-- |Push from edges back to the nodes
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

-- |Push to the edges from the nodes
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
convertProbsGG :: ColoredStateGraph -> ColoredStateGraph -> ColoredStateGraph
convertProbsGG oldg blankg = gmap copyProb (emap (\(EdgeInfo _ w)->(EdgeInfo 0.0 w)) blankg)
    where
        copyProb (adj1, node, NodeInfo a y, adj2) = (adj1, node, NodeInfo a p, adj2)
            where p = case (lab oldg node) of
                    Just (NodeInfo _ p) -> p 
                    Nothing -> 0.0

convertProbsVG :: SSD -> ColoredStateGraph -> ColoredStateGraph
convertProbsVG ssd blankg = gmap copyProb (emap (\(EdgeInfo _ w)->(EdgeInfo 0.0 w)) blankg)
    where
        copyProb (adj1, node, NodeInfo a y, adj2) = (adj1, node, NodeInfo a p, adj2)
            where p = ssd U.! node

convertProbsGV :: ColoredStateGraph -> SSD
convertProbsGV gr = U.fromList $ map (\(_,(NodeInfo _ p))->p) $ labNodes gr

buildKmaps :: ParseData -> KmapSet
buildKmaps = M.map buildKmap

buildKmap :: GeneInfo -> Kmap
buildKmap (GeneInfo n (Just val) _ _ _) = Kmap n [] (M.singleton [] (V 0 (valconv val)))
buildKmap (GeneInfo n _ _ [] []) = Kmap n [] (M.singleton [] Const)
buildKmap gi = Kmap (name gi) predictors (foldl updateMap initMap spways)
    where
        -- Sorted Pathways by ascending number of predictors
        spways = sortBy (comparing (\(Pathway x _ _ _)->length x)) $ pathways gi
        -- Sorted list of predictors
        predictors = sort.depends $ gi
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
            | head gene == '!' = filterPoses (not val) xs (tail gene)
            | otherwise = filter ((==(if val then 1 else 0)).(!!ind)) xs
                where Just ind = elemIndex gene predictors

-- Boolean to Double
valconv :: Bool -> Double
valconv val = if val then 1 else 0

fuzzGraphEdges :: ParseData -> KmapSet-> (Int, ColoredStateGraph) -> (Int, ColoredStateGraph)
fuzzGraphEdges pdata ks (x,orig) = (newX, mkGraph permedNodeStates edges)
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

kmapToStateGraph :: KmapSet -> ColoredStateGraph
kmapToStateGraph kset = mkGraph permedNodeStates edges
    where
        genes = M.keys kset
        nGenes = length genes
        nStates = 2^nGenes
        intStates = [0..nStates-1]
        permedStates = map (dec2bin nGenes) intStates
        stringStates = map (concatMap show) permedStates
        permedNodeStates =
            zip intStates (map (\x->NodeInfo x (1.0/(fromIntegral nStates))) stringStates)
        edges = concatMap genEdges permedStates
        genEdges :: [Int] -> [(Int,Int,EdgeInfo)]
        genEdges st = map (\(val,to)-> (from,bin2dec to,EdgeInfo 0.0 val)) tos
            where
                from = bin2dec st
                tos = [(1.0,[])] >>= foldr (>=>) return (stateKmapLus 0 st genes kset)
        
stateKmapLus :: Int -> [Int] -> [Gene] -> KmapSet -> [(Double,[Int]) -> [(Double,[Int])]]
stateKmapLus x st genes ks = map (evaluator.stateKmapLu x st genes) $ M.elems ks 

stateKmapLu :: Int -> [Int] -> [Gene] -> Kmap -> Double
stateKmapLu x st genes (Kmap g gs k) = case M.lookup loc k of
    Just X -> 0.5 --if x==0 then 0.5 else randUnsafe x
    Just (C _) -> 0.5 --if x==0 then 0.5 else randUnsafe x
    Just (V _ d) -> d
    Just Const -> fromIntegral $ st !! thisInd
    Nothing -> error ("Something wrong here " ++ show loc ++ show k ++ show inds ++ show gs ++ show genes)
    where
        loc = map (st!!) inds
        inds = concatMap (flip elemIndices genes) gs
        Just thisInd = elemIndex g genes 

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
            return y

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

-- Strips nodes that have no incoming edges
stripTransNodes :: Graph gr => gr NodeInfo b -> gr NodeInfo b
stripTransNodes orig = foldr strip orig (nodes orig)
    where   strip node gr
                | null $ inn orig node = delNode node gr
                | otherwise = gr
            prob i = let Just (NodeInfo _ p) = i in p 
   
-- Will only strip nodes that have no incoming edges and one outgoing edge
-- (deterministic transient nodes).
stripTransFuzzyNodes :: ColoredStateGraph -> ColoredStateGraph
stripTransFuzzyNodes orig = foldl' strip orig (nodes orig)
    where   
    strip gr node 
      | null inedges && null (tail outedges) = delNode node gr
      | otherwise = gr
      where 
      inedges = inn gr node
      outedges = out gr node

sumMass :: ColoredStateGraph -> Double
sumMass g = sum $ map ((\(NodeInfo _ pr)->pr).snd) (labNodes g)

calcSSAs :: Args String -> ParseData -> KmapSet -> ColoredStateGraph -> GeneMC
calcSSAs args p ks = zipNames . allgenes . ssaVec . simVec . seedList
  where
    n2 = getRequiredArg args "n2" -- Number of iterations per simulation
    n3 = getRequiredArg args "n3" -- Number of simruns
    n4 = getRequiredArg args "n4" -- Starting Seed
    pass n f = last.take n.iterate f 

    -- Create list of random seeds and the starting graph
    seedList :: ColoredStateGraph -> V.Vector (Int,ColoredStateGraph)
    seedList gr = V.fromList $ zip [n4..(n4+n3)] (repeat gr)

    -- Get a new randomized graph and simulate
    simVec :: V.Vector (Int,ColoredStateGraph) -> V.Vector (Int,ColoredStateGraph)
    simVec = G.map ((\(x,g)->(x,(pass n2 fullIteration g))).(fuzzGraphEdges p ks))  

    -- Now get the SSAs of these 
    ssaVec :: V.Vector (Int,ColoredStateGraph) -> V.Vector SSA
    ssaVec xv = G.map (genSSA.snd) xv `using` (parVector 2)

    -- Basically a transpose, to get a list of SSA for each gene instead of
    -- each simulation run
    allgenes :: V.Vector SSA -> V.Vector (U.Vector Double)
    allgenes ssas = G.map (\x-> G.convert $ G.map (G.!x) ssas) (V.fromList [0..(length $ M.keys p)-1])
    zipNames :: V.Vector (U.Vector Double) -> M.Map Gene (U.Vector Double)
    zipNames xv = M.fromList $ zip (M.keys p) (G.toList xv)

normalizeGraph :: ColoredStateGraph -> (Double, ColoredStateGraph)
normalizeGraph g = (tot, nmap (\(NodeInfo a pr)->(NodeInfo a (pr/tot))) g)
    where tot = sumMass g

printSSA :: SSA -> IO()
printSSA vec = do
    G.mapM_ (printf "%7.3f") vec
    putStrLn ""

-- Average the SSA over several graphs to avoid periodicities
genSSA :: ColoredStateGraph -> SSA
genSSA g = G.map (/denom) $ foldl1 (G.zipWith (+)) ssaList
    where
        denom = fromIntegral n
        n = 8
        ssaList = map genSSAForGraph glist
        glist = take n $ iterate fullIteration g

genSSAForGraph :: ColoredStateGraph -> SSA
genSSAForGraph g = G.foldr (G.zipWith (+)) (U.replicate len 0) filt2
    where   len = length $ name $ G.head filt1
            filt1 = V.fromList $ map snd (labNodes g)
            filt2 = G.map weight filt1 
            weight (NodeInfo a pr) = U.fromList $ map ((*pr).conv) a
            conv = (\x->if x=='1' then 1 else 0) 
            name (NodeInfo a pr) = a
