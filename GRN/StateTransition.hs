{-# LANGUAGE BangPatterns #-}

module GRN.StateTransition where

import Data.Graph.Inductive
import Data.Graph.Inductive.Graph
import qualified GRN.Graphviz as Gv
import GRN.Parse

import System.Cmd
import Control.Monad
import Data.List
import Data.Maybe
import qualified Data.Map as M
import Text.Printf
import System.Environment
import Data.Graph.Analysis

data NodeInfo = NodeInfo String !Double
type ColoredStateTGraph = Gr NodeInfo Double

instance Show NodeInfo where
    show (NodeInfo name num) = name ++ "\t" ++ (show num) ++ "\t"

type StateTGraph = Gr String String

dec2bin :: Int -> Int -> [Int]
dec2bin y len = (replicate (len-(length res)) 0) ++ res -- need to add padding
    where   res = reverse $ dec2bin' y
            dec2bin' 0 = []
            dec2bin' y = let (a,b) = quotRem y 2
                         in b : dec2bin' a
bin2dec :: [Int] -> Int
bin2dec st = foldl' (\a x->2*a+x) 0 st

fullIteration = {-# SCC "fullIteration" #-} iterateProbs2.iterateProbs 

iterateProbs2 :: ColoredStateTGraph -> ColoredStateTGraph
iterateProbs2 !sgraph = {-# SCC "iterate2" #-} gmap flowProbs sgraph
    where 
        flowProbs (adj1,node, !(NodeInfo a oldprob),adj2)=
            {-# SCC "flowProbs2" #-}(newIns, node, NodeInfo a newProb, newOuts)
            where   
                newProb = {-# SCC "newProb" #-} 
                    sum $ map (\(_,_,a)->a) (inn sgraph node)
                newOuts = {-# SCC "newOuts" #-} map (\(_,w)->(0,w)) adj2
                newIns = {-# SCC "newIns" #-} map (\(_,w)->(0,w)) adj1

iterateProbs :: ColoredStateTGraph -> ColoredStateTGraph
iterateProbs !sgraph = {-# SCC "iterate1" #-} gmap flowProbs sgraph
    where   
        flowProbs  (adj1, node, !(NodeInfo a oldprob), adj2) =
            {-# SCC "flowProbs1" #-} (newIns, node, NodeInfo a 0, newOuts)  
            where   
                newOuts = 
                    {-# SCC "newOuts1" #-} map (\(_,w)->(newSplitProb,w)) adj2
                newIns = {-# SCC "newIns1" #-} map updateIns adj1
                updateIns (l,n) = {-# SCC "updateIns1" #-}
                    let c@(_,_,(NodeInfo _ p),_) = context sgraph n 
                        num = fromIntegral $ toInteger $ outdeg' c         
                        in (p/num,n)
                newSplitProb = oldprob/outdegree
                outdegree = let c = context sgraph node
                    in fromIntegral $ toInteger $ outdeg' c

initProbs :: StateTGraph -> ColoredStateTGraph
initProbs sgraph = gmap initProb (emap (\x->0) sgraph) 
    where   
        initProb (adj1, node, a, adj2) = 
            (adj1, node, NodeInfo a (1.0/num), adj2)
        num = fromIntegral $ toInteger $ length $ nodes sgraph

buildStateTGraph :: ParseData -> StateTGraph
buildStateTGraph pdata = mkGraph permedNodeStates edges 
    where   permedStates = map (flip dec2bin nGenes) intStates
            permedNodeStates = 
                zip intStates (map (concatMap show) permedStates)
            intStates = [0..(2^nGenes)-1]
            nGenes = length (M.keys pdata)
            indexNames = zip [1..] (M.keys pdata)

            edges = concatMap genEdges permedStates 
            genEdges :: [Int] -> [(Int,Int,String)]
            genEdges st = map (\to-> (from,bin2dec to,"")) tos
                where   
                    from = bin2dec st
                    tos = [[]] >>= foldr (<=<) return (build $ length st)
                    build :: Int -> [([Int]->[[Int]])]
                    build 0 = [return.id]
                    build num = thisgene : build (num-1)
                        where thisgene = calcGene st num pdata

calcGene :: [Int] -> Int -> ParseData -> ([Int] -> [[Int]])
-- Given a state, and the gene number of interest, return a function
calcGene st prog pdata
    | hasKnockout = thisKnockout 
    | length thisDeps == 0 = thisCurr
    | otherwise = transInfo
    where   names = M.keys pdata
            hasKnockout = isJust $ knockout thisGI
            thisKnockout = case (fromMaybe False (knockout thisGI)) of
                False -> det0
                True  -> det1
            thisDeps = depends thisGI
            thisGI = getGI thisName
            thisPways = pathways thisGI
            thisName = names !! (prog-1)
            getInd n = let (Just ind) = elemIndex n names in ind
            getGI name = let (Just gi) = M.lookup name pdata in gi
            curr name = case (st !! (getInd name)) of
                0 -> det0
                1 -> det1
            thisCurr = curr thisName 

            transInfo 
                | null list     = rand
                | all (not) list = det0
                | all (id) list      = det1 
                | otherwise     = rand
                where   list = map snd . filter ((==maxlen).fst) $ tuplist
                        maxlen = maximum . map fst $ tuplist
                        tuplist = foldr searchPWs [] thisPways
                        searchPWs pw@(Pathway pres b1 b2 _) inlist
                            | and pwlist = (length pwlist,b2):inlist
                            | otherwise = inlist
                            where pwlist = searchPWs' pw
                        searchPWs' (Pathway [] b1 b2 _) = []
                        searchPWs' (Pathway (x:xs) b1 b2 _) 
                            | head x == '!' && opposite = go
                            | head x /= '!' && regular = go
                            | otherwise = stop
                            where 
                                query = st !! (getInd (x \\ "!"))
                                opposite = b1 == (0==query)
                                regular = b1 == (1==query)
                                stop = [True,False]-- This pathway doesn't apply
                                go = True : searchPWs' (Pathway xs b1 b2 [])
            rand xs = [xs++[0],xs++[1]]
            det0 xs = [xs++[0]]
            det1 xs = [xs++[1]]

buildBiGeneGraph :: ParseData -> StateTGraph
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

buildGeneGraph :: ParseData -> StateTGraph
buildGeneGraph pdata = mkGraph lnodes edges
    where   nodes = zipWith (\a (b,c)->(a,b,c)) [1..] (M.assocs pdata)
            lnodes = map (\(a,b,c)->(a,b)) nodes
            edges = concatMap edges1 nodes 
            edges1 (n,g,gi) = 
                zip3 (map n2N (allUpStream gi)) (repeat n) (repeat "")
            filt name = if head name == '!' then tail name else name
            n2N name = case lookup (filt name) name2NodeMap of
                Just x -> x
                Nothing -> error "no node mapping for this gene"
            name2NodeMap = zip (M.keys pdata) [1..] 
            allUpStream = concatMap (\(Pathway xs _ _ _)->xs) . pathways

stripTransNodes :: Graph gr => gr NodeInfo b -> gr NodeInfo b
stripTransNodes orig = foldr strip orig (nodes orig)
    where   strip node gr
                | (null $ inn orig node) && prob (lab gr node) /= 0.0 = error ("huh? non zero probability? "++show node)
                | null $ inn orig node = delNode node gr
                | otherwise = gr
            prob i = let Just (NodeInfo _ p) = i in p 
    

cus (NodeInfo name prob) = "[label=\""++name++
    "\",style=\"filled\",fillcolor=\"0.0 0.0 "++
    shade++"\""++",penwidth=\""++color++"\"]"
                where   shade = show $ 1 - (clamp.gammaCorr) prob
                        color = if name !! 4 == '1' then "1.0" else "1.0"
                        gammaCorr x = x**0.4
                        clamp x = if x>1 then 1 else if x<0 then 0 else x

genSSA :: Gr NodeInfo Double -> [Double]
genSSA g = map (/denom) $ foldl1 (zipWith (+)) ssaList
    where
        denom = fromIntegral n
        n = 8
        ssaList = map genSSA' glist
        glist = take n $ iterate fullIteration g

genSSA' :: Graph gr => gr NodeInfo Double -> [Double]
genSSA' g = foldr (zipWith (+)) (replicate len 0) filt2
    where   len = length $ name $ head filt1
            filt1 = map snd (labNodes g)
            filt2 = map weight filt1 
            weight (NodeInfo a pr) = map ((*pr).conv) a
            conv = (\x->if x=='1' then 1 else 0) 
            name (NodeInfo a pr) = a
            prob (NodeInfo a pr) = pr

genPrintSSA :: [Double] -> IO()
genPrintSSA vec = do
    mapM_ (printf "%7.3f") vec
    putStrLn ""

