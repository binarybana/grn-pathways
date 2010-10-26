module StateTransition where

import Data.Graph.Inductive
import Data.Graph.Inductive.Graph
import qualified Graphviz as Gv
import System.Cmd

import Control.Monad
import Data.List
import qualified Data.Map as M
import Parsea
import Data.Ratio

{--
myGraph :: Gr Char String
myGraph = mkGraph (zip [1..10] ['a'..]) (zip3 [1..10] (repeat 3) (repeat ""))

testDot :: IO ()
testDot = writeFile "tesths.dot" (graphviz myGraph "fgl" (8,8) (1,1) Portrait)
-- dot -Tpng -o tesths.png tesths.dot
--}

data NodeInfo = NodeInfo String Double 

instance Show NodeInfo where
    show (NodeInfo name num) = name ++ "\t" ++ (show num) ++ "\t"

type StateTGraph = Gr String String
type ColoredStateTGraph = Gr NodeInfo Double

dec2bin :: Int -> Int -> [Int]
dec2bin y len = (replicate (len-(length res)) 0) ++ res -- need to add padding
    where   res = reverse $ dec2bin' y
            dec2bin' 0 = []
            dec2bin' y = let (a,b) = quotRem y 2
                         in b : dec2bin' a
bin2dec :: [Int] -> Int
bin2dec st = foldl' (\a x->2*a+x) 0 st

fullIteration = iterateProbs2.iterateProbs 

iterateProbs2 :: ColoredStateTGraph -> ColoredStateTGraph
iterateProbs2 sgraph = gmap flowProbs sgraph
            where   flowProbs (adj1,node, (NodeInfo a oldprob) , adj2)=(newIns, node, NodeInfo a newProb, newOuts)  
                        where   newProb = sum $ map (\(_,p,a)->if p==node then a else 0.0) (labEdges sgraph)
                                newOuts = map (\(_,w)->(0.0,w)) adj2
                                newIns = map (\(_,w)->(0.0,w)) adj1
iterateProbs :: ColoredStateTGraph -> ColoredStateTGraph
iterateProbs sgraph = gmap flowProbs sgraph
            where   flowProbs  (adj1, node, (NodeInfo a oldprob), adj2) = (newIns, node, NodeInfo a 0.0, newOuts)  
                        where   newOuts = map (\(_,w)->(newSplitProb,w)) adj2
                                newIns = map updateIns adj1
                                updateIns (l,n) = let   c@(_,_,(NodeInfo _ p),_) = context sgraph n 
                                                        num = fromIntegral $ toInteger $ outdeg' c         
                                                in (p/num,n)
                                newSplitProb = oldprob/outdegree
                                outdegree = fromIntegral $ toInteger $ length adj2

initProbs :: StateTGraph -> ColoredStateTGraph
initProbs sgraph = gmap initProb (emap (\x->0.0) sgraph)
            where   initProb (adj1, node, a, adj2) = (adj1, node, NodeInfo a (1.0/num), adj2)
                    num = fromIntegral $ toInteger $ length $ nodes sgraph

buildStateTGraph :: ParseData -> StateTGraph
buildStateTGraph pdata = mkGraph permedNodeStates edges 
    where   permedStates = map (flip dec2bin nGenes) intStates
            permedNodeStates = zip intStates (map (concatMap show) permedStates)
            intStates = [0..(2^nGenes)-1]
            nGenes = length (M.keys pdata)
            indexNames = zip [1..] (M.keys pdata)

            edges = concatMap genEdges permedStates 
            genEdges :: [Int] -> [(Int,Int,String)]
            genEdges st = map (\to-> (from,bin2dec to,"")) tos
                where   from = bin2dec st
                        tos = [[]] >>= foldr (<=<) return (build $ length st)
                        build :: Int -> [([Int]->[[Int]])]
                        build 0 = [return.id]
                        build num = thisgene : build (num-1)
                            where thisgene = calcGene st num pdata
                                    --w = [[0],[1]] >>= app >>= app0 >>= app0
                                    --z = [[0],[1]] >>= foldr (<=<) return (replicate 3 app)

calcGene :: [Int] -> Int -> ParseData -> ([Int] -> [[Int]])
-- Given a state, and the gene number of interest, return a function
calcGene st prog pdata
    | length thisDeps == 0 = thisCurr
    | otherwise = transInfo
    where   names = M.keys pdata
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
                where   list = foldr searchPWs [] thisDeps
                        searchPWs depgene inlist = inlist ++ foldr match [] thisPways
                            where match (Pathway (pre:x) b1 b2 _) inlist = 
                                    if (not $ null x)
                                        then error "not implemented"
                                        else
                                            if pre==depgene && b1==((1 ==) $ st !! (getInd pre)) 
                                                then b2:inlist
                                                else inlist
            rand xs = [xs++[0],xs++[1]]
            det0 xs = [xs++[0]]
            det1 xs = [xs++[1]]

buildBiGeneGraph :: ParseData -> StateTGraph
buildBiGeneGraph pdata = mkGraph (map (\(a,b)->(b,a)) name2NodeMap) (edges ++ redges)
    where   edges = concatMap splitGis gis 
                where   gis = (M.elems pdata)
                        splitGis gi = concatMap (buildEdges (name gi)) (pathways gi)
            buildEdges down (Pathway ups b1 b2 _) = map buildEdge ups
                    where buildEdge up = case (b1,b2) of
                            (True,True)   -> (n2N up,n2N down,"")
                            (True,False)  -> (n2N up,negn2N down,"")
                            (False,True)  -> (negn2N up,n2N down,"")
                            (False,False) -> (negn2N up,negn2N down,"")
            negn2N name = n2N ('n':name)
            n2N name = case lookup name name2NodeMap of
                Just x -> x
                Nothing -> error "no node mapping for this gene"
            name2NodeMap = zip ((M.keys pdata) ++ ngenes) [1..]
            redges = zip3 (map n2N ngenes) (map (n2N.(\x->tail x))ngenes) (repeat "red")
            ngenes = map ('n':) (M.keys pdata)

buildGeneGraph :: ParseData -> StateTGraph
buildGeneGraph pdata = mkGraph lnodes edges
    where   nodes = zipWith (\a (b,c)->(a,b,c)) [1..] (M.assocs pdata)
            lnodes = map (\(a,b,c)->(a,b)) nodes
            edges = concatMap edges1 nodes 
            edges1 (n,g,gi) = zip3 (map n2N (allUpStream gi)) (repeat n) (repeat "")
            n2N name = case lookup name name2NodeMap of
                Just x -> x
                Nothing -> error "no node mapping for this gene"
            name2NodeMap = zip (M.keys pdata) [1..] 
            allUpStream = concatMap (\(Pathway xs _ _ _)->xs) . pathways

main = testNFKB

test0 = do 
        let nw = "nfkb-easy.pw"
        con <- readFile nw
        let p = parsePW con
        let g = (initProbs.buildStateTGraph) p
        return g

testNFKB = do
        let n = "nfkb-easy"
        let nw = n ++ ".pw"
        let nd = n ++ ".dot"
        let ns = n ++ ".svg"
        con <- readFile nw
        let p = parsePW con
        print $ M.keys p
        let initg = (initProbs.buildStateTGraph) p
        let gl = take 500 $ iterate fullIteration initg
        let g = last gl
        let prob (NodeInfo a pr) = pr
        print $ map (sum.(map (prob.snd)).labNodes) gl

        print g 
        writeFile nd (Gv.graphvizWithNodeFormatter cus g "fgl" (5,5) (1,1) Gv.Portrait)
        --rawSystem "neato" ["-Tsvg","-o","nfkb.svg","nfkb.dot","-Gsplines=True","-Goverlap=false"]
        rawSystem "dot" ["-Tsvg","-o",ns,nd]

cus (NodeInfo name prob) = "[label=\""++name++"\",style=\"filled\",fillcolor=\"0.0 0.0 "++shade++"\"]"
                where shade = show (1-prob)
