{-# LANGUAGE BangPatterns #-}

module GRN.Fuzzed where

import Data.Graph.Inductive hiding ((><))

import GRN.Parse
import GRN.StateTransition

import Data.Packed.Matrix
import Data.Packed.Vector

import qualified Data.Map as M

data EdgeInfo = EdgeInfo Double Double
type FuzzyStateGraph = Gr NodeInfo EdgeInfo

convert :: ColoredStateTGraph -> FuzzyStateGraph
convert gr = gmap initEdgeProb (emap (\x->EdgeInfo 0 x) gr)
    where 
        initProb (adj1, node, a, adj2) = 
            (newadj1, node, a, newadj2)
        procEdges xs = map som xs



testFun = do
    con <- readFile "./pws/nfkb-super.pw"
    print $ testPure con
    
testPure con = finalGraph
    where
        start = parsePW $ con ++ (unlines.words $ "tnf=1 ltbr=0 lps=0") 
        initg = (initProbs.buildStateTGraph) start
        finalGraph = (pass 100 fullIteration).
                    (pass 20 (stripTransNodes.fullIteration)) $ initg
        pass n f = last.take n.iterate f 



