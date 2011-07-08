module GRN.Matrix where

import Data.Graph.Inductive hiding ((><))

import GRN.Parse
import GRN.StateTransition

import Data.Packed.Matrix
import Data.Packed.Vector

import qualified Data.Map as M

convert :: ColoredStateTGraph -> (Vector Double, Matrix Double)
convert gr = (map (\(_,NodeInfo _ x)-> x) (labNodes gr), mapMatrixWithIndex indfn empty)
    where 
        forward = M.fromList (zip [0..] (nodes gr))
        backward = M.fromList (zip (nodes gr) [0..])
        n = noNodes gr
        empty = (n><n) (repeat 0.0) :: Matrix Double
        indfn :: (Double,Double) -> Double -> Double
        indfn (i,j) _ 
            | elem j mappedouts = 1/(fromIntegral outd)
            | otherwise = 0.0
            where
                outd = outdeg gr node
                outes = map (\(_,x,_)->x) $ out gr node
                mappedouts = map (\x -> M.findWithDefault (-1) x backward) outes
                node = (nodes gr) !! (round i)



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



