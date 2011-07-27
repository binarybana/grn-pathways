{-# Language RecordWildCards #-}
-- |
-- Module    : GRN.EM
-- Copyright : (c) 2011 Jason Knight
-- License   : BSD3
--
-- Maintainer  : jason@jasonknight.us
-- Stability   : experimental
-- Portability : portable
-- 
-- Expectation Maximization Algorithm at work
-- 

module GRN.EM where 

import GRN.StateTransition
import GRN.Render
import GRN.Sparse
import GRN.Types
import GRN.Utils

import qualified Data.Map as M
import Data.List
import Data.Maybe
import Control.Monad
import Text.Printf
import System.Console.ParseArgs

import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Generic as G

import Control.Parallel.Strategies
import Data.Vector.Strategies

import Data.Graph.Analysis.Algorithms.Common (componentsOf)

import Statistics.Distribution
import Statistics.Distribution.Normal

import Numeric.LBFGSB 

import Graphics.Gnuplot.Simple 

smoothingVar = 0.1 -- Variance of smoothing gaussian

ssaMap :: [Gene] -> SSA -> SSAMap
ssaMap glist s = M.fromList $ zip glist (U.toList s)

measureMap :: ParseData -> MeasureMap
measureMap p = M.map fromJust $ M.fromList $ filter (isJust.snd) $ map (\(x,gi)-> (x,measurement gi)) (M.toList p)

uncertainVals :: KmapSet -> [UncertainSpot]
uncertainVals ks = filter filtFun $ concatMap explode (M.elems ks)
  where explode (Kmap g upstream x) = map (\(pos,entry)->(g,pos,entry)) (M.toList x)

filtFun :: (Gene,[Int],Kentry) -> Bool
filtFun (g,pos,x) = case x of
  C _ -> True
  X -> True
  _ -> False

fillKentries :: [UncertainSpot] -> Theta -> KmapSet -> KmapSet
fillKentries spots mm ks = foldr fillKentry ks (zip spots (U.toList mm))

fillKentry :: (UncertainSpot,Double) -> KmapSet -> KmapSet
fillKentry ((g,pos,_),val) ks = M.adjust fill g ks 
  where fill (Kmap g up x) = Kmap g up (M.insert pos (V 0 val) x)

runAndSplit :: Args String -> KmapSet -> Attractors
runAndSplit args ks = map ((\(d,g)->(d, ssaMap (M.keys ks) (genSSA g), g)).normalizeGraph) raw
  where ssd = simulateDOK args.kmapToDOK ks $ 0
        raw = componentsOf $ convertProbsVG ssd stripG
        stripG = last . take 10.iterate stripTransNodes $ kmapToStateGraph ks

expectationRun :: EMData -> Theta -> Double
expectationRun EMData{..} new = -1/normSum * (G.sum (V.generate (length emAttractors) eachAttractor))
  where normSum = sum $ map (\(w,issam,gr)-> w * msLookup issam emMeasurements ) emAttractors
        newAgivenTheta = generateNewAttractors emArgs emUnfilledKMS emUncertainLocs new
        eachAttractor i = msLookup oldssam emMeasurements * w * log (newAgivenO * msLookup newssam emMeasurements)
          where (w,oldssam,gr) = emAttractors !! i 
                (newAgivenO, newssam, _) = newAgivenTheta !! i

--expectationRun EMData{..} new = G.map (*(1/normSum)) $ V.generate (length emAttractors) eachAttractor
--  where normSum = sum $ map (\(w,issam,gr)-> w * msLookup issam emMeasurements ) newAgivenTheta
--        newAgivenTheta = generateNewAttractors emArgs emUnfilledKMS emUncertainLocs new
--        eachAttractor i = msLookup newssam emMeasurements * newAgivenO
--          where (w,oldssam,gr) = emAttractors !! i 
--                (newAgivenO, newssam, _) = newAgivenTheta !! i

posteriorProb :: EMData -> Theta -> U.Vector Double
posteriorProb EMData{..} new = G.convert $ G.map (*(1/normSum)) (V.generate (length emAttractors) eachAttractor) 
  where normSum = sum $ map (\(w,issam,gr)-> w * msLookup issam emMeasurements ) newAgivenTheta
        newAgivenTheta = generateNewAttractors emArgs emUnfilledKMS emUncertainLocs new
        eachAttractor i = msLookup newssam emMeasurements * newAgivenO
          where (newAgivenO, newssam, _) = newAgivenTheta !! i

-- Given an unfilled KMS and a new theta set, fill the kmap, then generate an
-- Attractors list
generateNewAttractors :: Args String -> KmapSet -> [UncertainSpot] -> Theta -> Attractors
generateNewAttractors args ks us theta = runAndSplit args newks
  where newks = fillKentries us theta ks

-- Given a theta vector, ssa vector, and a set of measurements, multiply the
-- measurement lookups convolved with a measurement gaussian against each other
-- to get a joint lookup.
msLookup :: SSAMap -> MeasureMap -> Double
msLookup ssam mm = product $ M.elems $ M.mapWithKey draw mm
  where draw gene meas = truncatedNorm (ssam M.! gene) smoothingVar meas 

-- truncated Norm:: Mu  -> Var    -> Query
truncatedNorm :: Double -> Double -> Double -> Double
truncatedNorm mu var pt = 1/var * std ((pt-mu)/var) / denom
  where std = density standard 
        cstd = cumulative standard 
        denom = cstd ((1-mu)/var) - cstd ((0-mu)/var)

emRun :: Args String -> ParseData -> IO () 
emRun args p = do
  let
    ks = buildKmaps p 
    ms = measureMap p
    uvals = uncertainVals ks


    n = length uvals

    midt = U.replicate n 0.5
    --mapt = map (U.fromList.map fromIntegral.dec2bin n) [0..2^n-1]
    mapt = map (U.replicate n) [0,0.9,1.0]
    --mapt = map (U.fromList) [replicate n 0.0, replicate n 0.001]

    emdat = EMData (runAndSplit args (fillKentries uvals midt ks)) ks args ms uvals
    

    --prVal = map (expectationRun emdat) mapt
    prVal = map (posteriorProb emdat) mapt 
    prVal2 = map (expectationRun emdat) mapt `using` parList rdeepseq

    --filt = U.map (\x->if x<1e-3 then 0 else x)
    printList x = (map (printf "%7.3f") (G.toList x)) ++ [putStrLn ""]
    elOpt = optimize (expectationRun emdat) (U.replicate n 0.5) 
    ag = approxGrad (expectationRun emdat) (U.replicate n 0.5)
    ag2 = approxGrad2 (expectationRun emdat) (U.replicate n 0.5)
    xs = [0,0.1..1]

  print n
  --plotFunc3d [] [] xs xs (\x y -> expectationRun emdat (U.fromList [x,0.5,y]))
  --mapM_ (sequence.printList) [ag,ag2]
  --print $ last $ emAttractors emdat
  --putStrLn ""
  --putStrLn "Posterior Probs:"
  mapM_ (sequence.printList) prVal
  --putStrLn ""
  --putStrLn "Expectations:"
  --mapM_ (printf "%7.3f\n") prVal2
  --putStrLn ""
  --putStrLn "Optimal Posterior Prob:"
  --mapM_ (sequence.printList) [posteriorProb emdat elOpt]
  
  --putStrLn ""
  --printf "Optimal Expectation: %5.2f\n" (expectationRun emdat elOpt)
  --putStrLn ""
  --putStrLn "**** Optimal Theta ****"
  --mapM_ (sequence.printList) [elOpt]
  --
  --putStrLn "***********************\n"
 

