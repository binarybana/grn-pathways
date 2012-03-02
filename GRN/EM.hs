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
import Data.Char
import Data.Maybe
import Control.Monad
import Text.Printf
import System.Console.ParseArgs

import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Generic as G

import Control.Parallel.Strategies
import Data.Vector.Strategies

import Data.Graph.Inductive
import Data.Graph.Analysis.Algorithms.Common (componentsOf)

import Statistics.Distribution
import Statistics.Distribution.Normal

import Numeric.LBFGSB 

import Graphics.Gnuplot.Simple 

smoothingVar = 0.2 -- Variance of smoothing gaussian

ssaMap :: [Gene] -> SSA -> SSAMap
ssaMap glist s = M.fromList $ zip glist (U.toList s)

-- | Returns a map Gene->Double of any gene that has a measurement recorded.
measureMap :: ParseData -> MeasureMap
measureMap p = M.map fromJust $ M.fromList $ filter (isJust.snd) $ map (\(x,gi)-> (x,measurement gi)) (M.toList p)

-- | Returns a list of Uncertain locations.
uncertainVals :: KmapSet -> [UncertainSpot]
uncertainVals ks = filter filtFun $ concatMap explode (M.elems ks)
  where explode (Kmap g upstream x) = map (\(pos,entry)->(g,pos,entry)) (M.toList x)

filtFun :: (Gene,[Int],Kentry) -> Bool
filtFun (g,pos,x) = case x of
  C _ -> True
  X -> True
  _ -> False

fillKentries :: [UncertainSpot] -> KmapSet -> Theta -> KmapSet
fillKentries spots ks mm = foldr fillKentry ks (zip spots (U.toList mm))

fillKentry :: (UncertainSpot,Double) -> KmapSet -> KmapSet
fillKentry ((g,pos,_),val) ks = M.adjust fill g ks 
  where fill (Kmap g up x) = Kmap g up (M.insert pos (V 0 val) x)

-- truncated Norm:: Mu  -> Var    -> Query
truncatedNorm :: Double -> Double -> Double -> Double
truncatedNorm mu var pt = 1/var * std ((pt-mu)/var) / denom
  where std = density standard 
        cstd = cumulative standard 
        denom = cstd ((1-mu)/var) - cstd ((0-mu)/var)

-- Given a theta vector, ssa vector, and a set of measurements, multiply the
-- measurement lookups convolved with a measurement gaussian against each other
-- to get a joint lookup.
msLookup :: SSAMap -> MeasureMap -> Double
msLookup ssam mm = product $ M.elems $ M.mapWithKey draw mm
  where draw gene meas = truncatedNorm (ssam M.! gene) smoothingVar meas 

runAndSplit :: Args String -> KmapSet -> Masks -> Attractors
runAndSplit args ks masks = map ((\(d,s)->(d, ssaMap (M.keys ks) (ssdToSSA s))).normalizeSSD) raw
  where base = simulateCSC args csc (U.replicate n (1/(fromIntegral n)))
        dok = kmapToDOK ks 0
        csc = dokToCSC dok
        n = 2^(M.size ks) 
        raw = map (filtSSD base) masks
        filtSSD base mask = zero `U.update` inds 
          where zero = U.replicate n 0.0
                inds = U.map (\i->(i,base U.! i)) mask

runAndSplit2 :: Args String -> KmapSet -> Masks -> Attractors
runAndSplit2 args ks masks = map ((\(d,g)->(d, ssaMap (M.keys ks) (genSSA g))).normalizeGraph) raw
  where ssd = simulateDOKUnif args.kmapToDOK ks $ 0
        raw = componentsOf $ convertProbsVG ssd stripG
        stripG = last . take 10.iterate stripTransNodes $ kmapToStateGraph ks

-- Given an unfilled KMS and a new theta set, fill the kmap, then generate an
-- Attractors list
generateNewAttractors :: Args String -> KmapSet -> Masks -> [UncertainSpot] -> Theta -> Attractors
generateNewAttractors args ks masks us theta = runAndSplit args newks masks
  where newks = fillKentries us ks theta

expectationRun :: EMData -> Theta -> Double
expectationRun EMData{..} new = emNormSum * (G.sum (V.generate (length emAttractors) eachAttractor))
  where newAgivenTheta = generateNewAttractors emArgs emUnfilledKMS emCompNodes emUncertainLocs new
        eachAttractor i = msLookup oldssam emMeasurements * w * log (newAgivenO * msLookup newssam emMeasurements)
          where (w,oldssam) = emAttractors !! i 
                (newAgivenO, newssam) = newAgivenTheta !! i

normSum :: Attractors -> MeasureMap -> Double
normSum atracs mm = (-1) / (sum $ map (\(w,issam)-> w * msLookup issam mm) atracs)

posteriorProb :: EMData -> Theta -> U.Vector Double
posteriorProb EMData{..} new = G.convert $ G.map (*(1/normSum)) (V.generate (length emAttractors) eachAttractor) 
  where normSum = sum $ map (\(w,issam)-> w * msLookup issam emMeasurements ) newAgivenTheta
        newAgivenTheta = generateNewAttractors emArgs emUnfilledKMS emCompNodes emUncertainLocs new
        eachAttractor i = msLookup newssam emMeasurements * newAgivenO
          where (newAgivenO, newssam) = newAgivenTheta !! i

fullEM :: EMData -> Int -> Theta -> [(Theta,EMData)]
fullEM emd 0 oldtheta = [(oldtheta,emd)]
fullEM emd maxiter oldtheta 
  | thetaClose = (oldtheta,emd) : [(newtheta,newemd)]
  | otherwise = (oldtheta,emd) : fullEM newemd (maxiter-1) newtheta
  where newemd = emd { emAttractors=newAtracs, emNormSum=newNormSum } 
        newAtracs = generateNewAttractors (emArgs emd) (emUnfilledKMS emd) 
                        (emCompNodes emd) (emUncertainLocs emd) newtheta
        newNormSum = normSum newAtracs (emMeasurements emd)
        newtheta = optimize (expectationRun emd) oldtheta
        thetaClose = (<1e-5).U.sum.U.map (**2).U.zipWith (-) newtheta $ oldtheta

--data EMData = EMData {
--      emAttractors :: Attractors
--    , emUnfilledKMS :: KmapSet
--    , emComps :: [ColoredStateGraph]
--    , emCompNodes :: Masks
--    , emNormSum :: Double
--    , emArgs :: Args String
--    , emMeasurements :: MeasureMap
--    , emUncertainLocs :: [UncertainSpot] 
--    } 

emRun :: Args String -> ParseData -> IO () 
emRun args p = do
  let
    gen = gotArg args "generate"
    ks = buildKmaps p 
    ms = measureMap p
    uvals = uncertainVals ks


    ng = M.size ks
    n = length uvals
    printList :: [U.Vector Double] -> String
    printList x = concat $ intersperse "\n" (map (concatMap (printf "%7.3f").G.toList) $ x)

    --startTheta = U.fromList [0.379,0.5,1.0]
    startTheta = U.replicate n 0.5

    newks = fillKentries uvals ks startTheta
    stripG = componentsOf . last . take 10.iterate stripTransNodes $ kmapToStateGraph newks
    stripNodes = map (U.fromList.sort.nodes) stripG
    startAtracs = runAndSplit args newks stripNodes
    norms = normSum startAtracs ms

    emdat = EMData startAtracs ks stripG stripNodes norms args ms uvals
    

    mapt = map (U.fromList.map fromIntegral.dec2bin n) [0..2^n-1]
    --mapt = map (U.replicate n) [0,0.1,0.9,1.0]
    --mapt = map (U.fromList) [replicate n 0.0, replicate n 0.001]
    --prVal = map (expectationRun emdat) mapt
    prVal = map (posteriorProb emdat) mapt 
    prVal2 = map (expectationRun emdat) mapt 

    elOpt = optimize (expectationRun emdat) (U.replicate n 0.5) 
    path = fullEM emdat 10 startTheta

    xs = clamped [0,0.05..1]

    clamped xs = map (\x-> if (x>1.0) then 1.0 else (if (x<0.0) then 0.0 else x)) xs

    printSeg :: String -> String -> IO ()
    printSeg title st = do
          let l = 2 + length title
              upperTitle = map toUpper title
              width = 80
              half = (width-l) `div` 2
              stars = replicate half '*'
          printf "\n%s %s %s\n\n%s\n\n%s\n\n" stars upperTitle stars st (replicate width '*')


  --print n

  --forM_ mapt $ \t -> do
  --  let path = fullEM emdat 10 t
  printSeg "theta path" (printList.map (\(t,e)->t) $ path)
  printSeg "posterior probs" (printList.map (\(t,e)->posteriorProb e t) $ path)

  --printSeg "attractors" (concatMap (printf "%7s".map intToDigit.dec2bin ng.G.head) stripNodes)
  printSeg "attractors" (concatMap (printf "%7d".G.head) stripNodes)
  print $ M.keys ks

  --printf "%d Evaluations\n" ((length xs)^2)
  plotFunc3d [XLabel "X", YLabel "Y", XRange (0,1), YRange (0,1)] [] xs xs 
     (\x y -> (-1)*(expectationRun emdat (U.fromList [0.5,x,y])))
  
  --putStrLn "Expectations:"
  --mapM_ (printf "%7.3f") prVal2
  --putStrLn ""

  --print elOpt

  


  --putStrLn "Posterior Probs:"
  --mapM_ (sequence.printList.U.fromList.sort.U.toList) prVal
  --putStrLn ""
  --mapM_ (printf "%5.2f") (map U.sum prVal)
  --putStrLn ""

  --putStrLn ""
  --putStrLn ""
  --print $ generateNewAttractors args ks stripNodes uvals (last mapt)
  --
  --putStrLn ""
  --putStrLn "Expectations:"
  --mapM_ (printf "%7.3f\n") prVal2
  --putStrLn ""

  --putStrLn ""
  --putStrLn ""
  --putStrLn "**** Optimal Posterior Probs ****"
  --mapM_ (sequence.printList) [posteriorProb emdat elOpt]
  --mapM_ (printf "%7.3d".(U.head)) stripNodes
  --putStrLn ""
  --putStrLn "*********************************"
  --
  --putStrLn ""
  --printf "Optimal Expectation: %5.2f\n" (expectationRun emdat elOpt)
  --putStrLn ""
  --putStrLn "**** Optimal Theta ****"
  --mapM_ (sequence.printList) [elOpt]

  printSeg "uncertain locs" (show uvals)
  printSeg "measurements" (show ms)
 
  when gen $ do
    let newks2 = fillKentries uvals ks (fst.last $ path)
        finalGraph = simulate args $ kmapToStateGraph newks2
        stripG2 = componentsOf . last . take 10.iterate stripTransNodes $ finalGraph
        stripNodes2 = map (U.fromList.sort.nodes) stripG2

    printSeg "final attractors" (concatMap (printf "%7d".G.head) stripNodes2)
    drawStateGraph finalGraph args

    
    --let newks2 = fillKentries uvals ks (fst.head $ path)
    --    interGraph = kmapToStateGraph newks2
    --    ssd = simulateDOK args $ kmapToDOK newks2

    ---- WATCH OUT TODO I just found a huge bug: StateTransition.simulate does
    ---- not respect edge weights when moving probabilities around

