-- |
-- Module    : Main
-- Copyright : (c) 2011 Jason Knight
-- License   : BSD3
--
-- Maintainer  : jason@jasonknight.us
-- Stability   : experimental
-- Portability : portable
-- 
-- Simulates many different random values in the unknown Kmap entries and
-- estimates the density of this distribution for each gene.
--
-- Then we will convolve this with the 'measurement kernel' to provide a non
-- dirac prior for every gene (even deterministic ones).
-- 

module GRN.Density where 

import GRN.StateTransition
import GRN.Render
import GRN.Types

import qualified Data.Map as M
import Data.List
import Control.Monad
import Text.Printf
import System.Console.ParseArgs

import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Generic as G

import Control.Parallel.Strategies

import Numeric.FFT.Vector.Unitary

import Data.Graph.Analysis.Algorithms.Common (componentsOf)

import Statistics.KernelDensity
import Statistics.Sample
import Statistics.Constants
import Statistics.Distribution
import Statistics.Distribution.Normal
import Graphics.Gnuplot.Simple hiding (Points)



step = 0.01 -- Step size
nl = 101 -- Number of poitns
smoothingVar = 0.001 -- Variance of smoothing gaussian

points = U.enumFromStepN 0.00 step nl

normSmooth mu = G.map (density $ normalDistr mu smoothingVar) points
stdNorm = normSmooth 0.5

myKernelEst :: V.Vector Double -> U.Vector Double
myKernelEst sample = estimatePDF myKern (bandwidth gaussianBW sample) sample (Points points)

myKern f h p v = exp (-0.5 * u * u) * g
    where u = (v - p) / h
          g = f * 0.5 * m_2_sqrt_pi * m_1_sqrt_2 


myEst xv  | stdDev xv < 0.001 = normalize.normSmooth $ (G.sum xv / (fromIntegral.G.length $ xv))
          | otherwise = normalize.clamp.convolve stdNorm.myKernelEst $ xv
          
normalize xv = let scale = 1.0 / (step * G.sum xv) in G.map (*scale) xv
clamp xv = G.map (\x->if x<0 then 0 else x) xv


-- Assume that both xs and ys are the same size
convolve xs ys = G.slice start l uncentered
  where
    xsF = run dftR2C (xs G.++ zeros)
    ysF = run dftR2C (ys G.++ zeros)
    l = G.length xs
    start = floor $ (fromIntegral l)/2
    zeros = U.replicate l 0.0
    mult = G.zipWith (*) xsF ysF
    uncentered = run dftC2R mult



calcSSAs :: Args String -> ParseData -> [V.Vector (V.Vector Double)]
calcSSAs args p = map (simSSAs.snd.normalizeGraph) grComps
  where
    gen = gotArg args "generate"
    ks = buildKmaps p 
    grComps = componentsOf $ simulate args $ kmapToStateGraph p ks
    pass n f = last.take n.iterate f 
    n2 = getRequiredArg args "n2" -- Number of iterations per simulation
    n3 = getRequiredArg args "n3" -- Number of simruns
    n4 = getRequiredArg args "n4" -- Starting Seed

    seedList :: ColoredStateGraph -> V.Vector (Int,ColoredStateGraph)
    seedList gr = V.fromList $ zip [n4..(n4+n3)] (repeat gr)

    simVec :: V.Vector (Int,ColoredStateGraph) -> V.Vector (Int,ColoredStateGraph)
    simVec = G.map ((\(x,g)->(x,(pass n2 fullIteration g))).(reDrawEdges p ks)) 

    ssaVec :: V.Vector (Int,ColoredStateGraph) -> V.Vector SSA
    ssaVec = V.fromList . withStrategy (parList rpar) . G.toList . G.map (genSSA.snd)

    allgenes ssas = G.map (\x-> G.map (G.!x) ssas) (V.fromList [0..(length $ M.keys p)-1])

    simSSAs = allgenes . ssaVec . simVec . seedList 

calcEsts :: V.Vector (V.Vector Double) -> V.Vector (U.Vector Double)
calcEsts = G.map myEst

simRuns :: Args String -> ParseData -> IO () --U.Vector SSA
simRuns args p = do
  let
    titles = V.fromList $ map (\x->defaultStyle {lineSpec = CustomStyle [LineTitle x]} ) (M.keys p)
    titleAll = V.zip titles
    
    prepForPlot = G.toList.titleAll.G.map (G.toList.G.zip points)
    allgenes = calcSSAs args p
    allests = map calcEsts allgenes

  mapM_ (plotPathsStyle [XRange (0.0,1.0)] . prepForPlot) allests

  --when gen $ drawStateGraph snd.G.head $ simVec) args

