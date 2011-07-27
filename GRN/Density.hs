-- |
-- Module    : GRN.Density
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
import GRN.Sparse
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

myKernelEst :: U.Vector Double -> U.Vector Double
myKernelEst sample = estimatePDF myKern (bandwidth gaussianBW sample) sample (Points points)

myKern f h p v = exp (-0.5 * u * u) * g
    where u = (v - p) / h
          g = f * 0.5 * m_2_sqrt_pi * m_1_sqrt_2 

myEst :: U.Vector Double -> U.Vector Double
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

calcEsts :: GeneMC -> AttractorDensity
calcEsts = M.map myEst

simRuns :: Args String -> ParseData -> IO () --U.Vector SSA
simRuns args p = do
  let
    titles = map (\x->defaultStyle {lineSpec = CustomStyle [LineTitle x]} ) (M.keys p)
    titleAll = zip titles
    prepForPlot = titleAll.map (G.toList.G.zip points.snd).M.toList
    ks = buildKmaps p 
    dm = getRequiredArg args "dmode" -- Number of simruns

    vals = concatMap (\(Kmap _ _ x)->M.elems x) (M.elems ks)
    unknowns = [C x | C x <- vals] ++ [X | X <- vals]

  print unknowns

  when (dm == "m") $ do  
    let
      allssas = calcSSAsDOK args p ks
      allests = calcEsts allssas
    --print allssas
    plotPathsStyle [XRange (0.0,1.0)] (prepForPlot allests)


  when (dm == "g") $ do
    let
      grComps = componentsOf $ simulate args $ kmapToStateGraph ks
      grNComps = map normalizeGraph grComps
      allests = map (\(d,g)->(d, g, calcEsts.calcSSAs args p ks $ g)) grNComps

    mapM_ (plotPathsStyle [XRange (0.0,1.0)] . prepForPlot . (\(_,_,x)->x)) allests
    --when gen $ drawStateGraph snd.G.head $ simVec) args

