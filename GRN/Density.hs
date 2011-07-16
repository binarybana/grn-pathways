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

import Numeric.FFT.Vector.Unitary

import Statistics.KernelDensity
import Statistics.Sample
import Statistics.Constants
import Statistics.Distribution
import Statistics.Distribution.Normal
import Graphics.Gnuplot.Simple hiding (Points)

simRuns :: Args String -> ParseData -> IO () --U.Vector SSA
simRuns args p = do
  let 
    gen = gotArg args "generate"
    ks = buildKmaps p 
    gr = simulate args $ kmapToStateGraph p ks
    pass n f = last.take n.iterate f 
    n2 = getRequiredArg args "n2" -- Number of iterations per simulation
    n3 = getRequiredArg args "n3" -- Number of simruns
    n4 = getRequiredArg args "n4" -- Starting Seed

    simVec = V.iterateN n3 ((\(x,g)->(x,(pass n2 fullIteration g))).(reDrawEdges p ks)) (n4,gr)

    ssaVec :: V.Vector SSA
    ssaVec = G.map (genSSA.snd) simVec

    allgenes :: V.Vector (V.Vector Double)
    allgenes = G.map (\x-> G.map (G.!x) ssaVec) (V.fromList [0..(length $ M.keys p)-1])

    allvars = G.map stdDev allgenes

    step = 0.01 -- Step size
    nl = 101 -- Number of poitns
    smoothingVar = 0.01 -- Variance of smoothing gaussian

    points = U.enumFromStepN 0.00 step nl
    
    normSmooth mu = G.map (density $ normalDistr mu smoothingVar) points
    stdNorm = normSmooth 0.5

    myKernelEst :: V.Vector Double -> U.Vector Double
    myKernelEst sample = estimatePDF myKern (bandwidth gaussianBW sample) sample (Points points)

    myKern f h p v = exp (-0.5 * u * u) * g
        where u = (v - p) / h
              g = f * 0.5 * m_2_sqrt_pi * m_1_sqrt_2 

    allests :: V.Vector (U.Vector Double)
    allests = G.map myEst allgenes

    myEst xv  | stdDev xv < 0.001 = normalize.normSmooth $ (G.sum xv / (fromIntegral.G.length $ xv))
              | otherwise = normalize.clamp.convolve stdNorm.myKernelEst $ xv
              
    normalize xv = let scale = 1.0 / (step * G.sum xv) in G.map (*scale) xv
    clamp xv = G.map (\x->if x<0 then 0 else x) xv

    titles = V.fromList $ map (\x->defaultStyle {lineSpec = CustomStyle [LineTitle x]} ) (M.keys p)

    titleAll = V.zip titles

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

    prepForPlot = G.toList.titleAll.G.map (G.toList.G.zip points)

  mapM_ (printf "%7s") (M.keys p)
  putStrLn ""
  G.mapM_ printSSA ssaVec

  putStrLn ""
  putStrLn ""
  putStrLn "Std Deviations:"

  mapM_ (printf "%7s") (M.keys p)
  putStrLn ""
  G.mapM_ (printf "%7.4f") allvars

  putStrLn ""
  putStrLn ""
  putStrLn "Integrations:"

  mapM_ (printf "%7s") (M.keys p)
  putStrLn ""
  G.mapM_ (printf "%7.2f") (G.map ((*step).G.sum) allests)
  putStrLn ""
  
  plotPathsStyle [XRange (0.0,1.0)] $ prepForPlot allests

  when gen $ drawStateGraph (snd.G.last.G.take n3 $ simVec) args

