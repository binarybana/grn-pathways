{-# Language RecordWildCards #-}
-- |
-- Module    : GRN.Sparse
-- Copyright : (c) 2011 Jason Knight
-- License   : BSD3
--
-- Maintainer  : jason@jasonknight.us
-- Stability   : experimental
-- Portability : portable
-- 
-- Provides sparse matrix capabilities for an attempt at speeding up the
-- simulation of the Markov Chain.
--

module GRN.Sparse (
    calcSSAsDOK
  , simulateDOKUnif 
  , simulateDOK 
  , normalizeSSD 
  , simulateCSC 
  , kmapToDOK 
  , ssdToSSA 
  , weightedExpansion
  , getNSamples
  , module GRN.SparseLib
  ) where

import GRN.Types
import GRN.SparseLib
import GRN.Utils
import GRN.StateTransition

import qualified Data.Vector as V
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as U

import Data.Vector.Strategies

import Control.Monad
import Data.List
import Data.Bits
import Data.Ord (comparing)
import qualified Data.Map as M
import System.Console.ParseArgs
import System.Random.MWC (withSystemRandom, initialize, uniformVector)
import Data.Word (Word32)

simulateDOKUnif :: Args String -> DOK -> SSD
simulateDOKUnif args d@(DOK (m,n) _) = simulateDOK args d (U.replicate n (1.0/(fromIntegral n)))

simulateDOK :: Args String -> DOK -> SSD -> SSD
simulateDOK args d@(DOK (m,n) _) start = regSim (n1+n2) start
      where 
            n1 = getRequiredArg args "n1" :: Int
            n2 = getRequiredArg args "n2"
            p = getRequiredArg args "perturb"
            csc = dokToCSC d
            regSim 0 v = v
            regSim n v
              | p < 1e-15 = regSim (n-1) (multiplyCSCVM csc v)
              | otherwise = regSim (n-1) 
                  (perturb p . multiplyCSCVM csc $ v)

-- | Perturb an SSD by a small amount
perturb :: Double -> SSD -> SSD
perturb p sd = snd . normalizeSSD $ added where
  added = U.map (+p) sd

normalizeSSD :: SSD -> (Double,SSA)
normalizeSSD ssd = let fac=U.sum ssd in (fac, U.map (*(1/fac)) ssd)

simulateCSC :: Args String -> CSC -> SSD -> SSD
simulateCSC args csc start = regSim (n1+n2) start
      where 
            n1 = getRequiredArg args "n1" :: Int 
            n2 = getRequiredArg args "n2"
            avg = getRequiredArg args "avg" :: Int
            --norm xv = let fac = (U.sum xv) in (fac,U.map (*(1/fac)) xv)
            regSim 0 v = v
            regSim n v = regSim (n-1) (multiplyCSCVM csc v)

kmapToDOK :: KmapSet -> Int -> DOK
kmapToDOK kset x = DOK (nStates,nStates) (M.fromList edges)
    where
        genes = M.keys kset
        nGenes = length genes
        nStates = 2^nGenes 
        intStates = [0..nStates-1]
        permedStates = map (dec2bin nGenes) intStates 
        edges = concatMap genEdges permedStates 
        genEdges :: [Int] -> [((Int,Int),Double)]
        genEdges st = map (\(val,to)-> ((from,bin2dec to),val)) tos
            where   
                from = bin2dec st
                tos = [(1.0,[])] >>= foldr (>=>) return (stateKmapLus x st genes kset)

calcSSAsDOK :: Args String -> ParseData -> KmapSet -> GeneMC
calcSSAsDOK args p ks = zipNames . allgenes . ssaVec . simVec $ seedList
  where
    n3 = getRequiredArg args "n3" -- Number of simruns
    n4 = getRequiredArg args "n4" -- Starting Seed
    pass n f = last.take n.iterate f 

    -- Get a new randomized graph and simulate
    seedList = V.enumFromN (n4) (n3+n4)

    simVec :: V.Vector Int -> V.Vector SSD
    simVec xv = G.map (simulateDOKUnif args.kmapToDOK ks) xv `using` (parVector 2)

    -- Now get the SSAs of these 
    ssaVec :: V.Vector SSD -> V.Vector SSA
    ssaVec xv = G.map ssdToSSA xv `using` (parVector 2)

    -- Basically a transpose, to get a list of SSA for each gene instead of
    -- each simulation run
    allgenes :: V.Vector SSA -> V.Vector (U.Vector Double)
    allgenes ssas = G.map (\x-> G.convert $ G.map (G.!x) ssas) (V.fromList [0..(length $ M.keys p)-1])
    zipNames :: V.Vector (U.Vector Double) -> M.Map Gene (U.Vector Double)
    zipNames xv = M.fromList $ zip (M.keys p) (G.toList xv)

-- | Expands out all determistic networks for a given likelihood network.
-- Number of networks grows superexponentially as (N+1)^(2^N) 
-- ie. 4, 81, 65536, 1.5e11...
weightedExpansion :: DOK -> V.Vector (Double, CSR)
weightedExpansion dk@(DOK (m,n) _) = wtdlist where
  csr = toCSR dk
  cols = csrcolIndices csr
  rows = csrrowIndices csr
  vals = csrValues csr
  perrow = (rows U.! 1) - (U.head rows)
  cartesian = G.sequence . chunks perrow . G.convert $ cols
  cartvals = G.map G.product . G.sequence . chunks perrow . G.convert $ vals
  newnnz = G.length . G.head $ cartesian
  newrows = U.enumFromN 0 newnnz
  newvals = U.replicate newnnz 1.0
  wtdlist = V.zip cartvals (G.map (\x -> CSR (m,n) newvals (G.convert x) newrows) cartesian)

chunks :: Int -> V.Vector a -> V.Vector (V.Vector a)
chunks n x = case V.splitAt n x of
              (a,b) | V.null a -> V.empty
                    | otherwise -> V.cons a (chunks n b)

getNSamples :: Int -> DOK -> IO (V.Vector (Double,CSR))
getNSamples samps dk@(DOK (m,n) _) = do
  let
    csr = toCSR dk
    cols = csrcolIndices csr
    rows = csrrowIndices csr
    vals = csrValues csr
    perrow = (rows U.! 1) - (U.head rows)
    wtdvals = G.convert $ G.zip vals cols
    rowchunks = chunks perrow $ wtdvals
    newnnz = G.length rowchunks 
    newrows = U.enumFromN 0 newnnz
    newvals = U.replicate newnnz 1.0
    go gen = do
      unifs <- uniformVector gen m :: IO (V.Vector Double)
      let
        pickedcols = V.zipWith pick unifs rowchunks 
        weight = G.product . fst . G.unzip $ pickedcols
        network = CSR (m,n) newvals (G.convert . snd . G.unzip $ pickedcols) newrows
      return (weight, network)
    pick :: Double -> V.Vector (Double,Int) -> (Double,Int)
    pick t vec 
      | G.null vec = error "Pick used with empty list"
      | t <= w = (w,x)
      | otherwise = pick (t-w) xs 
      where 
        (w,x) = G.head vec
        xs = G.tail vec
  -- Now ,get a list of generators, one for each sample
  gen <- withSystemRandom (\gen -> initialize =<< (uniformVector gen 256 :: IO (U.Vector Word32)))
  gens <- V.replicateM samps (initialize =<< (uniformVector gen 256 :: IO (U.Vector Word32)))
  G.mapM go gens

ssdToSSA :: SSD -> SSA
ssdToSSA ssd = U.reverse $ U.generate ngenes count
  where
  n = U.length ssd
  ngenes = round $ logBase 2 (fromIntegral n)
  count i = U.sum $ U.ifilter (\ind _-> testBit ind i) ssd

