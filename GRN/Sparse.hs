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

module GRN.Sparse where

import Data.Vector.Unboxed (Vector)

import GRN.Types
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

dokToCSC :: DOK -> CSC
dokToCSC (DOK (m,n) d) = CSC (m,n) vals rows cumsum
  where
  cumsum = U.scanl' (+) 0 counts
  counts = (U.replicate n 0) U.// colGroups
  colGroups = map (\x->(head x,length x)) $ group $ map (snd.fst) xs
  vals = U.fromList (map snd xs)
  rows = U.fromList (map (fst.fst) xs)
  xs = sortBy sortCol $ M.toList d
  sortCol ((r1,c1),_) ((r2,c2),_)
    | c1 < c2 = LT
    | c1 > c2 = GT
    | c1 == c2 = compare r1 r2

multiplyCSCVM :: CSC -> Vector Double -> Vector Double
multiplyCSCVM CSC{..} x = {-# SCC "csc1" #-}  U.generate (U.length x) outer
  where outer i = {-# SCC "csc2" #-}  U.sum . U.map inner $ U.slice start (end-start) pre
          where inner j = {-# SCC "csc5" #-}  (cscValues ! j) * (x ! (rowIndices ! j))
                start   = colIndices ! i
                end     = colIndices ! (i+1)
                (!) a b = U.unsafeIndex a b
                pre     = U.enumFromN 0 (U.length rowIndices) :: U.Vector Int

simulateDOK :: Args String -> DOK -> SSD
simulateDOK args d@(DOK (m,n) _) = regSim (n1+n2) (U.replicate n (1.0/(fromIntegral n)))
      where 
            n1 = getRequiredArg args "n1" :: Int
            n2 = getRequiredArg args "n2"
            csc = dokToCSC d
            regSim 0 v = v
            regSim n v = regSim (n-1) (multiplyCSCVM csc v)

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

calcSSAsDOK :: Args String -> ParseData -> KmapSet -> GeneMC
calcSSAsDOK args p ks = zipNames . allgenes . ssaVec . simVec $ seedList
  where
    n3 = getRequiredArg args "n3" -- Number of simruns
    n4 = getRequiredArg args "n4" -- Starting Seed
    pass n f = last.take n.iterate f 

    -- Get a new randomized graph and simulate
    seedList = V.enumFromN (n4) (n3+n4)

    simVec :: V.Vector Int -> V.Vector SSD
    simVec xv = G.map (simulateDOK args.kmapToDOK ks) xv `using` (parVector 2)

    -- Now get the SSAs of these 
    ssaVec :: V.Vector SSD -> V.Vector SSA
    ssaVec xv = G.map ssdToSSA xv `using` (parVector 2)

    -- Basically a transpose, to get a list of SSA for each gene instead of
    -- each simulation run
    allgenes :: V.Vector SSA -> V.Vector (U.Vector Double)
    allgenes ssas = G.map (\x-> G.convert $ G.map (G.!x) ssas) (V.fromList [0..(length $ M.keys p)-1])
    zipNames :: V.Vector (U.Vector Double) -> M.Map Gene (U.Vector Double)
    zipNames xv = M.fromList $ zip (M.keys p) (G.toList xv)

ssdToSSA :: SSD -> SSA
ssdToSSA ssd = U.reverse $ U.generate ngenes count
  where
  n = U.length ssd
  ngenes = round $ logBase 2 (fromIntegral n)
  count i = U.sum $ U.ifilter (\ind _-> testBit ind i) ssd

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

--multiplyCRSMV :: CRS Double -> Vector Double -> Vector Double
--multiplyCRSMV CRS{..} x = generate (U.length x) outer
--  where outer i = U.sum . U.map inner $ U.enumFromN start (end-start)
--          where inner j = (crsValues ! j) * (x ! (colIndices ! j))
--                start   = rowIndices ! i
--                end     = rowIndices ! (i+1)
--                (!) a b = unsafeIndex a b


-- test = do
--     let --m = CRS (U.replicate n 1.0) (U.enumFromN 0 n) (U.enumFromN 0 (n+1))
--         d = DOK (n,n) (M.fromList [((7999,7999),1),((7999,4),1)])
--         ms = CSC (n,n) (U.replicate n 1.0) (U.enumFromN 0 n) (U.enumFromN 0 (n+1))
--         m = dokToCSC d
--         v = U.enumFromN 1.0 n
-- 
--     --print $ U.last $ multiplyCSRMV m v
--     print $ U.last $ multiplyCSCVM m v
