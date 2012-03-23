{-# Language RecordWildCards #-}
-- |
-- Module    : GRN.SparseLib
-- Copyright : (c) 2012 Jason Knight
-- License   : BSD3
--
-- Maintainer  : jason@jasonknight.us
-- Stability   : experimental
-- Portability : portable
-- 
-- Provides a sparse matrix library.
--

module GRN.SparseLib where

import Data.Vector.Unboxed (Vector)

import qualified Data.Vector as V
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as U

import qualified Data.Map as M
import Data.List (group, groupBy, sortBy) 

-- | A dictionary of values of (row,col) order 
data DOK = DOK (Int,Int) (M.Map (Int,Int) Double)
  deriving (Show,Eq,Ord)

-- | A compressed row storage (CSR) sparse matrix.
data CSR = CSR {
      csrmatShape :: (Int,Int)
    , csrValues :: U.Vector Double
    , csrcolIndices :: U.Vector Int
    , csrrowIndices :: U.Vector Int
    } deriving (Show,Eq,Ord)

-- | A compressed column storage (CSC) sparse matrix.
data CSC = CSC {
      cscmatShape :: (Int,Int)
    , cscValues :: U.Vector Double -- Values
    , cscrowIndices :: U.Vector Int -- length = nnz, 
    , csccolIndices :: U.Vector Int -- length = ncols+1
    } deriving (Show,Eq,Ord)

class SparseMatrix a where
  toCSC :: a -> CSC
  toCSR :: a -> CSR
  toDOK :: a -> DOK

instance SparseMatrix DOK where
  toDOK = id
  toCSC = dokToCSC
  toCSR = dokToCSR

instance SparseMatrix CSC where
  toDOK = cscToDOK
  toCSC = id
  toCSR = dokToCSR . cscToDOK

instance SparseMatrix CSR where
  toDOK = csrToDOK
  toCSC = dokToCSC . csrToDOK
  toCSR = id

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

dokToCSR :: DOK -> CSR
dokToCSR (DOK (m,n) d) = CSR (m,n) vals cols cumsum
  where
  cumsum = U.scanl' (+) 0 counts
  counts = (U.replicate n 0) U.// rowGroups
  rowGroups = map (\x->(head x,length x)) $ group $ map (fst.fst) xs
  vals = U.fromList (map snd xs)
  cols = U.fromList (map (snd.fst) xs)
  xs = sortBy sortRow $ M.toList d
  sortRow ((r1,c1),_) ((r2,c2),_)
    | r1 < r2 = LT
    | r1 > r2 = GT
    | r1 == r2 = compare c1 c2

cscToDOK :: CSC -> DOK
cscToDOK (CSC (m,n) vals rows cumsum) = DOK (m,n) d where
  d = M.fromList . U.toList . U.zip inds $ vals
  inds = U.zip rows (inverseCumSum cumsum)

csrToDOK :: CSR -> DOK
csrToDOK (CSR (m,n) vals cols cumsum) = DOK (m,n) d where
  d = M.fromList . U.toList . U.zip inds $ vals
  inds = U.zip (inverseCumSum cumsum) cols

inverseCumSum :: U.Vector Int -> U.Vector Int
inverseCumSum vec = G.concatMap (\(i,x) -> G.replicate x i) . G.indexed . 
                      G.zipWith (-) (G.tail vec) $ vec

multiplyCSCVM :: CSC -> Vector Double -> Vector Double
multiplyCSCVM CSC{..} x = {-# SCC "csc1" #-}  U.generate (U.length x) outer
  where outer i = {-# SCC "csc2" #-}  U.sum . U.map inner $ U.slice start (end-start) pre
          where inner j = {-# SCC "csc5" #-}  (cscValues ! j) * (x ! (cscrowIndices ! j))
                start   = csccolIndices ! i
                end     = csccolIndices ! (i+1)
                (!) a b = U.unsafeIndex a b
                pre     = U.enumFromN 0 (U.length cscrowIndices) :: U.Vector Int


checkStochasticDOK :: DOK -> Bool
checkStochasticDOK (DOK (n,_) mat) = all close rows where
  close x = if x-1 < 1e-12 then True else False
  rows = map (sum . map snd) $ groupBy (\((x1,_),_) ((x2,_),_) -> x1==x2) $ M.assocs mat

--multiplyCRSMV :: CRS Double -> Vector Double -> Vector Double
--multiplyCRSMV CRS{..} x = generate (U.length x) outer
--  where outer i = U.sum . U.map inner $ U.enumFromN start (end-start)
--          where inner j = (crsValues ! j) * (x ! (crcolIndices ! j))
--                start   = csrowIndices ! i
--                end     = csrowIndices ! (i+1)
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
