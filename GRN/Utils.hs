-- |
-- Module    : GRN.Utils
-- Copyright : (c) 2011 Jason Knight
-- License   : BSD3
--
-- Maintainer  : jason@jasonknight.us
-- Stability   : experimental
-- Portability : portable
-- 

module GRN.Utils (
    dec2bin
  , bin2dec
  , b2i
  , chunks) where

import Data.List
import qualified Data.Vector as V

dec2bin :: Int -> Int -> [Int]
dec2bin len y = (replicate (len-(length res)) 0) ++ res -- need to add padding
    where   res = reverse $ dec2bin' y 
            dec2bin' 0 = []
            dec2bin' y = let (a,b) = quotRem y 2
                         in b : dec2bin' a
bin2dec :: [Int] -> Int
bin2dec st = foldl' (\a x->2*a+x) 0 st
--

b2i :: Bool -> Int
b2i True = 1
b2i _ = 0

i2b :: Int -> Bool
i2b 1 = True
i2b 0 = False
i2b _ = error "Ternary and higher quatizations are not yet supported."

chunks :: Int -> V.Vector a -> V.Vector (V.Vector a)
chunks n x = case V.splitAt n x of
              (a,b) | V.null a -> V.empty
                    | otherwise -> V.cons a (chunks n b)

