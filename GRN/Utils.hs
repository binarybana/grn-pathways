-- |
-- Module    : GRN.Utils
-- Copyright : (c) 2011 Jason Knight
-- License   : BSD3
--
-- Maintainer  : jason@jasonknight.us
-- Stability   : experimental
-- Portability : portable
-- 

module GRN.Utils (dec2bin, bin2dec) where
import Data.List

dec2bin :: Int -> Int -> [Int]
dec2bin len y = (replicate (len-(length res)) 0) ++ res -- need to add padding
    where   res = reverse $ dec2bin' y 
            dec2bin' 0 = []
            dec2bin' y = let (a,b) = quotRem y 2
                         in b : dec2bin' a
bin2dec :: [Int] -> Int
bin2dec st = foldl' (\a x->2*a+x) 0 st
