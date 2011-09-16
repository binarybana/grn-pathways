-- |
-- Module    : GRN.Uncertainty
-- Copyright : (c) 2011 Jason Knight
-- License   : BSD3
--
-- Maintainer  : jason@jasonknight.us
-- Stability   : experimental
-- Portability : portable
-- 
-- Print out a set of boolean networks that represent the uncertainty class.
-- 

module GRN.Uncertainty where 

import GRN.StateTransition
import GRN.Render
import GRN.Sparse
import GRN.Types
import GRN.Utils
import GRN.EM

import qualified Data.Map as M
import Data.List
import Data.Char
import Data.Bits
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
import Text.PrettyPrint


-- Strips nodes that violate the input node values
stripInputNodes :: ParseData -> ColoredStateGraph -> ColoredStateGraph
stripInputNodes pdata orig = foldr strip orig (nodes orig)
    where   strip node gr
                | compareBin imask ivals node = gr
                | otherwise = delNode node gr
            inames = sort $ inputGeneNames pdata
            ivals = map (b2i . fromJust . knockout . fromJust . flip M.lookup pdata) inames
            imask = inputGenes pdata

-- | Returns true if the node matches the input parameters as determined by the
-- bitmask, and the values that go in those entries
compareBin :: [Bool] -> [Int] -> Int -> Bool
compareBin bm vals state = snd $ foldl helper (vals,True) (zip bm (dec2bin (length bm) state))
  where
    helper ([],judge) _ = ([], judge)
    helper ((x:xs),judge) (bmv, sv)
        | not judge = ([],False)
        | not bmv = (x:xs, True)
        | sv /= x = ([],False)
        | otherwise = (xs, True)

-- | Filter a binary number by a 'bitmask' of bools
binaryFilter :: [Bool] -> [Int] ->  [Int]
binaryFilter bm num = fst . unzip . filter (not . snd) $ zip num bm

-- | Find the 'bitmask' for input genes
inputGenes :: ParseData -> [Bool]
inputGenes pdata = map (inputGene pdata) (M.keys pdata)

-- | Find the names of the input genes
inputGeneNames :: ParseData -> [Gene]
inputGeneNames pdata = filter (inputGene pdata) (M.keys pdata)

-- | Use this to find input genes by the parse data.
inputGene :: ParseData -> Gene -> Bool
inputGene pdata g = null $ depends gi
  where
    Just gi = M.lookup g pdata

uncertaintyPrint :: Args String -> ParseData -> IO () 
uncertaintyPrint args p = do
  let
    {-gen = gotArg args "generate"-}
    ks = buildKmaps p 
    ms = measureMap p
    n1 = getRequiredArg args "n1" :: Int
    uvals = uncertainVals ks


    ng = M.size ks
    n = length uvals

    hypercorners = map (U.fromList.map fromIntegral.dec2bin n) [0..2^n-1]

    {-printList :: [U.Vector Double] -> String-}
    {-printList x = concat $ intersperse "\n" (map (concatMap (printf "%7.3f").G.toList) $ x)-}

    {---startTheta = U.fromList [0.379,0.5,1.0]-}
    {-startTheta = U.replicate n 0.5-}

    {-newks = fillKentries uvals ks startTheta-}
    {-stripG = componentsOf . last . take 10.iterate stripTransNodes $ kmapToStateGraph newks-}
    {-stripNodes = map (U.fromList.sort.nodes) stripG-}
    {-startAtracs = runAndSplit args newks stripNodes-}
    {-norms = normSum startAtracs ms-}

    {-emdat = EMData startAtracs ks stripG stripNodes norms args ms uvals-}
    
    {--- type fullEM :: EMData -> Int -> Theta -> [(Theta,EMData)]-}
    {-optimumTheta = fst . head . fullEM emdat 10 $ startTheta-}

    kmaps = map (fillKentries uvals ks) hypercorners
    graphs = map (stripInputNodes p . kmapToStateGraph) kmaps
    edgesl = map edges graphs

    bp = bin2dec . binaryFilter (inputGenes p) . dec2bin ng
    filtedges = map (map (\(x,y) -> (bp x, bp y))) 
    
    fedges = filtedges edgesl
    
  pprint fedges
  {-print "CHECK:"-}
  {-print $ map (length . nub . map fst) fedges-}
  {-print ""-}

  return ()
  

class Pretty a where
  prettyPrint :: a -> IO()
  prettyPrint = putStrLn . renderStyle (Style LeftMode 80 0) . prettyShow
  prettyShow :: a -> Doc

pprint :: Pretty a => a -> IO()
pprint = prettyPrint

instance Pretty Int where
  prettyShow = text . show
instance Pretty a => Pretty [a] where
  prettyShow = brackets . hcat . punctuate comma . map prettyShow
instance (Pretty a, Pretty b) => Pretty (a,b) where
  prettyShow (a,b) = prettyShow a <> colon <> prettyShow b 

