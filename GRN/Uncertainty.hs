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

import Graphics.Gnuplot.Simple 
import Text.PrettyPrint

-- | Collapses SSD onto space of non-input nodes
collapseSSD :: ParseData -> SSD -> SSD
collapseSSD pdata orig = snd $ foldl' splitAndSum (pdata,orig) inputSet where
  inputSet = inputGeneNames pdata

-- | Partition an SSD based on the value of a gene and sum the SSD (reducing
-- its size by 2 in the process).
splitAndSum :: (ParseData,SSD) -> Gene -> (ParseData,SSD)
splitAndSum (pd, sd) gene = (M.delete gene pd, U.zipWith (filtSum) sp1 sp2) where
  filtSum (_,x1) (_,x2) = x1+x2 -- Have to get rid of the indices
  (sp1,sp2) = U.partition evalGene (U.indexed sd) -- Wish there was a U.ipartition
  evalGene (ind,val) = testBit ind (numGenes-geneIndex-1)
  geneIndex = M.findIndex gene pd
  numGenes = M.size pd

-- | Find the names of the input genes
inputGeneNames :: ParseData -> [Gene]
inputGeneNames pdata = filter (inputGene pdata) (M.keys pdata)

-- | Use this to find input genes by the parse data.
inputGene :: ParseData -> Gene -> Bool
inputGene pdata g = null $ depends gi
  where
    Just gi = M.lookup g pdata

computeSSD :: IO ()
computeSSD = do
    let r file = liftM (U.fromList . head . read) $ readFile file :: IO (U.Vector Double)
    con1a <- r "moham_100.dat"
    con1b <- r "moham_111.dat"

    con2a <- r "moham_111.dat"
    con2b <- r "moham_101.dat"

    con3a <- r "moham_x1x.dat"
    con3b <- r "moham_1xx.dat"

    let ssd1 = U.sum . U.zipWith (\x y->(min x y)*0.5) con1a $ con1b
        ssd2 = U.sum . U.zipWith (\x y->(min x y)*0.5) con2a $ con2b
        ssd3 = U.sum . U.zipWith (\x y->(min x y)*0.5) con3a $ con3b

    print ssd1
    print ssd2
    print ssd3

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

    {-hypercorners = map (U.fromList.map fromIntegral.dec2bin n) [0..2^n-1]-}

    xs = clamped [0,0.5,1]
    center = replicate n 0.5
    clamped xs = map (\x-> if (x>1.0) then 1.0 else (if (x<0.0) then 0.0 else x)) xs
    grid = [[]] >>= foldr (>=>) return (replicate n (\x-> [x++[y] | y <- xs]))
    specialgrid = center : (grid \\ [center])
    hypergrid = map (U.fromList) specialgrid

    kmaps = map (fillKentries uvals ks) hypergrid
    ssds = parMap rdeepseq (collapseSSD p . simulateDOKUnif args . flip kmapToDOK 0) kmaps 

  writeFile "moham.dat" (pprint ssds)
  {-plotList [] (U.toList . head $ ssds)-}

  {-print $ map U.length ssds-}
  {-putStrLn ""-}
  print $ length ssds
  putStrLn ""
  {-print $ map (fst . normalizeSSD) ssds-}
  {-print $ "Optimum Theta: " ++ show optimumTheta-}

  return ()
  

class Pretty a where
  prettyPrint :: a -> String
  prettyPrint = renderStyle (Style LeftMode 80 0) . prettyShow
  prettyShow :: a -> Doc

pprint :: Pretty a => a -> String
pprint = prettyPrint

instance Pretty Int where
  prettyShow = text . show
instance Pretty Double where
  prettyShow = text . printf "%f"
instance Pretty a => Pretty [a] where
  prettyShow = brackets . hcat . punctuate comma . map prettyShow
instance (Pretty a, U.Unbox a) => Pretty (U.Vector a) where
  prettyShow = brackets . hcat . punctuate comma . map prettyShow . U.toList
instance (Pretty a, Pretty b) => Pretty (a,b) where
  prettyShow (a,b) = prettyShow a <> colon <> prettyShow b 

