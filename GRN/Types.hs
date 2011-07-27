{-# LANGUAGE BangPatterns,TypeSynonymInstances #-}
-- |
-- Module    : GRN.Types
-- Copyright : (c) 2011 Jason Knight
-- License   : BSD3
--
-- Maintainer  : jason@jasonknight.us
-- Stability   : experimental
-- Portability : portable
-- 
-- Handles all the type information
--

module GRN.Types where

--import Data.Graph.Inductive.PatriciaTree (Gr,UGr)
import Data.Graph.Inductive
import System.Console.ParseArgs
import Data.Map (Map)
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector as V
import Control.DeepSeq

type Gene = String
data Pathway = Pathway [Gene] Bool Bool [Gene] deriving (Show)

data GeneInfo = GeneInfo {
                    name        :: Gene,
                    knockout    :: Maybe Bool,
                    measurement :: Maybe Double,
                    depends     :: [Gene],
                    pathways    :: [Pathway] } deriving (Show)

type ParseData = Map Gene GeneInfo

data NodeInfo = NodeInfo String Double
data EdgeInfo = EdgeInfo Double Double
type ColoredStateGraph = Gr NodeInfo EdgeInfo
type DirectedGraph = Gr String String

instance NFData SSA where
    rnf vec = seq vec ()

-- An attractor contains the prior probability, the estimated densities for
-- each gene, and the normalized, reduced state graph.

type AttractorDensity = Map Gene (U.Vector Double)
type GeneMC = Map Gene (U.Vector Double)
type SSA = U.Vector Double
type SSD = U.Vector Double

type KmapSet = Map Gene Kmap 
data Kmap = Kmap Gene [Gene] (Map [Int] Kentry) 
    deriving Show
data Kentry = X | C Int | V Int Double | Const
    deriving Show 

instance Show NodeInfo where
    show (NodeInfo name num) = name ++ "\t" ++ (show num) ++ "\t"
instance Show EdgeInfo where
    show (EdgeInfo _ weight) = show weight

-- | A compressed row storage (CRS) sparse matrix.
--data CRS a = CRS {
--      crsValues :: Vector a
--    , colIndices :: Vector Int
--    , rowIndices :: Vector Int
--    } deriving (Show)

-- | A dictionary of values of (row,col) order 
data DOK = DOK (Int,Int) (Map (Int,Int) Double)
  deriving (Show)

-- | A compressed row storage (CRS) sparse matrix.
data CSC = CSC {
      matShape :: (Int,Int)
    , cscValues :: U.Vector Double
    , rowIndices :: U.Vector Int
    , colIndices :: U.Vector Int
    } deriving (Show)

data EMData = EMData {
      emAttractors :: Attractors
    , emUnfilledKMS :: KmapSet
    , emComps :: [ColoredStateGraph]
    , emCompNodes :: Masks
    , emNormSum :: Double
    , emArgs :: Args String
    , emMeasurements :: MeasureMap
    , emUncertainLocs :: [UncertainSpot] 
    } 

type Masks = [U.Vector Int]
type UncertainSpot = (Gene,[Int],Kentry)
type SSAMap = Map Gene Double
type MeasureMap = Map Gene Double
type Attractors = [(Double, SSAMap)]
type Theta = U.Vector Double

