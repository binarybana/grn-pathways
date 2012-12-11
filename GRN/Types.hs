{-# LANGUAGE BangPatterns,TypeSynonymInstances,FlexibleInstances #-}
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

import GRN.SparseLib
import Data.Graph.Inductive
import System.Console.ParseArgs
import Data.Map (Map)
import qualified Data.Vector.Unboxed as U
import Control.DeepSeq

type Gene = String
data Pathway = Pathway [Gene] Bool Bool [Gene] deriving (Show,Eq)

data GeneInfo = GeneInfo {
                    name        :: Gene,
                    knockout    :: Maybe Bool,
                    measurement :: Maybe Double,
                    depends     :: [Gene],
                    pathways    :: [Pathway] } deriving (Show,Eq)

type ParseData = Map Gene GeneInfo
data ParseControl = ParseControl { 
                    pctargets :: [(Gene,Double)]
                  , pccontrols :: [Gene] 
                  } deriving (Show, Eq)

data NodeInfo = NodeInfo String Double
data EdgeInfo = EdgeInfo Double Double
type ColoredStateGraph = Gr NodeInfo EdgeInfo
type DirectedGraph = Gr String String

-- | Right now we're assuming a stationary control policy
type ControlPolicy = [(Gene,Bool)]

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
    show (NodeInfo nam num) = nam ++ "\t" ++ (show num) ++ "\t"
instance Show EdgeInfo where
    show (EdgeInfo _ weight) = show weight

data EMData = EMData {
      emAttractors :: Attractors
    , emUnfilledKMS :: KmapSet
    ,nemComps :: [ColoredStateGraph]
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

