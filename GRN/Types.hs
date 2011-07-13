{-# LANGUAGE BangPatterns #-}
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

import Data.Graph.Inductive
import Data.Map (Map)

type Gene = String
data Pathway = Pathway [Gene] Bool Bool [Gene] deriving (Show)

data GeneInfo = GeneInfo {
                    name        :: Gene,
                    knockout    :: Maybe Bool,
                    depends     :: [Gene],
                    pathways    :: [Pathway] } deriving (Show)

type ParseData = Map Gene GeneInfo

-----------------------

data NodeInfo = NodeInfo String !Double
data EdgeInfo = EdgeInfo Double Double
type ColoredStateGraph = Gr NodeInfo EdgeInfo
type DirectedGraph = Gr String String

type KmapSet = Map Gene Kmap 
data Kmap = Kmap Gene [Gene] (Map [Int] Kentry) 
    deriving Show
data Kentry = X | C Int | V Int Double | Const
    deriving Show 

instance Show NodeInfo where
    show (NodeInfo name num) = name ++ "\t" ++ (show num) ++ "\t"
instance Show EdgeInfo where
    show (EdgeInfo _ weight) = show weight

