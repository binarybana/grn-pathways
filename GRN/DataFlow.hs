{-# LANGUAGE BangPatterns #-}
-- |
-- Module    : GRN.DataFlow
-- Copyright : (c) 2011 Jason Knight
-- License   : BSD3
--
-- Maintainer  : jason@jasonknight.us
-- Stability   : experimental
-- Portability : portable
-- 
-- Produces dataflow models from parsing info and then simulates these networks
-- if desired.
--
module GRN.DataFlow where

import GRN.Parse
import GRN.Types
import GRN.Utils
import GRN.Sparse

import qualified Data.Vector as V
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as U

import Data.Vector.Strategies
import Control.Parallel.Strategies

import System.Cmd
import Control.Monad
import Data.List
import Data.Bits
import Data.Ord (comparing)
import Data.Maybe
import qualified Data.Map as M
import Text.Printf
import System.Environment
import System.Console.ParseArgs
import Data.Graph.Analysis 
import System.Random
import System.IO.Unsafe


{-type ParseData = Map Gene GeneInfo-}
{-type Gene = String-}
{-data Pathway = Pathway [Gene] Bool Bool [Gene] deriving (Show,Eq)-}

--data GeneInfo = GeneInfo {
--                    name        :: Gene,
--                    knockout    :: Maybe Bool,
--                    measurement :: Maybe Double,
--                    depends     :: [Gene],
--                    pathways    :: [Pathway] } deriving (Show,Eq)
--

parseToDataFlow :: Args String -> ParseData -> DOK
parseToDataFlow args p = DOK (nStates,nStates) (M.fromList $ V.toList $ V.concatMap id pathPass)
    where
        genes = M.keys p
        nGenes = length genes
        nStates = 2^nGenes 
        startWeight = 1.0/(fromIntegral (nGenes+1))
        intStates = V.enumFromN 0 (nStates)
        startEdges = V.map genEdges intStates
        hammers = V.fromList $ (0 :) $ take nGenes $ map (2^) [0..]
        genEdges :: Int -> V.Vector ((Int,Int),Double)
        genEdges st = V.map (\h -> ((st,xor st h), startWeight)) hammers
        knockList = filter (isJust.snd) $ zip [0..] (reverse $ map knockout (M.elems p)) 
        pathList = concatMap pathways (M.elems p)
        knockfn = case getRequiredArg args "knock" of
                  "det" -> modKnockDet
                  "stoch"-> modKnockStoch gamma
        gamma = getRequiredArg args "gamma"
        knockPass = foldr (\k -> V.map (knockfn k)) startEdges knockList
        pathPass = foldr (\path -> V.map (modPath gamma genes path)) knockPass pathList

modKnockDet :: (Int,Maybe Bool) 
              -> V.Vector ((Int,Int),Double) 
              -> V.Vector ((Int,Int),Double)
modKnockDet (ind, Just val) group = normalize newgroup where
  newgroup = V.filter (\((_,st),_) -> val == testBit st ind) group

modKnockStoch :: Double 
                -> (Int,Maybe Bool) 
                -> V.Vector ((Int,Int),Double) 
                -> V.Vector ((Int,Int),Double)
modKnockStoch gamma (ind, Just val) group = normalize $ V.map modify group where
  modify c@((st1,st2),wt) 
    | val /= testBit st2 ind = ((st1,st2),wt*gamma)
    | otherwise = c

modPath :: Double -> [Gene] -> Pathway -> V.Vector ((Int,Int),Double) -> V.Vector ((Int,Int),Double)
modPath gamma glist (Pathway up pre post (affec:[])) group = normalize newgroup where
  newgroup = V.map applies group
  match val state ('!':gene) = val /= testBit state (getInd gene)
  match val state (gene) = val == testBit state (getInd gene)
  getInd g = fromJust $ elemIndex g (reverse glist) 
  -- ^^ reverse is for the right to left binary order, but l->r gene list
  applies c@((st1,st2),wt) 
    | all (match pre st1) up && not (match post st2 affec) = ((st1,st2),wt*gamma)
    | otherwise = c

normalize :: V.Vector ((Int,Int),Double) -> V.Vector ((Int,Int),Double)
normalize group = V.map (\((x,y),d) -> ((x,y),d*fac)) group where
  fac = (1.0 /) . V.sum . V.map snd $ group

dataFlowToGraph :: DOK -> ColoredStateGraph
dataFlowToGraph (DOK (n,m) mat) = mkGraph permedNodeStates edges where
  nGenes = floor $ logBase 2 (fromIntegral n)
  intStates = [0.. 2^nGenes-1]
  permedStates = map (dec2bin nGenes) intStates
  stringStates = map (concatMap show) permedStates
  permedNodeStates =
      zip intStates (map (\x->NodeInfo x (1.0/(fromIntegral n))) stringStates)
  edges = map (\((st1,st2),wt) -> (st1,st2,EdgeInfo 0.0 wt)) $ M.assocs mat

test :: Args String -> String -> IO DOK
test args x = do
  con <- fileNamePW x
  let dk = parseToDataFlow args con
  return dk

