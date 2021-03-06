{-# LANGUAGE BangPatterns, RecordWildCards #-}
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
import GRN.Render

import qualified Data.Vector as V
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as U

import Data.Vector.Strategies
import Control.Parallel.Strategies

import System.Cmd
import Control.Monad
import Control.Applicative 
import Control.Arrow
import Data.List
import Data.Bits
import Data.Ord (comparing)
import Data.Maybe
import qualified Data.Map as M
import qualified Data.Set as S
import Text.Printf
import System.Environment
import System.Console.ParseArgs
import Data.Graph.Analysis 
import System.Random
import System.IO.Unsafe

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

--type ControlPolicy = ControlPolicy [(Gene,Bool)]
--data GeneInfo = GeneInfo {
--                    name        :: Gene,
--                    knockout    :: Maybe Bool,
--                    measurement :: Maybe Double,
--                    depends     :: [Gene],
--                    pathways    :: [Pathway] } deriving (Show,Eq)
--
--type ParseData = Map Gene GeneInfo
--data ParseControl = ParseControl { 
--                    pctargets :: [(Gene,Double)]
--                  , pccontrols :: [Gene] 
--                  } deriving (Show, Eq)
--
calcCost :: Args String -> ParseData -> ParseControl -> CSR -> Double
calcCost pargs pdata ParseControl{..} net = sum diffs where
  simSSD = simulateDOKUnif pargs (toDOK net)
  simSSA = ssdToSSA simSSD
  controlNames = map fst pctargets
  filtSSA = filter (\(gene,val) -> gene `elem` controlNames) (zip (M.keys pdata) (G.toList simSSA))
  diffs = zipWith (\(_,val1) (_,val2)-> (val1 - val2)**2) filtSSA pctargets

enumSinglePolicies :: [Gene] -> [ControlPolicy]
enumSinglePolicies conts = map (:[]) $ (,) <$> conts <*> [True,False]
--lets assume single gene control for the moment
-- to do all combinations, we would then want the power set of policies...
-- wow

policyExpansions :: Args String -> ParseData -> ParseControl -> [(ControlPolicy,DOK)]
policyExpansions pargs pdata pc@ParseControl{..} = zip policies dataFlows where
  policies = enumSinglePolicies pccontrols
  -- TODO FIXME will only work with single policies
  pdatas = map ((\(gene,val)-> M.adjust (\gi->gi{knockout=Just val}) gene pdata) . head) policies
  -- For each policy, generate a new parsedata with that gene knocked out 
  dataFlows = map (parseToDataFlow pargs) pdatas

policyCosts :: Args String -> ParseData -> ParseControl -> [(ControlPolicy,Double)]
policyCosts pargs pdata pc@ParseControl{..} = weightedCosts where
  expansesWithPols = map (id *** weightedExpansion) $ policyExpansions pargs pdata pc
  weight (cp, vec) = (cp, V.sum . V.map (\(val,net) -> val * calcCost pargs pdata pc net) $ vec)
  weightedCosts = map weight expansesWithPols

exactControlPolicy :: Args String -> ParseData -> ParseControl -> ControlPolicy
exactControlPolicy pargs pdata pc =  fst . minimumBy (comparing snd) . policyCosts pargs pdata $ pc

-- Generate 2*n weighted expansions (or samplings) for each of the n control
-- genes (for single gene stationary control), then calculate the costs for 
-- each of these, then take an n-way 'zip'
-- and determine control policy for each network (they should all be the same
-- right?) and the robust policy.
--
-- make sure and error on deterministic knockout policy, as only  a stochastic
-- knockout policy will preserve the total network
--
-- When sampling these won't be the same however... hmm what should we do then?
-- we'll worry about that when we get there. I think I can probably use the
-- same starting seed for each one and get the same sequence of networks
-- (again, as long as I am using complete hamming networks I should be fine).
-- but even using the same seed they will be different network draws.

uniqueSamples :: V.Vector (Double,CSR) -> V.Vector (Double,CSR)
uniqueSamples = V.fromList . S.toList . S.fromList . V.toList

printFact :: Show a => String -> a -> IO()
printFact str a = do
  putStr str
  print a
  putStrLn ""

printStats :: [ControlPolicy] -> [V.Vector (Double,CSR)] -> IO ()
printStats pols nets = do
  printFact "Number of policies: " (length pols)
  printFact "Number of networks: " (map V.length nets)
  printFact "Number of unique networks: " (map (V.length . uniqueSamples) nets)
  printFact "Fraction of Sampled Space: " (map (V.sum . fst . V.unzip . uniqueSamples) nets)
  
simControl :: Args String -> ParseData -> ParseControl -> IO ()
simControl pargs pdata pcon = do
        -- "Ghetto rigg!!!"
        let method = getRequiredArg pargs "submode"

        when (method == "e") $ do
          let (pols,likelynets) = unzip $ policyExpansions pargs pdata pcon
              expanded = map (weightedExpansion) likelynets
          printStats pols expanded

          return ()

        when (method=="s") $ do
          let n4 = getRequiredArg pargs "n4"
              (pols,likelynets) = unzip $ policyExpansions pargs pdata pcon
          if n4==0 then print "Error: You must have at least 1 sample with n4." else return ()

          samples <- mapM (getNSamples n4) likelynets
          printStats pols samples

          {-let gen = gotArg pargs "generate"-}
          {-let finalSSD = simulateDOKUnif pargs net-}
              {-graph = dataFlowToGraph net-}
          {-when gen $ drawDataFlow graph pargs-}
          return ()

            {-(parseToDataFlow pargs pdata)-}
            
            {-robustPol = robustPolicy controlPols-}

            {-finalSSD = simulateDOKUnif args initm-}
            {-graph = convertProbsVG finalSSD $ dataFlowToGraph initm-}
            {-(DOK (n,_) m) = initm-}
        {-mapM_ (printf "%7s") (M.keys p)-}
        {-putStrLn ""-}
        {-printSSA $ ssdToSSA finalSSD-}
        {-when gen $ drawDataFlow graph args-}

test :: Args String -> String -> IO DOK
test args x = do
  con <- fileNamePW x
  let dk = parseToDataFlow args con
  return dk
