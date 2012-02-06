{-# LANGUAGE BangPatterns, OverloadedStrings #-}
-- |
-- Module    : GRN.Render
-- Copyright : (c) 2011 Jason Knight
-- License   : BSD3
--
-- Maintainer  : jason@jasonknight.us
-- Stability   : experimental
-- Portability : portable
-- 
-- Renders StateGraphs, GeneGraphs and draws Gnuplot graphs
--

module GRN.Render where

import Data.Graph.Inductive
import Data.Graph.Inductive.Graph
import Data.GraphViz
import Data.GraphViz.Attributes
import Data.GraphViz.Attributes.Complete

import GRN.Types
import GRN.Parse
import GRN.StateTransition

import System.Environment
import System.Cmd
import Text.Printf
import System.Console.ParseArgs
import System.FilePath
import System.Directory
import Control.Monad
import Data.Text.Lazy (pack)
import qualified Data.Text.Lazy.IO as T


-- data NodeInfo = NodeInfo String !Double
-- data EdgeInfo = EdgeInfo Double Double
-- type ColoredStateGraph = Gr NodeInfo EdgeInfo
-- type DirectedGraph = Gr String String

labelFn :: (Node,NodeInfo) -> Attributes
labelFn (_, (NodeInfo name prob)) = [Label . StrLabel . pack $ name, FillColor $ HSV 0 0 shade, 
        Style [SItem Filled []]]
  where   shade = 1 - (clamp $ prob**0.4) -- The 0.4 is a gamma correction factor
          clamp x = if x>1 then 1 else if x<0 then 0 else x

edgeLabel :: (Node, Node, EdgeInfo) -> Attributes
edgeLabel (_,_,(EdgeInfo _ w)) = [Label . StrLabel . pack $ (printf "%3.1f" w)]

drawGeneGraph :: DirectedGraph -> Args String -> IO ()
drawGeneGraph gr args = do
    let outImg = getRequiredArg args "output"
        prog   = getRequiredArg args "prog"
        outDot = (dropExtension outImg) ++ ".dot"
        dotCons = graphToDot nonClusteredParams gr
        con = printDotGraph dotCons
    T.writeFile outDot con
    rawSystem prog ["-Tsvg","-o",outImg,outDot]
    when (gotArg args "open") $ do
        rawSystem "xdg-open" [outImg]
        return ()
    --removeFile outDot

drawStateGraph :: ColoredStateGraph -> Args String -> IO ()
drawStateGraph gr args = do
    let outImg = getRequiredArg args "output"
        prog   = getRequiredArg args "prog"
        outDot = (dropExtension outImg) ++ ".dot"
        dotCons = graphToDot nonClusteredParams{ fmtNode = labelFn } gr
        con = printDotGraph dotCons
    T.writeFile outDot con
    rawSystem prog ["-Tsvg","-o",outImg,outDot]
    when (gotArg args "open") $ do
        rawSystem "xdg-open" [outImg]
        return ()
    --removeFile outDot

