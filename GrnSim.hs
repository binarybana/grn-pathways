-- |
-- Module    : Main
-- Copyright : (c) 2011 Jason Knight
-- License   : BSD3
--
-- Maintainer  : jason@jasonknight.us
-- Stability   : experimental
-- Portability : portable
-- 
-- See the help information by executing the generated executable with no 
-- parameters.
-- 

module Main where 

import GRN.Parse
import GRN.StateTransition
import GRN.Render
import GRN.Types
import GRN.Density

import qualified Data.Map as M
import Data.List
import Control.Monad
import Text.Printf
import System.Environment
import System.Cmd
import System.Console.ParseArgs
import System.FilePath
import System.Directory

import Data.Graph.Analysis

import Data.Vector (Vector)
import qualified Data.Vector as V
import qualified Data.Vector.Generic as G
import Statistics.KernelDensity
import Statistics.Sample
import Graphics.Gnuplot.Simple

argList :: [Arg String]
argList = 
 [ Arg "output" (Just 'o') (Just "output") 
    (argDataDefaulted "FILE" ArgtypeString "out.svg") "Output file."
 , Arg "extra" (Just 'e') (Just "extra") 
    (argDataDefaulted "PW_ARGS" ArgtypeString "") 
    "Extra Pathway arguments; space delimited."
 , Arg "reduce" (Just 'r') (Just "reduce") Nothing 
    "Reduce output graph to attractor cycles." 
 , Arg "generate" (Just 'g') (Just "generate") Nothing
    "Generate an image"
 , Arg "open" (Just 'x') (Just "open") Nothing 
    "Open generated image automatically"
 , Arg "n1" Nothing (Just "n1") 
    (argDataDefaulted "INT" ArgtypeInt 15)
    "Default first elimination pass amount"
 , Arg "n2" Nothing (Just "n2")
    (argDataDefaulted "INT" ArgtypeInt 140)
    "Default second elimination pass amount"
 , Arg "n3" Nothing (Just "n3")
    (argDataDefaulted "INT" ArgtypeInt 140)
    "Runs for Density estimation ignored if m!=d"
 , Arg "n4" Nothing (Just "n4")
    (argDataDefaulted "INT" ArgtypeInt 0)
    ("If mode=d, then this is the seed, otherwise 0:deterministic equal valued "++
    "outgoing edge probabilities.")
 , Arg "mode" (Just 'm') (Just "mode")
    (argDataDefaulted "MODE" ArgtypeString "s")
    "s:state transition graph, p:pathway diagram, sm:modified stg, d:density estimation."
 , Arg "prog" (Just 'p') (Just "prog")
    (argDataDefaulted "PROG" ArgtypeString "dot")
    "Graphviz Draw Program to generate output"
 , Arg "input" Nothing Nothing
    (argDataRequired "FILE" ArgtypeString) "Input File."
 ]

main = do
    args <- parseArgsIO ArgsComplete argList
    let inFile = getRequiredArg args "input"
        n1 = getRequiredArg args "n1" :: Int
        n2 = getRequiredArg args "n2" :: Int
        mode = getRequiredArg args "mode"
        gen = gotArg args "generate"
    if (length (["s","p","sm","d"]\\[mode]) == 4)
        then error $ usageError args ("Choose a mode: s-state, p-pathway, "
            ++ "sm-state modified, d-density estimation.")
        else return ()
    con <- readFile inFile 
    let start = parsePW $ con ++ (unlines.words $ "tnf=0 ltbr=0 lps=0")
        p = parsePW $ con ++ (unlines.words $ getRequiredArg args "extra")

    when (mode == "s") $ do
        let initg = kmapToStateGraph p.buildKmaps $ p
            finalGraph = simulate args initg
        mapM_ (printf "%7s") (M.keys p)
        putStrLn ""
        printSSA $ genSSA finalGraph
        when gen $ drawStateGraph finalGraph args
        return ()

    when (mode == "d") $ do
        simRuns args p
        return ()

    when (mode == "sm") $ do
        let initg = kmapToStateGraph p.buildKmaps $ start
            secondg = convertProbs (simulate args initg) (kmapToStateGraph p.buildKmaps $ p)
            finalGraph = simulate args secondg
        mapM_ (printf "%7s") (M.keys p)
        putStrLn ""
        printSSA $ genSSA finalGraph
        when gen $ drawStateGraph finalGraph args
        return ()

    when (mode == "p") $ do
        let finalGraph = mkSimple $ buildGeneGraph p
        when gen $ drawGeneGraph finalGraph args
        return ()

