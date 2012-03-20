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
import GRN.Sparse
import GRN.EM
import GRN.Uncertainty
import GRN.DataFlow

import qualified Data.Map as M
import Data.List
import Data.Ord
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
import Statistics.Sample.KernelDensity
import Statistics.Sample
import Graphics.Gnuplot.Simple

argList :: [Arg String]
argList = 
 [ Arg "output" (Just 'o') (Just "output") 
    (argDataDefaulted "FILE" ArgtypeString "out.svg") "Output file."
 , Arg "extra" (Just 'e') (Just "extra") 
    (argDataDefaulted "PW_ARGS" ArgtypeString "") 
    "Extra Pathway arguments; space delimited."
 , Arg "extra2" (Just 'f') (Just "extra2") 
    (argDataDefaulted "PW_ARGS" ArgtypeString "") 
    "Extra Pathway arguments; space delimited For SM runs."
 , Arg "reduce" (Just 'r') (Just "reduce") Nothing 
    "Reduce output graph to attractor cycles." 
 , Arg "generate" (Just 'g') (Just "generate") Nothing
    "Generate an image"
 , Arg "open" (Just 'x') (Just "open") Nothing 
    "Open generated image automatically"
 , Arg "avg" Nothing (Just "avg") 
    (argDataDefaulted "INT" ArgtypeInt 8)
    "Averaging sums to perform"
 , Arg "n1" Nothing (Just "n1") 
    (argDataDefaulted "INT" ArgtypeInt 15)
    "First elimination (reduction) pass amount"
 , Arg "n2" Nothing (Just "n2")
    (argDataDefaulted "INT" ArgtypeInt 140)
    "Second pass (only simulation) amount"
 , Arg "n3" Nothing (Just "n3")
    (argDataDefaulted "INT" ArgtypeInt 140)
    "Runs for Density estimation ignored if m!=d"
 , Arg "n4" Nothing (Just "n4")
    (argDataDefaulted "INT" ArgtypeInt 0)
    ("If mode=d, then this is the seed, otherwise 0:deterministic equal valued "++
    "outgoing edge probabilities.")
 , Arg "mode" (Just 'm') (Just "mode")
    (argDataRequired "MODE" ArgtypeString)
    ("s:state transition graph, ss:state transition matrix, p:pathway diagram, \
    \sm:modified stg, d:density estimation, em:expectation maximization, \
    \m:Mohammad output, df: DataFlow, c: Control")
 , Arg "dmode" Nothing (Just "dmode")
    (argDataDefaulted "DMODE" ArgtypeString "m")
    "m: Matrix multiplication, g: graph"
 , Arg "prog" (Nothing) (Just "prog")
    (argDataDefaulted "PROG" ArgtypeString "dot")
    "Graphviz Draw Program to generate output"
 , Arg "knock" (Just 'k') (Just "knock")
    (argDataDefaulted "KNOCK" ArgtypeString "stoch")
    "Knockouts in Dataflow mode: det(erministic), stoch(astic)"
 , Arg "gamma" Nothing (Just "gamma") 
    (argDataDefaulted "DOUBLE" ArgtypeDouble 0.1)
    "Amount to modify edges in dataflow mode."
 , Arg "perturb" (Just 'p') (Just "perturb") 
    (argDataDefaulted "DOUBLE" ArgtypeDouble 0.0)
    "In Mohammad output, set the perturbation probability."
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
        modes = ["s","p","sm","d","ss","em","m","df","c"]
    if (length (modes\\[mode]) == length mode)
        then error $ usageError args ("Choose a valid mode.")
        else return ()
    con <- readFile inFile 
    let start = parsePW $ (unlines.words $ getRequiredArg args "extra2") ++ con
        p = parsePW $ (unlines.words $ getRequiredArg args "extra") ++ con 

    when (mode == "s") $ do
        let initg = kmapToStateGraph.buildKmaps $ p
            finalGraph = simulate args initg
        mapM_ (printf "%7s") (M.keys p)
        putStrLn ""
        printSSA $ genSSA finalGraph
        when gen $ drawStateGraph finalGraph args
        return ()

    when (mode == "ss") $ do
        let initm = (kmapToDOK.buildKmaps $ p)  0
            finalSSD = simulateDOKUnif args initm
            (DOK (n,_) m) = initm
        mapM_ (printf "%7s") (M.keys p)
        putStrLn ""
        --print $ maximum $ map (snd) $ M.keys m
        --print $ G.length $ colIndices $ dokToCSC initm
        --print $ G.length $ rowIndices $ dokToCSC initm
        --print $ G.length $ cscValues $ dokToCSC initm
        printSSA $ ssdToSSA finalSSD
        --Cant draw a graph when we have a matrix... well we could, but we'd
        --have to convert:
        --when gen $ drawStateGraph finalGraph args
        return ()

    when (mode == "df") $ do
        let initm = parseToDataFlow args p
            finalSSD = simulateDOKUnif args initm
            graph = convertProbsVG finalSSD $ dataFlowToGraph initm
            (DOK (n,_) m) = initm
        mapM_ (printf "%7s") (M.keys p)
        putStrLn ""
        printSSA $ ssdToSSA finalSSD
        when gen $ drawDataFlow graph args
        return ()

    when (mode == "c") $ do
        let pcon = parseControl con
        simControl args p pcon

    when (mode == "d") $ do
        simRuns args p
        return ()

    when (mode == "m") $ do
        uncertaintyPrint args p
        return ()

    when (mode == "sm") $ do
        let initm = (kmapToDOK.buildKmaps $ start)  0
            secm = (kmapToDOK.buildKmaps $ p) 0
            interSSD = simulateDOKUnif args initm
            (DOK (n,_) m) = initm
            finalSSD = simulateDOK args secm interSSD
        mapM_ (printf "%7s") (M.keys p)
        putStrLn ""
        printSSA $ ssdToSSA finalSSD
        return ()

    when (mode == "p") $ do
        let finalGraph = mkSimple $ buildGeneGraph p
        when gen $ drawGeneGraph finalGraph args
        return ()

defArgs :: Args String
defArgs = parseArgs ArgsComplete argList "./grnsim" ["-m","df","pws/doesntmatter.pw"]
