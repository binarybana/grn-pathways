module Main where 

import GRN.Parse
import GRN.Graphviz 
import GRN.StateTransition

import Data.Graph.Inductive hiding (Portrait)
import Data.Graph.Analysis

import qualified Data.Map as M
import Data.List
import Control.Monad
import Text.Printf
import System.Environment
import System.Cmd
import System.Console.ParseArgs
import System.FilePath
import System.Directory

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
 , Arg "mode" (Just 'm') (Just "mode")
    (argDataDefaulted "MODE" ArgtypeString "s")
    "s:state transition graph, p:pathway diagram, sm:modified stg."
 , Arg "prog" (Just 'p') (Just "prog")
    (argDataDefaulted "PROG" ArgtypeString "dot")
    "Graphviz Draw Program to generate output"
 , Arg "input" Nothing Nothing
    (argDataRequired "FILE" ArgtypeString) "Input File."
 ]

main = do
    args <- parseArgsIO ArgsComplete argList
    let inFile = getRequiredArg args "input"
        n1 = getRequiredArg args "n1"
        n2 = getRequiredArg args "n2"
        mode = getRequiredArg args "mode"
        gen = gotArg args "generate"
    if (length (["s","p","sm"]\\[mode]) == 3)
        then error $ usageError args "Choose a mode: s-state, p-pathway, sm-state modified."
        else return ()
    con <- readFile inFile 
    let start = parsePW $ con ++ (unlines.words $ "tnf=0 ltbr=0 lps=0")
        p = parsePW $ con ++ (unlines.words $ getRequiredArg args "extra")

    when (mode == "s") $ do
        let initg = (initProbs.buildStateTGraph) p
            reduceGraph = gotArg args "reduce"
            finalGraph = if reduceGraph 
                then (pass n2 fullIteration).
                    (pass n1 (stripTransNodes.fullIteration)) $ initg
                else pass (n1+n2) fullIteration initg
            pass n f = last.take n.iterate f 
        mapM_ (printf "%7s") (M.keys p)
        putStrLn ""
        genPrintSSA $ genSSA finalGraph
        when gen $ drawStateGraph finalGraph args
        return ()

    when (mode == "sm") $ do
        let initg = (initProbs.buildStateTGraph) start
            reduceGraph = gotArg args "reduce"
            interGraph = (pass n2 fullIteration).
                    (pass n1 (stripTransNodes.fullIteration)) $ initg
            secondg = convertProbs interGraph (buildStateTGraph p)
            finalGraph = if reduceGraph 
                then (pass n2 fullIteration).
                    (pass n1 (stripTransNodes.fullIteration)) $ secondg
                else pass (n1+n2) fullIteration secondg
            pass n f = last.take n.iterate f 
        mapM_ (printf "%7s") (M.keys p)
        putStrLn ""
        genPrintSSA $ genSSA finalGraph
        when gen $ drawStateGraph finalGraph args
        return ()

    when (mode == "p") $ do
        let finalGraph = mkSimple $ buildGeneGraph p
        when gen $ drawGeneGraph finalGraph args
        return ()

drawGeneGraph gr args = do
    let outImg = getRequiredArg args "output"
        prog   = getRequiredArg args "prog"
        outDot = (dropExtension outImg) ++ ".dot"
    writeFile outDot
        (GRN.Graphviz.graphviz gr "fgl" (5,5) (1,1) Portrait)
    rawSystem prog ["-Tsvg","-o",outImg,outDot]
    when (gotArg args "open") $ do
        rawSystem "xdg-open" [outImg]
        return ()
    --removeFile outDot

drawStateGraph gr args = do
    let outImg = getRequiredArg args "output"
        prog   = getRequiredArg args "prog"
        outDot = (dropExtension outImg) ++ ".dot"
    writeFile outDot 
        (graphvizWithNodeFormatter cus gr "fgl" (5,5) (1,1) Portrait)
    --rawSystem "neato" ["-Tsvg","-o",outImg,outDot]
    rawSystem prog ["-Tsvg","-o",outImg,outDot]
    when (gotArg args "open") $ do
        rawSystem "xdg-open" [outImg]
        return ()
    --removeFile outDot
