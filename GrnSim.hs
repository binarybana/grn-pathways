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
        n1 = getRequiredArg args "n1"
        n2 = getRequiredArg args "n2"
        n3 = getRequiredArg args "n3"
        n4 = getRequiredArg args "n4"
        mode = getRequiredArg args "mode"
        gen = gotArg args "generate"
    if (length (["s","p","sm","d"]\\[mode]) == 4)
        then error $ usageError args "Choose a mode: s-state, p-pathway, sm-state modified, d-density estimation."
        else return ()
    con <- readFile inFile 
    let start = parsePW $ con ++ (unlines.words $ "tnf=0 ltbr=0 lps=0")
        p = parsePW $ con ++ (unlines.words $ getRequiredArg args "extra")

    when (mode == "s") $ do
        let initg = kmapToStateGraph p.buildKmaps $ p
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

    when (mode == "d") $ do
        let ks = buildKmaps p
            sim = (pass n1 (stripTransNodes.fullIteration)) $ kmapToStateGraph p ks
            initg = kmapToStateGraph p ks
            reduceGraph = gotArg args "reduce"
            firstGraph = if reduceGraph 
                then (pass n2 fullIteration).
                    (pass n1 (stripTransNodes.fullIteration)) $ initg
                else pass (n1+n2) fullIteration initg
            pass n f = last.take n.iterate f 

            inflist = iterate ((\(x,g)->(x,(pass 100 fullIteration g))).(reDrawEdges p ks)) (n4,firstGraph)
            infssas = map (genSSA.snd) inflist
            allgenes = map (\x-> map (!!x) (take n3 infssas)) [0..(length $ M.keys p)-1]
            allvars = map (stdDev.(V.fromList)) allgenes
            allests :: [[(Double, Double)]]
            allests = map ((\(pts,est)->G.toList (G.zip (fromPoints pts) est)).gaussianPDF n3.V.fromList) allgenes
            titleallests = zip titles allests
            titles = map (\x->defaultStyle {lineSpec = CustomStyle [LineTitle x]} ) (M.keys p)
            filtests = filter ((<200).maximum.map snd.snd) titleallests

            integrate xs = sum $ zipWith (*) (diff pts) (tail ys)
                where (pts,ys) = unzip xs

            diff (x:y:xs) = y-x:(diff (y:xs))
            diff _ = []

            correlate1 :: [Double] -> [Double] -> Double
            correlate1 ks = sum . zipWith (*) ks

            correlate :: [Double] -> [Double] -> [Double]
            correlate ks [] = []
            correlate ks xs = correlate1 ks xs : correlate ks (tail xs)

            convolve :: [Double] -> [Double] -> [Double]
            convolve ks = correlate (reverse ks)
            
        mapM_ (printf "%7s") (M.keys p)
        putStrLn ""
        mapM_ genPrintSSA (take n3 infssas)

        putStrLn ""
        putStrLn ""
        putStrLn "Std Deviations:"

        mapM_ (printf "%7s") (M.keys p)
        putStrLn ""
        mapM_ (printf "%7.4f") allvars

        putStrLn ""
        putStrLn ""
        putStrLn "Integrations:"

        mapM_ (printf "%7s") (M.keys p)
        putStrLn ""
        mapM_ (printf "%7.2f") (map ((\x->if x>1000 then 0 else x).integrate) allests)
        putStrLn ""

        --putStrLn ""
        --let (x,y) = unzip (head allests)
        --print $ (diff x)
        --putStrLn ""
        --print y

        plotPathsStyle [XRange (0.0,1.0)] filtests

        when gen $ drawStateGraph (snd.last.take n3 $ inflist) args
        return ()

    --when (mode == "d") $ do
    --    let ks = buildKmaps p
    --        initg = kmapToStateGraph p ks
    --        reduceGraph = gotArg args "reduce"
    --        firstGraph = if reduceGraph 
    --            then (pass n2 fullIteration).
    --                (pass n1 (stripTransNodes.fullIteration)) $ initg
    --            else pass (n1+n2) fullIteration initg
    --        pass n f = last.take n.iterate f 
    --        inflist = iterate ((pass 100 fullIteration).(reDrawEdges p ks)) firstGraph
    --        infssas = map genSSA inflist
    --    mapM_ (printf "%7s") (M.keys p)
    --    putStrLn ""
    --    mapM_ genPrintSSA (take 20 infssas)

    --    when gen $ drawStateGraph firstGraph args
    --    return ()

    when (mode == "sm") $ do
        let initg = kmapToStateGraph p.buildKmaps $ p
            reduceGraph = gotArg args "reduce"
            interGraph = (pass n2 fullIteration).
                    (pass n1 (stripTransNodes.fullIteration)) $ initg
            secondg = convertProbs interGraph (kmapToStateGraph p.buildKmaps $ p)
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

