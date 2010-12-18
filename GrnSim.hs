module Main where 

import GRN.Parse
import GRN.Graphviz 
import GRN.StateTransition

import Data.Graph.Inductive hiding (Portrait)

import qualified Data.Map as M
import System.Cmd
import Control.Monad
import Text.Printf
import System.Environment

main = do
    args <- getArgs
    let n = head args
        nw = n ++ ".pw"
        nd = n ++ ".dot"
        ns = n ++ ".svg"
    con <- readFile nw
    let p = parsePW con
        initg = (initProbs.buildStateTGraph) p
        gl = take 15 $ iterate (stripTransNodes.fullIteration) initg
        gl2 = take 135 $ iterate fullIteration (last gl)
        g = last gl2
    mapM_ (printf "%7s") (M.keys p)
    putStrLn ""
    genPrintSSA $ genSSA g
    writeFile nd 
        (graphvizWithNodeFormatter cus g "fgl" (5,5) (1,1) Portrait)
    rawSystem "dot" ["-Tsvg","-o",ns,nd]


n = "nfkb-easy"
nw = n ++ ".pw"
nd = n ++ ".dot"
ns = n ++ ".svg"

testGraph = do 
        con <- readFile nw
        let p = parsePW con
        let g = (initProbs.buildStateTGraph) p
        let gl = take 30 $ iterate fullIteration g
        return $ last gl 

testPrint = do
    g <- testGraph 
    genPrintSSA $ genSSA g

testDraw = do
        con <- readFile nw
        let p = parsePW con
        print $ M.keys p
        let initg = (initProbs.buildStateTGraph) p
        let gl = take 30 $ iterate fullIteration initg
        let g = last gl
        writeFile nd 
            (graphvizWithNodeFormatter cus g "fgl" (5,5) (1,1) Portrait)
        --rawSystem "neato" ["-Tsvg","-o","nfkb.svg","nfkb.dot",
            --"-Gsplines=True","-Goverlap=false"]
        rawSystem "dot" ["-Tsvg","-o",ns,nd]
        return g

