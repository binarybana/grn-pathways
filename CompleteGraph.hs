import Data.Graph.Inductive
import Data.Graph.Inductive.Graph
import qualified Graphviz as Gv
import System.Cmd

data NodeInfo = NodeInfo String Int 

instance Show NodeInfo where
    show (NodeInfo name num) = name ++ (show num)

myGraph :: Int -> Gr NodeInfo String 
myGraph n = mkGraph (zip [1..n] (repeat $ NodeInfo "name" 1)) [(i,j,"") | i<-[1..n], j<-[1..n]]

test n = do
        --print $ labEdges $ myGraph n
        writeFile "complete.dot" (Gv.graphviz (myGraph n) "fgl" (5,5) (1,1) Gv.Portrait)
        rawSystem "dot" ["-Tsvg","-o","complete.svg","complete.dot"]
