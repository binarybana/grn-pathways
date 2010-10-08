module Plot (plot) where

import qualified Graphics.Gnuplot.Plot.TwoDimensional as Plot2D
import qualified Graphics.Gnuplot.Graph.TwoDimensional as Graph2D
import qualified Graphics.Gnuplot.Advanced as Advanced 
import qualified Graphics.Gnuplot.Terminal.X11 as X11

list2d :: [Int] -> Plot2D.T Int Int 
list2d vec =
   Plot2D.list Graph2D.listPoints vec

plot :: [Int] -> IO ()
plot vec = do
      Advanced.plot X11.cons (list2d vec)
      return ()
