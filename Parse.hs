module Parse where

import Text.Parsec
import Text.Parsec.String
import Control.Monad
import Data.Maybe
import qualified Data.Map as M

type Gene = String
data Pathway = Pathway [Gene] Bool Bool [Gene] deriving (Show)

data GeneInfo = GeneInfo {
                    name        :: Gene,
                    depends     :: [Gene],
                    pathways    :: [Pathway] } deriving (Show)

type ParseData = M.Map Gene GeneInfo

parsePass :: Parser (Maybe a) -> Parser [a]
parsePass fn =  liftM catMaybes $ sepEndBy line (char '\n')
                where line = try fn <|> ((skipMany (noneOf "\n")) >> return Nothing)

deps :: Parser (Maybe (Gene,[Gene]))
deps = do
        g <- many1 alphaNum
        char '('
        deps <- sepBy (many1 alphaNum) (char ',')
        char ')'
        return $ Just (g,deps)

paths :: Parser (Maybe Pathway)
paths = do
        gpre <- sepBy1 (many1 alphaNum) (string "&&")
        char '-'
        pre <- digit
        char ','
        post <- digit
        string "->"
        gpost <- sepBy1 (many1 alphaNum) (string "&&")
        let i2b x = if x == '0' then False else True
        return $ Just (Pathway gpre (i2b pre) (i2b post) gpost)

addPath (Pathway _ _  _  []) initMap =  initMap
addPath (Pathway a b1 b2 (p:ps)) initMap = addPath (Pathway a b1 b2 ps) postMap 
        where   postMap = M.adjust (\x -> x{pathways=pw:(pathways x)}) p initMap
                pw = Pathway a b1 b2 [p]

parsePW input = foldr addPath initMap pathList
        where   initMap = M.fromList $ map create depsList
                create (g,d) = (g, GeneInfo g d [])
                depsList = execute deps "Dependencies"
                pathList = execute paths "Pathways" 
                execute fn passname = case parse (parsePass fn) passname input of
                            Left err -> error ("Parsing error: " ++ (show err))
                            Right ret -> ret

fileNamePW :: String -> IO()
fileNamePW filename = do
        con <- readFile filename
        print $ parsePW con
