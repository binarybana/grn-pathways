module GRN.Parse where

import Text.Parsec
import Text.Parsec.String
import Control.Monad
import Data.Maybe
import qualified Data.Map as M

type Gene = String
data Pathway = Pathway [Gene] Bool Bool [Gene] deriving (Show)

data GeneInfo = GeneInfo {
                    name        :: Gene,
                    knockout    :: Maybe Bool,
                    depends     :: [Gene],
                    pathways    :: [Pathway] } deriving (Show)

type ParseData = M.Map Gene GeneInfo

parsePass :: Parser (Maybe a) -> Parser [a]
parsePass fn =  liftM catMaybes $ sepEndBy line (char '\n')
                where line = try fn <|> ((skipMany (noneOf "\n")) >> return Nothing)

ss = skipMany space
tillEnd = skipMany (noneOf "\n")

deps :: Parser (Maybe (Gene,[Gene]))
deps = do
        ss
        g <- many1 alphaNum
        ss
        char '('
        ss
        deps <- sepBy (many1 alphaNum) (char ',')
        ss
        char ')'
        tillEnd
        return $ Just (g,deps)

paths :: Parser (Maybe Pathway)
paths = do
        let ampSep = try (ss >> string "&&" >> ss)
            geneID = (many1 $ noneOf "#-& \n") <?> "GeneID"
        gpre <- sepBy1 geneID ampSep
        ss
        char '-'
        pre <- digit
        char ','
        post <- digit
        string "->"
        ss
        gpost <- sepBy1 geneID ampSep
        tillEnd
        let i2b x = if x == '0' then False else True
        return $ Just (Pathway gpre (i2b pre) (i2b post) gpost)

knocks = do
        ss
        gene <- many1 alphaNum
        ss
        char '='
        ss
        kstate <- digit
        tillEnd
        let kbool = if kstate == '0' then False else True
        return $ Just (gene,kbool)

addKnock (gene, kbool) initMap = M.adjust (\x -> x{knockout=Just kbool}) gene initMap

addPath (Pathway _ _  _  []) initMap =  initMap
addPath (Pathway a b1 b2 (p:ps)) initMap = addPath (Pathway a b1 b2 ps) postMap 
        where   postMap = M.adjust (\x -> x{pathways=pw:(pathways x)}) p initMap
                pw = Pathway a b1 b2 [p]


parsePW input = foldr addKnock pathMap knockList
        where   pathMap = foldr addPath initMap pathList
                initMap = M.fromList $ map create depsList
                create (g,d) = (g, GeneInfo g Nothing d [])
                knockList = execute knocks "Knockouts"
                depsList = execute deps "Dependencies"
                pathList = execute paths "Pathways" 
                execute fn passname = case parse (parsePass fn) passname input of
                            Left err -> error ("Parsing error: " ++ (show err))
                            Right ret -> ret

fileNamePW :: String -> IO()
fileNamePW filename = do
        con <- readFile filename
        print $ parsePW con
