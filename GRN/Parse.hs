-- |
-- Module    : GRN.Parse
-- Copyright : (c) 2011 Jason Knight
-- License   : BSD3
--
-- Maintainer  : jason@jasonknight.us
-- Stability   : experimental
-- Portability : portable
-- 
-- Parses a pathway file *.pw to produce the ParseData construction that can
-- then be processed further. See the example pathway files in pws folder 
-- for an example of the syntax.
--
module GRN.Parse where

import GRN.Types
import Text.Parsec
import Text.Parsec.String
import Text.ParserCombinators.Parsec.Char
import Text.Parsec.Token
import Numeric
import Control.Applicative ((<$),empty)
import Control.Monad
import Data.Maybe
import qualified Data.Map as M


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
        ss
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

knocks :: Parser (Maybe (Gene,Bool))
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

measurements :: Parser (Maybe (Gene,Double))
measurements = do
        ss
        gene <- many1 alphaNum
        ss
        string "=measured="
        ss
        measure <- p_float
        tillEnd
        return $ Just (gene,measure)

p_float :: CharParser () Double
p_float = do  s<- getInput
              case readSigned readFloat s of
                [(n,s')] -> n <$ setInput s'
                _        -> empty

validLine :: Parser (Maybe ())
validLine = do
            ss
            alphaNum <|> char '!'
            tillEnd
            return $ Just ()

addKnock (gene, kbool) initMap = M.adjust (\x -> x{knockout=Just kbool}) gene initMap

addMeasure :: (Gene, Double) -> ParseData -> ParseData
addMeasure (gene, measure) initMap = M.adjust (\x -> x{measurement=Just measure}) gene initMap

addPath (Pathway _ _  _  []) initMap =  initMap
addPath (Pathway a b1 b2 (p:ps)) initMap = addPath (Pathway a b1 b2 ps) postMap 
        where   postMap = M.adjust (\x -> x{pathways=pw:(pathways x)}) p initMap
                pw = Pathway a b1 b2 [p]


parsePW :: String -> ParseData
parsePW input = if not allUsed
                    then error "Syntax error in your pathway file."
                    else foldr addMeasure knockMap measureList
        where   knockMap = foldr addKnock pathMap knockList
                pathMap = foldr addPath initMap pathList
                initMap = M.fromList $ map create depsList
                create (g,d) = (g, GeneInfo g Nothing Nothing d [])
                knockList = execute knocks "Knockouts"
                measureList = execute measurements "Measurements"
                depsList = execute deps "Dependencies"
                pathList = execute paths "Pathways" 
                validList = execute validLine "Valids"
                allUsed = length validList == length pathList
                           + length depsList + length knockList + length measureList
                execute fn passname = case parse (parsePass fn) passname input of
                            Left err -> error ("Parsing error: " ++ (show err))
                            Right ret -> ret

fileNamePW :: FilePath -> IO ParseData
fileNamePW filename = do
        con <- readFile filename
        return $ parsePW con
