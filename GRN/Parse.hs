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
import Control.Applicative ((<$),empty,(<*))
import Control.Monad
import Data.Maybe
import qualified Data.Map as M


data ParseLine = ParseDependency (Gene, [Gene])
               | ParsePathway Pathway
               | ParseKnockout (Gene, Bool)
               | ParseMeasurement (Gene, Double)
                  deriving (Show,Eq)

dataLayers :: [ParseLine -> ParseData -> ParseData]
dataLayers = [addPath, addKnock, addMeasure]

parseLayers :: [Parser (Maybe ParseLine)]
parseLayers = [deps, paths, knocks, measurements]

parseLine :: Parser (Maybe ParseLine)
parseLine = do
    ss
    try (comment >> return Nothing) 
      <|> choice (map try parseLayers)

ss = skipMany space
tillNext = skipMany (noneOf "\n") >> ss
comment = ss >> char '#' >> tillNext
geneID = (many1 $ noneOf "#-&<>=: \n") <?> "GeneID"
p_float :: CharParser () Double
p_float = do  s<- getInput
              case readSigned readFloat s of
                [(n,s')] -> n <$ setInput s'
                _        -> empty

deps :: Parser (Maybe ParseLine)
deps = do
        g <- many1 alphaNum
        ss
        char '('
        ss
        deps <- sepBy (many1 alphaNum) (char ',')
        ss
        char ')'
        tillNext
        return $ Just $ ParseDependency (g, deps)

paths :: Parser (Maybe ParseLine)
paths = do
        let ampSep = try (ss >> string "&&" >> ss)
        gpre <- sepBy1 geneID ampSep
        ss
        char '-'
        pre <- digit
        char ','
        post <- digit
        string "->"
        ss
        gpost <- sepBy1 geneID ampSep
        tillNext
        let i2b x = if x == '0' then False else True
        return $ Just (ParsePathway $ Pathway gpre (i2b pre) (i2b post) gpost)

knocks :: Parser (Maybe ParseLine)
knocks = do
        gene <- many1 alphaNum
        ss
        char '='
        ss
        kstate <- digit
        tillNext
        let kbool = if kstate == '0' then False else True
        return $ Just $ ParseKnockout (gene, kbool)

measurements :: Parser (Maybe ParseLine)
measurements = do
        gene <- many1 alphaNum
        ss
        string "=measured="
        ss
        measure <- p_float
        tillNext
        return $ Just $ ParseMeasurement (gene, measure)


addKnock :: ParseLine -> ParseData -> ParseData
addKnock (ParseKnockout (gene, kbool)) initMap = M.adjust (\x -> x{knockout=Just kbool}) gene initMap
addKnock _ initMap = initMap

addMeasure :: ParseLine -> ParseData -> ParseData
addMeasure (ParseMeasurement (gene, measure)) initMap = M.adjust (\x -> x{measurement=Just measure}) gene initMap
addMeasure _ initMap = initMap

addPath :: ParseLine -> ParseData -> ParseData
addPath (ParsePathway (Pathway a b1 b2 (p:ps))) initMap = addPath (ParsePathway $ Pathway a b1 b2 ps) postMap 
        where   postMap = M.adjust (\x -> x{pathways=pw:(pathways x)}) p initMap
                pw = Pathway a b1 b2 [p]
addPath _ initMap = initMap

parseText :: Parser [ParseLine]
parseText = do ss
               plines <- many parseLine
               ss
               return (catMaybes plines)


parsePW :: String -> ParseData
parsePW input = foldr (\f dict -> foldr f dict plines) initMap dataLayers
        where   
                initMap = M.fromList $ map create depsList
                create (g,d) = (g, GeneInfo g Nothing Nothing d [])
                depsList      = [ x | ParseDependency x <- plines] 
                plines = case parse parseText "Input text" input of
                            Left err -> error ("Parsing error: " ++ (show err))
                            Right ret -> ret

fileNamePW :: FilePath -> IO ParseData
fileNamePW filename = do
        con <- readFile filename
        return $ parsePW con

parseControlText :: Parser ParseControl
parseControlText = do 
        ss
        manyTill anyChar (try (string "<control>"))
        ss
        targets <- between (string "<targets>" >> ss) (string "</targets>" >> ss) $ many1 $ do 
          g <- geneID
          ss
          char '='
          m <- p_float
          ss
          return (g,m)
        controls <- between (string "<controls>" >> ss) (string "</controls>" >> ss) $ many1 (geneID <* ss)
        string "</control>"
        ss
        return ParseControl { pctargets = targets, pccontrols = controls }

parseControl :: String -> ParseControl
parseControl input = pcontrol
    where   
        pcontrol = case parse parseControlText "Input Control text" input of
                    Left err -> error ("Parsing control error: " ++ (show err))
                    Right ret -> ret

