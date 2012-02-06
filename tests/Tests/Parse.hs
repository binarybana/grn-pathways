{-# OPTIONS_GHC -fno-warn-orphans #-}
{-# LANGUAGE ScopedTypeVariables  #-}
-- Required for Param
{-# LANGUAGE FlexibleInstances    #-}
{-# LANGUAGE OverlappingInstances #-}
module Tests.Parse (
    parseTests
  ) where

import Test.Framework                       (Test,testGroup)
import Test.Framework.Providers.QuickCheck2 (testProperty)
import Test.QuickCheck         as QC
import Test.QuickCheck.Monadic as QC
import Text.Printf

import GRN.Parse
import GRN.Types

import Prelude hiding (catch)
import Tests.Helpers

import qualified Data.Map as M

-- | Tests for all distributions
parseTests :: Test
parseTests = testGroup "Tests for Pathway Parsing"
  [ unitTests
  , testProperty "quickcheck test" revProp
  ]


revProp :: [Double] -> Bool
revProp x = x == (reverse . reverse $ x)

----------------------------------------------------------------
-- Unit tests
----------------------------------------------------------------

unitTests :: Test
unitTests = testGroup "Unit tests"
  [ testAssertion "testing" $ 1 =~ (1 + 1e-100)
  , testAssertion "reverse" $ (reverse [1,2,3]) == [3,2,1]
  , testAssertion "parsepw" $ (parsePW "d() \n a(d) \n d -1,0-> a\n") == 
    M.fromList [("d", GeneInfo "d" Nothing Nothing [] [])
               ,("a", GeneInfo "a" Nothing Nothing ["d"] [Pathway ["d"] True False ["a"]])]
  ]
