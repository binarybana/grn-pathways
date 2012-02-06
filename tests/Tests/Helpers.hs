-- | Helpers for testing
module Tests.Helpers where

import Data.Typeable

import qualified Test.HUnit      as HU
import Test.Framework
import Test.Framework.Providers.HUnit

import Numeric.MathFunctions.Constants (m_epsilon)
----------------------------------------------------------------
-- Helpers
----------------------------------------------------------------

data T a = T
-- | String representation of type name
typeName :: Typeable a => T a -> String
typeName = show . typeOf . typeParam
  where
    typeParam :: T a -> a
    typeParam _ = undefined

-- | Approximate equality for 'Double'. Doesn't work well for numbers
--   which are almost zero.
eq :: Double                    -- ^ Relative error
   -> Double -> Double -> Bool
eq eps a b 
  | a == 0 && b == 0 = True
  | otherwise        = abs (a - b) <= eps * max (abs a) (abs b)

-- | Approximately equal up to 1 ulp
(=~) :: Double -> Double -> Bool
(=~) = eq m_epsilon


----------------------------------------------------------------
-- HUnit helpers
----------------------------------------------------------------

testAssertion :: String -> Bool -> Test
testAssertion str cont = testCase str $ HU.assertBool str cont

testEquality :: (Show a, Eq a) => String -> a -> a -> Test
testEquality msg a b = testCase msg $ HU.assertEqual msg a b
