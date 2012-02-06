import Test.Framework       (defaultMainWithArgs)

import Tests.Parse

main :: IO ()
main = defaultMainWithArgs [ parseTests
                   ] ["--color"]
