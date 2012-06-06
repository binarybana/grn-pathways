GRN-Pathways
============
Welcome to the GRN-Pathways program. This is a program to produce biological 
genetic regulatory network models from biological pathway data. Simply create a 
pathway file (see the examples in the pws folder for syntax), and run the 
program with the appropriate options to generate pathway, state space, or 
conditional distribution graphs. Also enables setting various conditions and 
knockout conditions as well. 

NB: Must have gnuplot, the fftw3 library, and graphviz installed on the PATH in 
order for most of the commands to work. 

Installing
----------
Assuming you have Haskell installed and the appropriate prerequisites:

> cabal configure
> cabal build
> ./grnsim -m d -r -e "tnf=0 lps=1 ltbr=0" "pws/nfkb-super.pw" --n4 1

Or run ./simruns (a short python script) and compare to gold.txt to make sure 
everything is running alright (I know, its a poor man's way to test, but I've 
got paper deadlines to worry about!).

Future Work
-----------
You can find the paper published from some of this code in the [IEEE 
Transactions on Biomedical Engineering] [1].

I have moved on to other techniques due to the super-exponential behavior of 
some of the Hamming graph code in there (there is not much documentation, but 
if you have the slightest interest, please contact me!).

Problems
--------
Just let me know if you run into any problems or have any questions: 
(jason@jasonknight.us)

  [1]: http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=06176008 "IEEE"
