import os

from pySMOKEPostProcessor import elements_balance

# -------------------------------------------------------------------------------------
#  - - - Elements balance analysis - - -
# For an OpenSMOKE run, this returns one subplot per element, showing which species
# carries that element, and how much of it, across the whole run.
#
# All plotting utilities are built inside the post-processor function itself,
# calling elements_balance(...) does the parsing and shows the plots directly
#
# How to use:
#  kineticFolder    -- folder containing the mech (kinetics.xml is read), read the
#                      elements composition per species
#  resultsFolder    -- folder containing Output.xml (a single run), or a folder of
#                      Case*/Output.xml (a sweep, e.g. over temperature)
#  elements_list    -- elements symbols to plot, in a list of string as:["C", "O", "H"]
#                      any case is accepted, if not present in any species is skipped
#  threshold        -- *optional* puts a threshold at which a species is individually
#                      shown, under it is put in the "Others", base value is 5%
# -------------------------------------------------------------------------------------

kineticFolder = os.path.join("..", "data", "ROPA", "kinetics")
resultsFolder = os.path.join("..", "data", "ROPA", "Output")

elements_balance(
    kineticFolder,
    resultsFolder,
    elements_list=["C", "H", "O"],
    threshold=0.05,
)
