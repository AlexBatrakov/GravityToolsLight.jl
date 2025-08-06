using Revise
using GravityToolsLight
using Measurements
using PyPlot
using JLD
pygui(true)

eosname = :WFF1
path_to_grids = "/Users/abatrakov/Documents/Work/PhD/computed_grids_fine"

grid = read_DEFGrid(eosname, path_to_grids)
mgrid = interpolate_DEFMassGrid(grid, -0.00123, -4.123)
interpolate_NS(mgrid, 1.5)

