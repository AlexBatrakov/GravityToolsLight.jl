using GravityTools

eosname = :MPA1
path_to_grids = "/Users/abatrakov/Documents/PhD_work/projects/computed_grids"

grid = read_DEFGrid(eosname, path_to_grids)
mgrid = interpolate_DEFMassGrid(grid, -0.00123, -4.123)
interpolate_NS(mgrid, 1.5)


theory = DEF()
eosname = :MPA1
bnsys = BinarySystem(:NS, :WD)
sets = Settings("/Users/abatrakov/Documents/PhD_work/projects/computed_grids")
pf = DEFPhysicalFramework(theory, eosname, bnsys, sets)

pf.eosname = :MPA1
reinitialize!(pf)

pf.bnsys.PSR.mass = 1.5