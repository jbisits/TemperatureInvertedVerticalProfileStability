using .VerticalProfileStability
using MAT

goship_files = glob("*.mat", GOSHIP_DATADIR)

goship_data = matopen(goship_files[1])
keys(goship_data)
close(goship_data)

vars = matread(goship_files[1])
keys(vars["D_reported"])

lons = vec(vars["D_reported"]["lonlist"][2])
lats = vec(vars["D_reported"]["latlist"][2])

scatter(lons, lats)

vars["D_reported"]["deplist"][1]
## Measurements I want
["CTDprs", "CTDCT", "CTDSA", "lonlist", "latlist", "deplist"]
