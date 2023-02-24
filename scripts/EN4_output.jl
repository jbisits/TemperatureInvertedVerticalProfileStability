using .VerticalProfileStability

en4_data = readdir(EN4_DATADIR)[2]
en4_data_2021 = glob("*.nc", joinpath(EN4_DATADIR, en4_data))

ds = NCDataset(en4_data_2021; aggdim = "N_PROF")
lat = ds["LATITUDE"][:]
lon = ds["LONGITUDE"][:]
θ = ds["POTM_CORRECTED"][:, :]
Sₚ = ds["PSAL_CORRECTED"][:, :]
z = -ds["DEPH_CORRECTED"][:, :]
