# Aims of this project

The underlying question for this project is **does cabbeling shape the nature of vertical profiles in the gloabal ocean?**
To see if there is evidence for this happening we are using a highly idealised one-dimensional water column model, and global circulation model output from [ECCOv4r4](https://www.ecco-group.org/products-ECCO-V4r4.htm).
We may also look at observations.

Modules have been written for the one-dimensional water column model (`watercolmodel.jl`) and analysing the ECCO data (`maxdensitydiff.jl`).

## One dimensional Water column model

This model is used to investigate the stability of a water column to non-linear mixing processes between cold/fresh surface water and warm/salty deep water.
To do this we setup a two layer model and look at various initial conditions of salinity in the mixed (upper) layer to see if non-linear mixing processes are triggered at a certain threshold salinity value.

### Model setup

The model is two mdae up of two layers,

1. a mixed layer that is 100m in depth (0 - 100m) with initial conditions of constant temperature (set to $-1.85^{\circ}C$) and varying salinity values
2. a lower layer that is 400m in depth (100 - 500m) with initial conditions of constant temperature (set to $0.5^{\circ}C$) and increasing salinity.

The temperature of $-1.85^{\circ}C$ in the upper layer is chosen to be close to the freezing point of water as typically winter surface water in the marginal ice zone and other parts of the Southern Ocean will hover around this temperature.
We then set the salinity at the interface of the layers to be 34.7g/kg.
This gives us a warm/salty deep water mass sitting directly below a cold/fresh surface water mass.

The model is then set to use the TEOS-10 equation of state (EOS) to evolve the buoyancy in the model.
Note that using a non-linear EOS is critical for this experiments as we wish to observe non-linear mixing processes that are not present when a linear EOS is used.
To simulate convection that occurs if denser water is sitting on top of less dense water we use the convective adjusment diffusivity which applies a background (very small) diffusivity when the water column is gravitationallyu stable, that is if $\frac{\partial b}{\partial z} \geq 0$.
Where the water column is not gravtitationally stable ($\frac{\partial b}{\partial z} < 0$) a much stronger diffusivity is applied to simulate the effect of convection.

### Parameters

The only parameters that are set in the model are:

- the temperature of the mixed layer, $-1.85^{\circ}C$;
- the temperature of the lower layer, $0.5^{\circ}C$; and
- the salinity at the interface, 34.7 g/kg.

With these values we calculate everything else required.
We use a reference density for the equation of state defined by the density at the interface referenced to sea surface (as this is how we plot).
Using `GibbsSeaWater` we have `gsw_rho(34.7, 0.5, 100) ≈ 1028.18` which we take as the reference density for the model.
Other initial conditions for the salinity and temperature can be entered as `kwargs` (see below).

We then plot this reference density (though refernce to 0 dbar) in $T-S$ space and linearise about the water mass at the interface using
$$
T = T_{l} - \frac{\beta(S_{l}, T_{l}, 0)}{\alpha(S_{l}, T_{l}, 0)}\left(S - S_{l}\right)
$$
This tangent line to the isopycnal provides the critical condition for salinity.
We then choose values for the salinity either side of the tangent to set in the mixed layer.

The model is then run to see if the small amount of background mixing forms denser water trigger convection and if this is the case is it a cabbeling or static instability that causes this.

All parameters are set as defaults so the model (with the default parameters) is run by callin `run_model()`.
The default values for the parameters are

```julia
                    Tᵤ = -1.85,
                    Tₗ = 0.5,
                    Sᵤ = 34.5,
                    Sₗ = 34.7,
                    Sₘ = 34.75,
                    Sᵣ = 34.578,
               num_ics = 11,
                    Nz = 100, # number of levels
                    Lz = 500, # ovearll depth
           ref_density = 1028.18,
          convective_κ = 1.0,
          background_κ = 1e-5,
            convective = true,
                    Δt = 20 / 60,   # timestep in minutes, i.e. this is a 20 second timestep
              savepath = "Output",
            sim_length = 60,  # in days
             save_freq = 30   # in minutes
```

A subscript `ᵤ` indicates upper layer (i.e. the mixed layer) and a subscript `ₗ` is the lower layer.
The subscript `ₘ` is the maximum salinity for the lower layer and the subscript `ᵣ` is the maximum of the range of values for the salinity initial conditions in the mixed layer.

The simulation is setup to run a number of separate one dimensional water column models, each corresponding to an initial condition in the mixed layer.
The number of initial conditions is set by `num_ics` and has a default value of 11.
The salinity values used as initial conditions in the mixed layer are then derived by

```julia
sal_ics = range(Sᵤ, Sᵣ, length = num_ics)
```

Each model is saved to a **separate** `.jld2` file.
The module `DensityDiff` has functions to extract, reformat and save timeseries of the saved output.

**Note** due to the way this simulation is setup a single water column model cannot be run.

### Stability of the water column

To diagnose the static and cabbeling instabilities we find the maximum density difference between two levels of the water column that have a temperature difference greater than some chosen (by us) threshold $\Delta\Theta_{thres}$.  
The static density difference is defined as $\Delta \rho_{s} = \rho(S_{u}, T_{u}, \overline{p}) - \rho(S_{l}, T_{l}, \overline{p})$ (computed using [TEOS10](https://www.teos-10.org/pubs/gsw/html/gsw_contents.html)) where subscripts indicate upper and lower levels of the water column.
The cabbeling density difference is defined as $\Delta \rho_{c} = \rho_{0}(-\alpha(S_{l}, T_{l}, p_{l})(T_{u} - T_{l}) + \beta(S_{l}, T_{l}, p_{l})(S_{u} - S_{l}))$.
The maximum static and cabbeling density difference is then computed at each timestep so that we track the most unstable part of the water column.

## Two layer model

After runnning this highly resolved model an explicit two layer model was made.
See `run_model_two_layer()`, much of the model is the same (where possible) and the idea is to further experiment how the cabbeling instability could be parameterised in this two layer model.
