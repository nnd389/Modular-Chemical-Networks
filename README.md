
# Modular-Chemical-Networks
This has code that implements a variety of chemical networks using Julia!



In particular, there are three main networks implemented here. A note on chemical network ODEs: Writing up the system of ODEs from a chemical network is a simple but tedious task. You can find tutorials on how to do this online, or see the general form in my poster. (See Cosmic AI Horizons Conference.pdf) All the networks were implemented using Julia's Catalyst.jl, which does the work of converting the reactions into ODE form for you. But if you'd like to see the ODE form of a network, see https://docs.sciml.ai/SciMLBenchmarksOutput/stable/AstroChem/nelson/. You can find the script for that page in /scripts/Nelson/Nelson_network_fordocs.jl. 

## Small: Nelson and Langer (1999)
This network has 23 reactions, 14 chemical species and comes from the paper "On the Stability and Evolution of Isolated BOK Globules" by Nelson and Langer written in 1999. 
See /scripts/Nelson/Nelson_network_fordocs.jl for the network in ODE form and /scripts/Nelson/Nelson_catalyst.jl for the network in Catalyst (reaction) form. 

## Medium: Glover et. al. (2010)
This network has 219 reactions and 32 chemical species and comes from the paper "Modelling CO formation in the turbulent interstellar medium" by Glover, Federrath, Mac Low, and Klessen in 2010. 
See /scripts/Glover/Glover.jl

## Large: The UMIST Database for Astrochemistry
At the time of implementing it, this network had 6173 reactions and 468 chemical species. It comes from the paper "The UMIST database for astrochemistry 2022" by Millar, Walsh, Van de Sande, and Markwick in 2022. They quite recently updated the network to have around 8000 reactions, and it just so happens that at the time of implementing it, I narrowly missed the update. So the version you will find here is an older version of UMIST. Millar and others if you are reading this, thank you for putting such great effort to create a great network. You can imagine my surprise when I found out I barely missed the update. Thankfully, I didn't do it all by hand. See /scripts/UMIST/Umist_xlsx_to_julia.py for a script that takes in the UMIST database in excel form and converts it to Julia code. So jokes on you, I came prepared! (I still only implemented the 6000 reaction version though hehe)

### Some lamps to light your way
- At the top of most of the julia or python scripts you might find a "Note from Nina". Please refer to that note to get a brief idea of what each script is for. 
- Timing: Please be careful using the Julia @time macro, nust know that when you use it, the first time estimate is moot, then the others will work. Meaning for each function in your julia script that you want to time, you have to put the time macro at least twice and only refer to the second estimate.
- Many of the scripts were used to create data that is used in my most recent poster, meaning all networks were solved with two  different  sets  of  initial conditions (u0). First, “Fiducial u0” with initial values H2 = 0.5, He = 0.1, C+ = 1e-9, C = 5e-9, O = 1e-8, HCO+ = 9e-9, CO = 9.9e-5, H2O = 8e-5, O2 = 1e-8, H3+ = 1e-8, and e = 0. Second, “Neutral u0” where initial abundances were changed  but  atom  amounts were conserved with values H2 = 0.5, He = 0.1, C+ = 0.0,  C = 9.91e-5, O = 1.79e-4, HCO+ = 0.0, CO = 0.0, H2O = 8e-5, O2 = 1e-8, H3+ = 0.0, and e = 2e-8. 
- I tried to replace the word "ionized" with "fiducial" but might have missed a few, if you see an "ionized" something I mean fiducial.
- 
