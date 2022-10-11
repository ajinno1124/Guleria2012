cd ..
julia ./run/run_LY.jl
cd plot
python BindingEnergyL.py
python Plot_BE_LY.py
