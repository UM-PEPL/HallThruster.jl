# %%
import hallthruster as het
import matplotlib.pyplot as plt
import numpy as np

# %%

config = {
    "thruster": {
        "name": "SPT-100",
        "geometry": {
            "channel_length": 0.025,
            "inner_radius": 0.0345,
            "outer_radius": 0.05,
        },
        "magnetic_field": {
            "file": "bfield_spt100.csv"
        }
    },
    "propellants": [
        {
            "gas": "Xenon",
            "flow_rate_kg_s": 5e-6,
            "max_charge": 3,
        }
    ],
    "discharge_voltage": 300.0,
    "domain": (0.0, 0.08),
    "anom_model": {
        "type": "TwoZoneBohm",
        "c1": 0.00625,
        "c2": 0.0625,
    },
}

simulation = {
    "dt": 5e-9,
    "adaptive": True,
    "grid": {
        "type": "EvenGrid",
        "num_cells": 100,
    },
    "num_save": 100,
    "duration": 1e-3,
}

postprocess = {
    "output_file": "output.json",
    "save_time_resolved": False,
    "average_start_time": 5e-4,
}

# %% 

input = {"config": config, "simulation": simulation, "postprocess": postprocess}

solution = het.run_simulation(input, jl_env = "../")

# %%

f, axes = plt.subplots(1, 2, figsize = (10,5))
xlabel = 'Axial position [cm]'
axes[0].set_xlabel(xlabel)
axes[1].set_xlabel(xlabel)
axes[0].set_ylabel('\\$u_i\\$ [km/s]')
axes[1].set_ylabel('\\$T_e\\$ [eV]')
axes[0].set_title('Ion velocity')
axes[1].set_title('Electron temperature')

config = solution['input']['config']
propellant = config['propellants'][0]
avg = solution['output']['average']
z_cm = np.array(avg['z']) * 100

for j in range(propellant["max_charge"]):
    ui_km_s = np.array(avg["ions"]["Xe"][j]["u"]) / 1000
    axes[0].plot(z_cm, ui_km_s, label = f"{propellant['gas']} {j+1}+")

axes[1].plot(z_cm, avg['Tev'])
axes[0].legend()
plt.tight_layout()

plt.show()