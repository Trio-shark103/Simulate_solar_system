from utils import read_json_config, accn, numerical_integration, at_aphelion,plot_simulated_data
import numpy as np


def setup_simulation(config):
    planet_name = config['planet_info']['name']
    color_at_perihelion = config['planet_info']['perihelion_color']
    color_at_aphelion = config['planet_info']['aphelion_color']
    initial_position = np.array(config['initial_conditions']['position_at_perihelion']) * 1e9
    initial_velocity = np.array(config['initial_conditions']['velocity_at_perihelion']) * 1e3
    time_step = config['time_settings']['time_step']
    max_time = config['time_settings']['simulation_time'] * 24 * 3600
    method_integration = config['numerical_integration']['method']

    t = np.arange(0, max_time, time_step)
    r = np.empty((len(t), 2))
    v = np.empty((len(t), 2))

    # FIX: Proper initialization
    r[0] = initial_position
    v[0] = -initial_velocity  # reverse y-direction (if needed)

    return planet_name, color_at_perihelion, color_at_aphelion, r, v, t, time_step, method_integration


#read config.json
config= read_json_config("config.json")


# Constants 
G = 6.6743e-11
M_Sun = 1.989e30 #kg


# setup simulation
planet_name, color_at_perihelion, color_at_aphelion, \
    r,v,t,time_step,method_integration = setup_simulation(config)

# Call numerical integration
numerical_integration(G, M_Sun, r, v, accn, time_step, method=method_integration)

# Get data of earth at its aphelion 
arg_aphelion, vel_aphelion, pos_aphelion = at_aphelion(r, v) 


# Plot the simulated data
plot_simulated_data( r, method_integration, arg_aphelion,vel_aphelion,pos_aphelion,
                            planet_name, color_at_perihelion, color_at_aphelion)