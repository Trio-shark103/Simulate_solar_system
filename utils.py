import matplotlib.pyplot as plt
import numpy as np
import json


def read_json_config(file_path):
     """
          Read a JSON configuration file and return the configuration data.
          
          Parameters:
          -file_path(str): the path to the JSON configuration file
          
          
          Returns:
          -dict : a dictionary cotaining the configuration data.
          """

     # Load the JSON file for configuration
     with open(file_path, 'r') as file:
          config_data = json.load(file)
     return config_data


#Define the function that gets us the accn vector when passeed in the position vector
def accn(G, M_sun, r):
     """This is the formular that used to calculate the gravitational acceleration
     a= (-G * M_sun / |r|^3) * r_from_sun_to_earth

     where :
     -G is the universal Gravitational constant,
     -M_sun i s the mass of the sun ,
     |r| is the magnitude of the position vector r, and 
     -r_from_sun_to_earth is the unit vector pointing from the sun to earth.
     
     Note:
     -the negative sigm indicates that the acceleration is directed towards the sun"""
     return (-G*M_sun / np.linalg.norm(r)**3) * r


#Euler Integration 
def euler_method(G, M_sun, r,v,accn,dt):

     """Ordinary differnetial equation(ode)  for position
     dr/dt =v
     r_new = r_old  + V*dt

     # ode for velecity 
     dr/dt =a
     v_new = v_old  + a (r_old)*dt

     #Parameter
     #r: empty array for position of size t
     #v: empty array for velocity of size t
     #accn: func to calculate the accn at given position
     #dt: time step for the simulation
     This function will update the empty arrays for rand v with the simulated data"""

     for i in range(1, len(r)):
          r[i] = r[i-1] + v[i-1] * dt
          v[i] = v[i-1] + accn(G, M_sun, r[i-1]) *dt  
     
# RK4 Integration
def rk4_method(G, M_sun, r, v, accn, dt):
     """Ordinary differnetial equation(ode)  for position
     dr/dt =v
     r_new = r_old  + dt/6 (k1r + 2*k2r +2*k3r +k4r)

     # ode for velecity 
     dr/dt =a
     v_new = v_old  + dt/6(k1v +2*k2V +2*k3v +k4v)

     #Parameter
     #r: empty array for position of size t
     #v: empty array for velocity of size t
     #accn: func to calculate the accn at given position
     #dt: time step for the simulation
     This function will update the empty arrays for rand v with the simulated data
     Mrthod to calculate the steps 
     step1: - 0
     k1v = accn(r[i-1])
     k1r = v[i-1]

     step2: - dt/2 using step 1
     k2v = accn(r[i-1] + k1r *dt/2 )
     k2r =  v[i-1] + k1v * dt/2

     step3: - dt/2 using step 2
     k3v =accn(r[i-1] + k2r * dt/2)
     k3r = v[i-1] + k2v * dt/2
      
     step 4 :- dt using step 3
     k4v = accn (r[i-1] + k3r * dt)
     k4r = v[i-1] + k3v * dt

     this function will update the empthy arrays for r and v with the simulated data

     """

     for i in range (1 , len(r)): 
         # step1: - 0
          k1v = accn(G, M_sun, r[i-1])
          k1r = v[i-1]

          #step2: - dt/2 using step 1
          k2v = accn(G, M_sun, r[i-1] + k1r *dt/2 )
          k2r =  v[i-1] + k1v * dt/2

          #step3: - dt/2 using step 2
          k3v =accn(G, M_sun, [i-1] + k2r * dt/2)
          k3r = v[i-1] + k2v * dt/2
          
          #step 4 :- dt using step 3
          k4v = accn (G, M_sun, r[i-1] + k3r * dt)
          k4r = v[i-1] + k3v * dt

          # update the r and v
          v[i] = v[i-1]+ dt/6 * (k1v +2*k2v + 2*k3v + k4v )

          r[i] = r[i-1]+ dt/6 * (k1r +2*k2r + 2*k3r + k4r)


def numerical_integration (G, M_sun, r, v, accn, dt, method):
     """this function will apply the numerical)integratrion based on the method choosen 
     if the method is euler or rk4, the respective method will be implemented
     eles it will raise an error or exception
     Parameters
     #r: empty array for position of size t
     #v: empty array for velocity of size t
     #accn: func to calculate the accn at given position
     #dt: time step for the simulation"""

     if method.lower()== "euler":
          euler_method(G, M_sun,r,v,accn,dt)
     elif method.lower()=="rk4":
          rk4_method(G, M_sun,r,v,accn,dt)
     else:
          raise Exception (f'You can either choose "euler" or "rk4". Your current input for the method is:- {method}')

def at_aphelion(r, v):
     sizes = np.array([np.linalg.norm(position) for position in r])
     arg_aphelion = int(np.argmax(sizes))  # Ensure it's an int
     pos_aphelion = sizes[arg_aphelion]
     vel_aphelion = np.linalg.norm(v[arg_aphelion])

     return arg_aphelion, vel_aphelion, pos_aphelion


"""
     Determine the aphelion position and the associated parameters from simulate data
     
     parameters:
     -r (numpy.ndarray): simulated data array fpr the positions of the planet
     -v (numpy.ndarray): Simulated data array for the velocities of the planet.
     
     Returns:
     Tuple[int, float, float]L A tuple containing:
     -arg_aphelion(int): the ondex at which the planet is at its aphelion
     -vel_aphelion(float): The velocity of the planet at its aphelion
     -pos_aphelion(float): The position of the planet at its aphelion."""


# Find the point at which Earth is at its Aphelion ( point at ehich the earth is farthest away from the sun)


def plot_simulated_data(r, method_integration, arg_aphelion, vel_aphelion, pos_aphelion, name_planet, color_peri, color_ap):
     plt.style.use('dark_background')
     plt.figure(figsize=(7,12))
     plt.subplot(projection='3d')

     #Add Suptitle
     suptitle_str = 'RK4' if method_integration.lower()=='rk4' else 'Euler'
     plt.suptitle(suptitle_str + ' Method', color='r', fontsize =18, weight ='bold')

     # Add title
     title_str = f'At Aphelion, the {name_planet} is {round(pos_aphelion/1e9, 1)} million km away from the sun\nMoving at the speed at the speed of {round(vel_aphelion/1e3, 1)}km/s'
     plt.title(title_str, fontsize=14, color='orange')

     # Plot the orbits, sun, earth at perihelion and aphelion
     # plt.plot(r[:, 0], r[:, 1], color='tab:pink', lw=2, label='Orbit')
     plt.scatter(0, 0, color='yellow', s=1000, label='Sun')
     plt.scatter(r[0,0], r[0,1], s=200, label=f'{name_planet} at its Perihelion', color=color_peri)
     plt.scatter(r[arg_aphelion,0], r[arg_aphelion,1], s=200,
                    label=f'{name_planet} at its Aphelion', color=color_ap)

# Add Legend and Customize it
     legend = plt.legend(loc='lower right', frameon=False)
     for i in range(1, 4):
          legend.legend_handles[i]._sizes = [150] if i == 1 else [80]

# Turn off the axis and Display the result
     plt.axis('off') 
     plt.show()
