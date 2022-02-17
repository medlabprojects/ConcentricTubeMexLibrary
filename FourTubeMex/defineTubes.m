function tubes = defineTubes()

tube1.OD = 1.1e-9; %outer diameter of inner tube [m]
tube1.ID = 1e-9; %inner diameter of inner tube [m]
tube1.k = 00; %curvature of the inner tube [m^-1]
tube1.L = 150e-3; %total length of the inner tube [m]
tube1.Lt = 100e-3; %length of inner tube straight section [m]
tube1.E = 50e9; %Young's Modulus of inner tube[Pa]
tube1.G = tube1.E/2/(1.4); %Shear Modulus of inner tube[Pa]

tube2.OD = 1.4e-3; %outer diameter of the middle tube [m]
tube2.ID = 1.3e-3; %inner diameter of the middle tube [m]
tube2.k = 30; %curvature of the middle tube [m^-1]
tube2.L = 100e-3; %total length of the middle tube [m]
tube2.Lt = 70e-3; %length of middle tube straight section [m]
tube2.E = 50e9; %Young's Modulus of the middle tube [Pa]
tube2.G = tube2.E/2/(1.4); %Shear modulus of the middle tube [Pa]

tube3.OD = 1.7e-3; %Outer diameter of the outer tube [m]
tube3.ID = 1.6e-3; %Inner diameter of the outer tube [m]
tube3.k = 0; %Curvature of the outer tube [m^-1]
tube3.L = 50e-3; %Total length of the outer tube [m]
tube3.Lt = 25e-3; %Length of outer tube straight section [m]
tube3.E = 50e9; %Young's modulus of the outer tube [Pa]
tube3.G = tube3.E/2/(1.4); %Shear modulus of the outer tube [Pa]

tube4.OD = 1.7e-3; %Outer diameter of the outer tube [m]
tube4.ID = 1.6e-3; %Inner diameter of the outer tube [m]
tube4.k = 0; %Curvature of the outer tube [m^-1]
tube4.L = 50e-3; %Total length of the outer tube [m]
tube4.Lt = 25e-3; %Length of outer tube straight section [m]
tube4.E = 50e9; %Young's modulus of the outer tube [Pa]
tube4.G = tube4.E/2/(1.4); %Shear modulus of the outer tube [Pa]


tubes.tube1 = tube1;
tubes.tube2 = tube2;
tubes.tube3 = tube3;
tubes.tube4 = tube4;
end
