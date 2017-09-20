%
% Finite element code for 1D wave propagation boundary condition problem
% Written by: Naeem Khoshnevis (15 October 2015)
% Email: nkhshnvs@memphis.edu
% Last update: 26 April 2016

clc
clear all;
close all;

%% 
disp('                                   * * *                                ');
disp(' Finite element code for 1D wave propagation boundary condition problem ');
disp('     Analysis of site response through equivalent linear method         ');
disp('         Written by: Naeem Khoshnevis (15 October 2015)                 ');
disp('                         nkhshnvs@memphis.edu                           ');
disp('                     Last update: 9 February 2016                       ');
disp('                                   * * *                                ');


t1=tic;
%% Soil Parameters

% Soil layers
% Define type of soils and rock (for each type (e.g. soil) there should be
% soil_ggmax.txt, soil_damping.txt, and soil_prop.txt)
% soil_ggmax.txt ----> 2 columns c1: strain   c2: g/gmax
% soil_damping.txt ----> 2 columns c1: strain  c2: damping
% soil_prop.txt  ----> 4 columns c1: Vmax, c2: Gmax, c3: rho c4:damping 
% soil_prop.txt  ----> 4 columns c1: Vmax, c2: rho,  c3:damping 

rock_soil_type = {'rock','sand','sand2','sand3'};

%% Soil layers

% soil_layers --> 4 Columns : C1: soil type, C2: thick, C3: Max element
% size, C4: 1 = do equivalent linear process, 0= don't do eq linear process
% Bedrock should have two elements to be considered.


DRM_e_size = 8;
DRM_e_vs   = 640;
DRM_e_mu   = 1105920000;

soil_layers = [ 2 80                             8          1
                2 120                            8          1
                2 286                            8          1
                2 10                             8          0
                2 DRM_e_size*2                   DRM_e_size 0 
                 ];

depth_results=8:8:496; % Input depth that you want waveform for them.  

%% Simulation Parameters

sim_time      = 8;
dt            = 0.001;
input_acceleration = 'input_acc/Zeros.txt';
num_it        = 1;      % Number of iteration for equivalent linear method.
g             = 9.80665;
max_value_acc = 1; % coefficient for maximum value of the input as % of g.
solution_type = 'acc';  % acceleration (acc) will force the mass, displacement (disp) will dislocate the base node.
force_coeff   = (DRM_e_mu/DRM_e_size)*0.01; % temporal force coeffitient (mu/h)
timeshift     = (DRM_e_size/DRM_e_vs/dt); %temporal time shift (h/vs/dt)


% Damping options
% SRD   ==> Simplified Rayleigh (1 frequency Rayleigh Damping)
% FIRD  ==> Frequency Independent Rayleigh Damping (2 frequency Rayleigh
% Damping) %todo: misnomer, will be changed to RD2
% FDRD  ==> Frequency Dependent Rayleigh Damping 
% BKT   ==> Based on Bielak, Karaoglu and Taborda (2011) (Q from table)
% BKT2  ==> Based on Bielak, Karaoglu and Taborda (2011) (2 Maxwell elements)
% BKT3  ==> Based on Taborda, Huda, Khoshnevis and Bielak (2017) (3 Maxwell elements)
% BKT3F ==> Frequency dependent BKT3
% None  ==> Without damping model.
% Nonlinear ==> Will use nonlinear approach.

use_damping   = 'BKT2';


% * FIRD damping in this program is not frequency independent damping
% implemented in deepsoil. I will rename this damping.
%% Running the simulation

run_fem


