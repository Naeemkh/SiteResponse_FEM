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

rock_soil_type = {'Bedrock','rock','clay','sand'};

%% Soil layers

%soil_layers --> 4 Columns : C1: soil type, C2: thick, C3: Max element
%size, C4: 1 = do equivalent linear process, 2= don't do eq linear process
% Bedrock should have two elements to be considered. 

soil_layers = [ 
                4 99       1    1
                4 1        0.1  1
                2 0.2      0.1  1
                                ];

depth_results=[1 5]; % Input depth that you want waveform for them.  

%% Simulation Parameters

sim_time      = 4;
dt            = 0.0001;
use_damping   = 2;      % 1-Simplified Rayleigh 2-Freq-Independent Rayleigh  3- Freq-dependent Rayleigh 4-BKT 5-None
input_acceleration = 'input_acc/ricker_10Hz.txt';
num_it        = 1;      % Number of iteration for equivalent linear method.
g             = 9.81;
max_value_acc = 0.1;    % coefficient for maximum value of the input as % of g.
solution_type = 'acc';  % acceleration (acc) will force the mass, displacement (disp) will dislocate the base node.

%% Running the simulation

run_fem


