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

soil_pro=read_soil_input(rock_soil_type); % read the input files.

%% Soil layers

%soil_layers --> 4 Columns : C1: soil type, C2: thick, C3: Max element
%size, C4: 1 = do equivalent linear process, 2= don't do eq linear process

soil_layers = [ 
                4 100      1   1
                2 2        1   1
                               ];

depth_results=[1 5]; % Input depth that you want waveform for them.          
          
%% Simulation Parameters

sim_time      = 4;
dt            = 0.0001;
use_damping   = 1;      % 1-Simplified Rayleigh 2-Freq-Independent Rayleigh  3-BKT 4-None
input_acceleration = 'input_acc/ricker_10Hz.txt';
num_it        = 1;      % Number of iteration for equivalent linear method.
g             = 9.81;
max_value_acc = 0.1;    % coefficient for maximum value of the input as % of g.
solution_type = 'acc';  % acceleration (acc) will force the mass, displacement (disp) will dislocate the base node.


% Simulation Name
serial_no = load('serial_no.txt');
serial_no = serial_no +1;
sim_name = sprintf('%s%s%s%s%s%s%s%s%s','sim_',num2str(serial_no),'_',num2str(use_damping),'_',num2str(sim_time),'_',num2str(dt),'_',num2str(max(max(soil_layers(:,3))))); % sim_usedamping_simtime_dt_elemntsize;
save('serial_no.txt','serial_no','-ascii');
nn = sprintf('%s%s','-------> Simulation Number: ',num2str(serial_no));
disp(nn);

acc_vec_1 = load(input_acceleration);
acc_vec_1(:,2)=(acc_vec_1(:,2)/abs(max(acc_vec_1(:,2))))*max_value_acc*g;

%% Building Material Matrix

material_mat=build_material_mat(soil_pro,soil_layers);


%% Number of elements Meshing

element_index = meshing_domain(material_mat);

%% Initilizing the outcome matrix.

strmat=[];

%% Reporting simulation parameters

output = sim_parameters(dt,sim_time,sim_name,solution_type,...
                       rock_soil_type,soil_pro,soil_layers,material_mat,...
                       use_damping,input_acceleration,acc_vec_1,max_value_acc,...
                       num_it,...
                       element_index);
                   
                
%%                   

for eq_it=1:num_it % Equivalent linear method's iteration
%% updating damping and stiffness of the elements based on the strain level.

element_index = update_material_pro(output,eq_it);

%% Generate M and K matrix

M_mat = mass_mat_gen(element_index);
K_mat = stiffness_mat_gen(element_index);


%% Generating Damping Matrix

C = damping_mat_gen(M_mat,K_mat,material_mat,element_index,use_damping);

%% Boundary condition
% 
[C,M_mat,K_mat,element_index] = boundary_condition(C,M_mat,K_mat,element_index);

M_inv = inv(M_mat);
K = K_mat;


%% Time solution

output = solving_time(output,M_inv,M_mat,K,C,element_index,acc_vec_1,solution_type);

%% Report acceleration, velocity, and displacement and simulation params

output = waveform_gen(output,element_index,depth_results,t1);

%% Extract strain

output=extract_strain_stress(output);

%% saving the value of G/Gmax and damping 

output = strain_damping_ggmax(output,eq_it);

disp('----------------------------------------------')
nn = sprintf('%s%s','-------> End of Iteration : ',num2str(eq_it),'.');
disp(nn);
disp('----------------------------------------------')

% for loop for iteration should end here.
end

%% saving the run summary for quick look
run_summary(output);

%% saving the whole simulation results
sname = sprintf('%s%s%s','simulation_results/simulation_',sim_name,'.mat');
save(sname,'output','-v7.3');

disp('-------> Done!')


