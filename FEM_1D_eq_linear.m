%
% Finite element code for 1D wave propagation boundary condition problem
% Written by: Naeem Khoshnevis (15 October 2015)
% Email: nkhshnvs@memphis.edu
% Last update: 9 February 2016

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

rock_soil_type = {'sand','rock','softsoil'};

soil_pro=read_soil_input(rock_soil_type); % read the input files.

%% Soil layers

%soil_layers --> 4 Columns : C1: soil type, C2: thick, C3: Max element
%size, C4: 1 = do equivalent linear process, 2= don't do eq linear process

soil_layers = [ 3 500   10   1
                1 500   10   1
                3 500   10   1
                2 500   10   1
              ];

depth_results=[10 50 ]; % Input depth that you want waveform for them.          
          
%% Simulation Parameters

sim_time      = 5;
dt            = 0.0001  ;
use_damping   = 2;      % 1-Simplified Rayleigh 2-Freq-Independent Rayleigh  3-BKT 4-None
input_acceleration = 'acc_mexican_hat_10hz.txt';


sim_name = sprintf('%s%s%s%s%s%s%s','sim_',num2str(use_damping),'_',num2str(sim_time),'_',num2str(dt),'_',num2str(max(max(soil_layers(:,3))))); % sim_usedamping_simtime_dt_elemntsize;

acc_vec_1 = load(input_acceleration);
acc_vec_1(:,2)=acc_vec_1(:,2)*20;

%% Building Material Matrix

material_mat=build_material_mat(soil_pro,soil_layers);
% num_material = size(material_mat,1);



%% Number of elements Meshing

element_index = meshing_domain(material_mat);



%% Generate M and K matrix



M_mat = mass_mat_gen(element_index);
K_mat = stiffness_mat_gen(element_index);

K = K_mat;
% K = K_mat *(le/n_e)* Mu;
M_inv = inv(M_mat);

% for ii=1:size(M_inv,1)
%     for jj=1:size(M_inv,1)
% 
%         if ii~=jj && abs(ii-jj)>1
%         
%             M_inv(ii,jj)=0;
%             
%         end
%     end
% end

%% Generating Damping Matrix

C = damping_mat_gen(M_mat,K_mat,material_mat,element_index,use_damping);

%% Boundary condition
C = boundary_condition(C);

%% Reporting simulation parameters

sim_p = sim_parameters(dt,sim_time,sim_name,...
                       rock_soil_type,soil_pro,soil_layers,material_mat,...
                       use_damping,input_acceleration,acc_vec_1);

%% Time solution


u = solving_time(sim_time,dt,M_inv,M_mat,K,C,element_index,acc_vec_1,'acc');


%% Report acceleration, velocity, and displacement and simulation params
output = waveform_gen(sim_p,u,element_index,depth_results,t1);

%% Extract strain
output=extract_strain_stress(output);





sname = sprintf('%s%s%s','simulation_results/simulation_',sim_name,'.mat');
save(sname,'output');




