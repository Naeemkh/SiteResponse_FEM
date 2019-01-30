%
% Scripts of running the simulaiton
%

%% loading the soil properties
soil_pro=read_soil_input(rock_soil_type); % read the input files.

%% Simulation params

% Simulation Name
serial_no = load('serial_no.txt');
serial_no = serial_no +1;
sim_name = sprintf('%s%s%s%s%s%s%s%s%s','sim_',num2str(serial_no),'_',num2str(use_damping),'_',num2str(sim_time),'_',num2str(dt),'_',num2str(max(max(soil_layers(:,3))))); % sim_usedamping_simtime_dt_elemntsize;
save('serial_no.txt','serial_no','-ascii');
nn = sprintf('%s%s','-------> Simulation Number: ',num2str(serial_no));
disp(nn);

if max_value_acc == -1
acc_vec_1 = load(input_acceleration);
else
acc_vec_1 = load(input_acceleration);
acc_vec_1(:,2)=(acc_vec_1(:,2)/abs(max(acc_vec_1(:,2))))*max_value_acc*g;
end

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
output.track_G          = zeros(size(element_index,1),num_it);
output.track_eff_strain = zeros(size(element_index,1),num_it);

for eq_it=1:num_it % Equivalent linear method's iteration
%% updating damping and stiffness of the elements based on the strain level.

element_index = update_material_pro(output,eq_it);

%% Generating Damping Matrix

[C,output] = damping_mat_gen(output,element_index);

%% Generate M and K matrix

% since in BKT2 I modify the shear modulus this part of code added here for
% now. I will polish it later.

element_index = output.element_index;

M_mat = mass_mat_gen(element_index);
K_mat = stiffness_mat_gen(element_index);




%% Boundary condition
% 
[C,M_mat,K_mat,element_index] = boundary_condition(C,M_mat,K_mat,element_index);

M_inv = inv(M_mat);
K = K_mat;


%% Time solution

output = solving_time(output,M_inv,M_mat,K,C,element_index,acc_vec_1,solution_type);

%% Report acceleration, velocity, and displacement and simulation params

output = waveform_gen(output,element_index,depth_results,t1,eq_it);

%% Extract strain

output=extract_strain_stress(output);

%% saving the value of G/Gmax and damping 

output = strain_damping_ggmax(output,eq_it);
output = transfer_function(output,eq_it);

disp('----------------------------------------------')
nn = sprintf('%s%s','-------> End of Iteration : ',num2str(eq_it),'.');
disp(nn);
disp('----------------------------------------------')


 
 output.track_eff_strain(:,eq_it) = output.results.effective_strain;
 output.track_G(:,eq_it) = output.element_index(:,7);


% for loop for iteration should end here.
end

%% saving the run summary for quick look
run_summary(output);

%% saving the whole simulation results
sname = sprintf('%s%s%s','simulation_results/simulation_',sim_name,'.mat');
save(sname,'output','-v7.3');

disp('-------> Done!')


