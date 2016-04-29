function output = sim_parameters(dt,sim_time,sim_name,solution_type,...
                       rock_soil_type,soil_pro,soil_layers,material_mat,...
                       use_damping,input_acceleration,acc_vec_1,max_value_acc,...
                       num_it,...
                       element_index)
                   
   
        output.simulationparams.dt=dt;
        output.simulationparams.sim_time = sim_time;
        output.simulationparams.n_timest = int64(ceil(sim_time/dt));
        output.simulationparams.rock_soil_type = rock_soil_type;
        output.simulationparams.soil_pro = soil_pro;
        output.simulationparams.soil_layers = soil_layers;
        output.simulationparams.material_mat = material_mat;
        output.simulationparams.damping.use_damping = use_damping;
        output.simulationparams.input_acceleration=input_acceleration;
        output.simulationparams.input_acc_vector = acc_vec_1;
        output.simulationparams.sim_name = sim_name; 
        output.simulationparams.num_it   = num_it;
        output.simulationparams.max_value_acc = max_value_acc;
        output.simulationparams.solution_type = solution_type;
        output.element_index = element_index;
                   
end