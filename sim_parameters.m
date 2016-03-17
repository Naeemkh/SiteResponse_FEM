function sim_p = sim_parameters(dt,sim_time,sim_name,...
                       rock_soil_type,soil_pro,soil_layers,material_mat,...
                       use_damping,input_acceleration,acc_vec_1,max_value_acc,...
                       num_it)
                   
   
        sim_p.dt=dt;
        sim_p.sim_time = sim_time;
        sim_p.rock_soil_type = rock_soil_type;
        sim_p.soil_pro = soil_pro;
        sim_p.soil_layers = soil_layers;
        sim_p.material_mat = material_mat;
        sim_p.use_damping = use_damping;
        sim_p.input_acceleration=input_acceleration;
        sim_p.input_acc_vector = acc_vec_1;
        sim_p.sim_name = sim_name; 
        sim_p.num_it   = num_it;
        sim_p.max_value_acc = max_value_acc;
                   
end