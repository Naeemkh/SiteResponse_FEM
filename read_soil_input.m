function soil_pro=read_soil_input(rock_soil_type)
% Function for reading the soil properties.
% Input: soil types arrays
% Output: soil properties structure.

disp('-------> Loading soil properties ...');

n_soil_type= size(rock_soil_type,2);

for i=1:n_soil_type

   
    
    temp_ggmax = sprintf('%s%s%s','input_soil_pro/',rock_soil_type{i},'_ggmax.txt');
    temp_damping = sprintf('%s%s%s','input_soil_pro/',rock_soil_type{i},'_damping.txt');
    temp_prop = sprintf('%s%s%s','input_soil_pro/',rock_soil_type{i},'_prop.txt');
    
    
    t_ggmax   = load(temp_ggmax);
    t_damping = load(temp_damping);
    t_prop    = load(temp_prop);
    
    
    F1=sprintf('%s%s%s','soil_pro.s_r_type_',num2str(i),'.ggmax = t_ggmax;');
    F2=sprintf('%s%s%s','soil_pro.s_r_type_',num2str(i),'.damping = t_damping;');
    F3=sprintf('%s%s%s','soil_pro.s_r_type_',num2str(i),'.prop = t_prop;');
 
    eval(F1);
    eval(F2);
    eval(F3);

end

nn = sprintf('%s%2.0f','-------> Properties of',n_soil_type,' material(s) are loaded.');
disp(nn);
end