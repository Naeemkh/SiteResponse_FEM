function material_mat=build_material_mat(soil_pro,soil_layers)

% material_mat = [No, density, module, x1, x2, maximum element, damping,
% eq_process]
disp('-------> Building material matrix ...');

n_layers=size(soil_layers,1);
material_mat=zeros(1,7);

for i=1:1:n_layers
    if i == 1
        
        
        % extracting the properties of the soil layer.
        
        f1=sprintf('%s%s%s','soil_property = soil_pro.s_r_type_',num2str(soil_layers(i,1)),'.prop;');
        eval(f1);
        
        material_mat(1,1) = 1;
        material_mat(1,2) = soil_property(1,2);
        material_mat(1,3) = soil_property(1,1).^2*soil_property(1,2);
        material_mat(1,4) = 0;
        material_mat(1,5) = soil_layers(1,2);
        material_mat(1,6) = soil_layers(1,3);
        material_mat(1,7) = soil_layers(1,4);
        
    elseif i > 1
        
        % extracting the properties of the soil layer.
        
        
        f1=sprintf('%s%s%s','soil_property = soil_pro.s_r_type_',num2str(soil_layers(i,1)),'.prop;');
        eval(f1);
        
        material_mat(i,1) = i;
        material_mat(i,2) = soil_property(1,2);
        material_mat(i,3) = soil_property(1,1).^2*soil_property(1,2);
        material_mat(i,4) = material_mat(i-1,5);
        material_mat(i,5) = material_mat(i,4)+soil_layers(i,2);
        material_mat(i,6) = soil_layers(i,3);
        material_mat(i,7) = soil_layers(i,4);
        
       
    end
    
end



