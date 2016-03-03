function element_index = meshing_domain(material_mat)

disp('-------> Meshing the simulation domain  ...');


num_material = size(material_mat,1);
end_meshing = 1;
i=1;
cur_material=1;


while end_meshing ~= 0
    
    
    for j=1:num_material
        
        if i==1
            cur_material=1;
        else
            
            if element_index(i-1,5) >= material_mat(j,4) && ...
                    element_index(i-1,5) <  material_mat(j,5)
                
                if abs(element_index(i-1,5) - material_mat(j,5))<0.0000001
                    if j==num_material
                        cur_material = j;
                    else
                        cur_material=j+1;
                        break
                    end
                else
                    cur_material=j;
                    break
                end
            end
        end
    end
    
        
    
    element_index(i,1) = i;    % element index
    element_index(i,2) = i;    % left node index
    element_index(i,3) = i+1;  % right node index
    element_index(i,6) = material_mat(cur_material,2); % density
    element_index(i,7) = material_mat(cur_material,3); % Module
    element_index(i,8) = material_mat(cur_material,7); % damping
    element_index(i,9) = material_mat(cur_material,9); % Material_type  based on material mat matrix
    
    
    if i==1
        space_inc = min(material_mat(cur_material,6),material_mat(cur_material,5));
    else
        space_inc = min(material_mat(cur_material,6),material_mat(cur_material,5)-element_index(i-1,5));
    end
    
    
    if (i==1) % node actual place.
        element_index(i,4) = 0;
        element_index(i,5) = element_index(i,4) + space_inc;
    else
        element_index(i,4) = element_index(i-1,5);
        element_index(i,5) = element_index(i,4) + space_inc;
    end
    
    
    
    if element_index(i,5) == material_mat(num_material,5)
        end_meshing = 0;
        element_index(i,:) = [];
        
    end
    
    
    i=i+1;
end

nn = sprintf('%s%s','-------> Number of total elements: ',num2str(size(element_index,1)));
disp(nn);


end