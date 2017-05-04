function element_index = update_material_pro(output,eq_it)

element_index = output.element_index;
if eq_it == 1
    
    % just skip this function
    
else
    
    effective_strain = output.results.effective_strain *100;
    
    
    for ii=1:output.simulationparams.n_element
        
        if element_index(ii,12) == 0  % this element is not allowed to be updated.
            continue % the elemnt material won't be updated.
        end
        
        f1 = sprintf('%s%s%s','soil_details = output.simulationparams.soil_pro.s_r_type_',num2str(element_index(ii,9)),';');
        eval(f1)
        
        
        damping = interp1(soil_details.damping(:,1),soil_details.damping(:,2),effective_strain(ii,1));
        module_factor = interp1(soil_details.ggmax(:,1),soil_details.ggmax(:,2),effective_strain(ii,1));
        
        if isnan(damping)==1
            
            if effective_strain(ii,1) <= soil_details.damping(1,1)
                
                damping = soil_details.damping(1,2);
                
            elseif effective_strain(ii,1) >= soil_details.damping(end,1)
                
                damping = soil_details.damping(end,2);
                
            end
        end
        
        if isnan(module_factor)==1
            
            if effective_strain(ii,1) <= soil_details.ggmax(1,1)
                
                module_factor = soil_details.ggmax(1,2);
                
            elseif effective_strain(ii,1) >= soil_details.ggmax(end,1)
                
                module_factor = soil_details.ggmax(end,2);
                
            end
            
        end
        
        
        damping = damping;
        module  = module_factor * element_index(ii,11);
        
        element_index(ii,7)=module;
        element_index(ii,8)=damping; 
        
        
        
        
        
    end
    
end


end