function output=extract_strain_stress(output)


node_displacement_time=output.nodetime;
element_index = output.element_index;

number_element  = output.simulationparams.n_element;
number_timest = output.simulationparams.n_timest;

strain_matrix=zeros(number_element,number_timest);

element_size = element_index(:,5) - element_index(:,4);

for i=1:number_timest

    d1=node_displacement_time(:,i);
    d2=d1;
    d3=d1;
    d2(number_element+1)=[];
    d3(1)=[];
    
    s=(d3-d2)./element_size;
    
    strain_matrix(:,i)=s;
    
     
end

effective_strain=max(strain_matrix,[],2);

output.results.strain_matrix=strain_matrix;
output.results.effective_strain = effective_strain;


end