function [strain_matrix]=extract_strain_stress_nonlin(u,output)

element_index = output.element_index;

n_e = size(element_index,1);

strain_matrix=zeros(n_e,1);

element_size = element_index(:,5) - element_index(:,4);

d1 = u;

    d2=d1;
    d3=d1;
    d2(n_e+1)=[];
    d3(1)=[];
    
    s=(d3-d2)./element_size;
        
    strain_matrix(:,1)=s;

end