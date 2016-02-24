function K_mat = stiffness_mat_gen(element_index)

% stiffness matrix
% K = K_mat * h * Mu (Module)
disp('-------> Generating the stiffness matrix ...');
n_e = size(element_index,1);
K_mat = zeros(n_e+1,n_e+1);

for i=1:n_e
    h  = element_index(i,5)-element_index(i,4); % defining h here will help the to define elements with different size.
    Mu = element_index(i,7);
    
    k_local(1,1)=1;
    k_local(1,2)=-1;
    k_local(2,1)=-1;
    k_local(2,2)=1;
    
    k_local_m=k_local*(1/h)*Mu;
    
    K_mat(i,i)     = K_mat(i,i)     +  k_local_m(1,1);
    K_mat(i,i+1)   = K_mat(i,i+1)   +  k_local_m(1,2);
    K_mat(i+1,i)   = K_mat(i+1,i)   +  k_local_m(2,1);
    K_mat(i+1,i+1) = K_mat(i+1,i+1) +  k_local_m(2,2);
    
end

end