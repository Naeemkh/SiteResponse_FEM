function M_mat = mass_mat_gen(element_index)

% generating the M matrix

% M = M_mat * h (element_size) * density;

disp('-------> Generating the mass matrix  ...');

n_e = size(element_index,1);
M_mat = zeros(n_e+1,n_e+1);

for i=1:n_e
    h  = element_index(i,5)-element_index(i,4); % defining h here will help the to define elements with different size.
    rho = element_index(i,6);
    
    m_local(1,1)=0.3333;
    m_local(1,2)=0.1666;
    m_local(2,1)=0.1666;
    m_local(2,2)=0.3333;
    
    
%     m_local(1,1)=0.5;
%     m_local(1,2)=0.0;
%     m_local(2,1)=0.0;
%     m_local(2,2)=0.5;
    
    
    
    m_local_m=m_local*(h)*rho;
    
    
    M_mat(i,i)     = M_mat(i,i)     + m_local_m(1,1);
    M_mat(i,i+1)   = M_mat(i,i+1)   + m_local_m(1,2);
    M_mat(i+1,i)   = M_mat(i+1,i)   + m_local_m(2,1);
    M_mat(i+1,i+1) = M_mat(i+1,i+1) + m_local_m(2,2);
        
end


end