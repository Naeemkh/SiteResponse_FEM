function C1=boundary_condition(C)

% Absorbing boundary condition
size_k = size(C,1);

damping_value=C(size_k,size_k);

if abs(damping_value) < 0.1
   
    damping_value = 10000;
end


for ik=size_k-5:size_k-1
    
    C(ik,ik)    = 20 * damping_value;
    C(ik,ik+1)  = -10 * damping_value;
    C(ik+1,ik)  = -10 * damping_value;
    C(ik+1,ik+1)= 20 * damping_value;
    
end

C1=C;

end