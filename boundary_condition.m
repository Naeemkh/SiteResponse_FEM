function C1=boundary_condition(C,M_mat,element_index)

% Absorbing boundary condition
v_s = sqrt(element_index(end,7)/element_index(end,6));
C(end,end) = M_mat(end,end)*v_s*0.6;
C1=C;
end