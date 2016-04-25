function [C1,M_mat1,K_mat1,element_index1]=boundary_condition(C,M_mat,K_mat,element_index)

% Absorbing boundary condition
% 
% C(:,end)=[];
% C(end,:)=[];
% 
% M_mat(:,end)=[];
% M_mat(end,:)=[];
% 
% K_mat(:,end)=[];
% K_mat(end,:)=[];
% 
% element_index(end,:)=[];

% K_mat(end,end)=0;


M_mat1=M_mat;
K_mat1=K_mat;
element_index1=element_index;


coef = 1;

v_s = sqrt(element_index(end,7)/element_index(end,6));
rho = element_index(end,6);

% C(end,end) = C(end,end) + rho*v_s*coef;
%   C(end,end) = rho*v_s*coef;
% C(end-1,end)=0;
% C(end,end-1)=0;


%  C(end,end) = C(end,end) + rho*v_s*coef;
% C(end,end) =  M_mat(end,end)*v_s*coef;

  C(end,end) =rho*v_s*coef;

C1=C;


end