function u = solving_time(sim_time,dt,M_inv,M_mat,K,C,element_index,acc_vec_1,input)

% Input could be acc : acceleration  disp: displacement.

disp('-------> Solving in the time domain ...');

pause(0.01); % Give a small break to remove figures and unnecessary variables.

n_e = size(element_index,1);
nt_step = floor(sim_time/dt);
u=zeros(n_e+1,nt_step);


% Initial condition.
% Inorder to define a initial condition for second order problem,
% we need to have 2 initial conditions.
% For the case of acceleration these two values are velocity and
% displacement.

% I don't define the velocity inside time solution. I just define the
% initial velocity and from that velocity I calculate the second time step
% of displacement.

   
   % resampling input file to be according to the simulation dt 
   
%     dt_data = acc_vec_1(2,1)-acc_vec_1(1,1);
          
        data_final_time = acc_vec_1(end,1);
        data_new_time_vec = 0:dt:data_final_time;
        acc_vec_new = interp1(acc_vec_1(:,1),acc_vec_1(:,2),data_new_time_vec);
        acc_vec = [data_new_time_vec' acc_vec_new'];
        
  % zero pading input file incase that there is no enough data to read.      
  if acc_vec(end,1) < sim_time
      
      acc_vec_temp_1 = acc_vec(end,1):dt:sim_time;
      acc_vec_temp_2 = zeros(size(acc_vec_temp_1,2),1);
      acc_vec_temp = [acc_vec_temp_1' acc_vec_temp_2];
      acc_vec = [acc_vec;acc_vec_temp];
      
  end
        
unit_vec=ones(n_e+1,1);


if strcmp(input,'acc')==1
    
    for tt = 3 : nt_step
        
            a1= -1*K*u(:,tt-1)-C*(u(:,tt-1)-u(:,tt-2))/dt;
            a2= M_mat*unit_vec*acc_vec(tt,2);
            a3= M_inv*(dt^2)*(a2+a1);
            
            u(:,tt)= a3+2*u(:,tt-1)-u(:,tt-2);

    end
    
elseif strcmp(input,'disp')==1
    
    for tt = 3 : nt_step
            a1= -1*K*u(:,tt-1)-C*(u(:,tt-1)-u(:,tt-2))/dt;
            a3= M_inv*(dt^2)*(a1);
            
            u(:,tt)= a3+2*u(:,tt-1)-u(:,tt-2);
            u(end,tt)=acc_vec(tt,2);
    end       
       
end


disp('-------> End of solving in time domain.');


end