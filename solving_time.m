function output = solving_time(output,M_inv,M_mat,K,C,element_index,acc_vec_1,solution_type)

% Input could be acc : acceleration  disp: displacement.

disp('-------> Solving in the time domain ...');

pause(0.01); % Give a small break to remove figures and unnecessary variables.

n_e = size(element_index,1);


nt_step  = output.simulationparams.n_timest;
sim_time = output.simulationparams.sim_time;
dt       = output.simulationparams.dt;

u=zeros(n_e+1,nt_step);
v=zeros(n_e+1,nt_step);
a=zeros(n_e+1,nt_step);
v_t=zeros(n_e+1,nt_step);
u_t=zeros(n_e+1,nt_step);
F  =zeros(n_e+1,nt_step);


use_damping = output.simulationparams.damping.use_damping;

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
    
    acc_vec_temp_1 = acc_vec(end,1)+dt:dt:sim_time-dt;
    acc_vec_temp_2 = zeros(size(acc_vec_temp_1,2),1);
    acc_vec_temp = [acc_vec_temp_1' acc_vec_temp_2];
    acc_vec = [acc_vec;acc_vec_temp];
    
end

output.acc_temp = acc_vec;

% Mu..+cu.+ku=-mug..;
acc_vec(:,2)=acc_vec(:,2)*+1;


input_acc = acc_vec * 1;
acc_gr    = input_acc(:,2)';
vel_gr    = cumtrapz(acc_gr)*dt;
disp_gr   = cumtrapz(vel_gr)*dt;

if strcmp(use_damping,'BKT2')==1 || strcmp(use_damping,'BKT3')==1 || strcmp(use_damping,'BKT3F')==1
    
    unit_vec=ones(n_e+1,1);
    
    
    switch use_damping
        
        case 'BKT2'
            
            % generating alpha, beta, gamma for each element. (Q= 1 / (2*damping))
            
            damping_zeta = element_index(:,8);  % damping
            Q = 1./ (2*damping_zeta/100);
            
            f_max=10;

            gamma_1 = 0.0373 * 2 * 3.14 * f_max;
            gamma_2 = 0.3082 * 2 * 3.14 * f_max;

            alpha_1 = (-2.656  * Q.^-0.8788 + 1.677)./Q;
            alpha_2 = (-0.5623 * Q.^-1.03   + 1.262)./Q;

            beta = (0.1876*Q.^(-0.9196)+0.6137)./(Q*2*3.14*f_max);
            
            % this part does not make sense, check it later
            
            alpha_1=[alpha_1;alpha_1(end,1)];
            alpha_2=[alpha_2;alpha_2(end,1)];
            beta=[beta;beta(end,1)];
            
            % Memory allocation
            
            sai_1=zeros(n_e+1,nt_step);
            sai_2=zeros(n_e+1,nt_step);
            
          
            for tt = 3 : nt_step
                
                F  = M_mat*unit_vec*acc_vec(tt,2);
                
%                 Storing the whole matrix is unnecessary.
                sai_1(:,tt-1) = (dt/2)*((1-dt*gamma_1)*u(:,tt-1)+u(:,tt-2))+exp(-gamma_1*dt)*(sai_1(:,tt-2));
                sai_2(:,tt-1) = (dt/2)*((1-dt*gamma_2)*u(:,tt-1)+u(:,tt-2))+exp(-gamma_2*dt)*(sai_2(:,tt-2));
                
                KU = K * (u(:,tt-1) + beta .* (u(:,tt-1)-u(:,tt-2))./dt -...
                    alpha_1.*gamma_1.*sai_1(:,tt-1)-...
                    alpha_2.*gamma_2.*sai_2(:,tt-1));
                
                CU_dot = C*(u(:,tt-1)-u(:,tt-2))/dt;
                
                u(:,tt) = M_inv * (F  - KU - CU_dot) * dt^2 +2*u(:,tt-1) - u(:,tt-2);
                
                
            end
            
        case 'BKT3'
            
        case 'BKT3F'
            
    end
    
else
    
%     fem_sol='implicit';
    fem_sol='explicit';
    
    unit_vec=ones(n_e+1,1);
    
    if strcmp(solution_type,'acc')==1
        
        if strcmp(fem_sol,'explicit')==1
         
            
            for tt = 3 : nt_step
                
                
                a1= -1*K*u(:,tt-1)-C*(u(:,tt-1)-u(:,tt-2))/dt;
                a2= M_mat*unit_vec*acc_vec(tt,2);
                a3= M_inv*(dt^2)*(a2+a1);
                u(:,tt)= a3+2*u(:,tt-1)-u(:,tt-2);
                
                %              u(end,tt)=0;
            end

        elseif strcmp(fem_sol,'implicit')==1
            
            % Newmark values
            beta  = 1/4;
            gamma = 1/2;
            
            a(:,1)=unit_vec*acc_vec(1,2);
            
            for tt = 2 : nt_step
                F(:,tt)   = M_mat*unit_vec*acc_vec(tt,2);
                u_t(:,tt) = u(:,tt-1)+dt*v(:,tt-1) + (1/2)*dt^2*(1-2*beta)*a(:,tt-1);
                v_t(:,tt) = v(:,tt-1)+dt*a(:,tt-1)*(1-gamma);
                RH        = F(:,tt) - C*v_t(:,tt)-K*u_t(:,tt);
                W         = M_mat+C*gamma*dt+K*beta*dt^2;
                a(:,tt)   = W\RH;
                v(:,tt)   = v(:,tt-1)+dt*((1-gamma)*a(:,tt-1)+gamma*a(:,tt));
                u(:,tt)   = u(:,tt-1) + dt*v(:,tt-1)+ (dt^2)*(1/2)*((1-2*beta)*a(:,tt-1)+2*beta*a(:,tt));
                
            end
            
            
            
        end
        
        
    elseif strcmp(solution_type,'disp')==1
        
        for tt = 3 : nt_step
            a1= -1*K*u(:,tt-1)-C*(u(:,tt-1)-u(:,tt-2))/dt;
            a3= M_inv*(dt^2)*(a1);
            u(:,tt)= a3+2*u(:,tt-1)-u(:,tt-2);
            u(end,tt)=disp_gr(1,tt);
        end
        
    end
    
    
end

output.nodetime=u;

disp('-------> End of solving in time domain.');


end