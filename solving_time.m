function output = solving_time(output,M_inv,M_mat,K,C,element_index,acc_vec_1,solution_type,Force,nl_e_n)

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

if strcmp(use_damping,'BKT')==1 || strcmp(use_damping,'BKT2')==1 || strcmp(use_damping,'BKT3')==1 || strcmp(use_damping,'BKT3F')==1
    
    unit_vec=ones(n_e+1,1);
    
    
    switch use_damping
        
        case 'BKT'
            
            % generating alpha, beta, gamma for each element. (Q= 1 / (2*damping))
            
            damping_zeta = element_index(:,8);  % damping
            Q = 1./ (2*damping_zeta);
            
            %                Qo,      alpha_1,      alpha_2,      gamma_1,      gamma_2,       beta
            Qtable = [
                5.00    0.211111102  0.236842104  0.032142857  0.271428571  0.1400000
                6.25    0.188888889  0.184210526  0.039893617  0.336879433  0.1015200
                8.33    0.157777778  0.139473684  0.045000000  0.380000000  0.0700000
                10.00   0.137777765  0.121052630  0.032942899  0.278184480  0.0683000
                15.00   0.097777765  0.081052630  0.032942899  0.278184480  0.0450000
                20.00   0.078139527  0.060526314  0.031409788  0.277574872  0.0342250
                25.00   0.064285708  0.049999999  0.031578947  0.285714286  0.0266000
                30.00   0.053658537  0.044736842  0.026640676  0.246913580  0.0230850
                35.00   0.046341463  0.038157895  0.027098480  0.251156642  0.0196690
                40.00   0.040487805  0.034210526  0.025949367  0.240506329  0.0173800
                45.00   0.036585366  0.028947368  0.031393568  0.290964778  0.0143660
                50.00   0.032926829  0.026315789  0.032488114  0.301109350  0.0126200
                60.00   0.027900000  0.022300000  0.027500000  0.254500000  0.0114000
                70.00   0.024000000  0.019000000  0.032488114  0.301109350  0.0083000
                80.00   0.020700000  0.017400000  0.025100000  0.232600000  0.0088000
                90.00   0.018700000  0.015400000  0.024400000  0.225600000  0.0079000
                100.00  0.017000000  0.014000000  0.028021016  0.288966725  0.0062810
                120.00  0.014200000  0.011500000  0.028000000  0.270000000  0.0052000
                150.00  0.011400000  0.009400000  0.024000000  0.231600000  0.0047000
                200.00  0.008500000  0.007050000  0.022603978  0.226039783  0.0035392
                250.00  0.006900000  0.005500000  0.026900000  0.259600000  0.0027000
                300.00  0.005700000  0.004700000  0.027072758  0.279187817  0.0021276
                350.00  0.004800000  0.004000000  0.024200000  0.233900000  0.0020000
                400.00  0.004300000  0.003600000  0.021425572  0.214255718  0.0017935
                450.00  0.003900000  0.003000000  0.028000000  0.271000000  0.0015000
                500.00  0.003500000  0.002850000  0.023408925  0.241404535  0.0013670
                
                ];
            
               
            % loop over all elements 
            
            alpha_1 = zeros(size(Q,1),1);
            alpha_2 = zeros(size(Q,1),1);
            gamma_1 = zeros(size(Q,1),1);
            gamma_2 = zeros(size(Q,1),1);
            beta    = zeros(size(Q,1),1);
            
            for ie = 1 : size(Q,1)
            
                q_index = Search_Quality_Table(Qtable,Q(ie,1));
                
                if q_index == -1
                  alpha_1(ie,1) = 0; 
                  alpha_2(ie,1) = 0; 
                  gamma_1(ie,1) = 0;
                  gamma_2(ie,1) = 0;
                  beta(ie,1)    = 0;
                else
                  alpha_1(ie,1) = Qtable(q_index,2);  
                  alpha_2(ie,1) = Qtable(q_index,3);   
                  gamma_1(ie,1) = Qtable(q_index,4);  
                  gamma_2(ie,1) = Qtable(q_index,5); 
                  beta(ie,1)    = Qtable(q_index,6); 
                    
                end
                
                               
            end
            
                        % this part does not make sense, check it later
            
            alpha_1=[alpha_1;alpha_1(end,1)];
            alpha_2=[alpha_2;alpha_2(end,1)];
            gamma_1=[gamma_1;gamma_1(end,1)];
            gamma_2=[gamma_2;gamma_2(end,1)];
            beta=[beta;beta(end,1)];

            
            % Memory allocation
            
            %             sai_1=zeros(n_e+1,nt_step);
            %             sai_2=zeros(n_e+1,nt_step);
            
            sai_1=zeros(n_e+1,1);
            sai_2=zeros(n_e+1,1);
            
            
            for tt = 3 : nt_step
                
                %                 F  = M_mat*unit_vec*acc_vec(tt,2);
                % adding direct force
                
                Force_1=Force.Force_1;
                Force_2=Force.Force_2;
                
                F = zeros(n_e+1,1);
                F(n_e+1,1) = Force_1(tt,1);
                F(n_e,1) = -Force_2(tt,1);
                
                %                 Storing the whole matrix is unnecessary.
                
                %                 sai_1(:,tt-1) = (dt/2)*((1-dt*gamma_1)*u(:,tt-1)+u(:,tt-2))+exp(-gamma_1*dt)*(sai_1(:,tt-2));
                %                 sai_2(:,tt-1) = (dt/2)*((1-dt*gamma_2)*u(:,tt-1)+u(:,tt-2))+exp(-gamma_2*dt)*(sai_2(:,tt-2));
                
                sai_1 = (dt/2)*((1-dt*gamma_1).*u(:,tt-1)+u(:,tt-2))+exp(-gamma_1*dt).*(sai_1);
                sai_2 = (dt/2)*((1-dt*gamma_2).*u(:,tt-1)+u(:,tt-2))+exp(-gamma_2*dt).*(sai_2);
                
                
                
                %                 KU = K * (u(:,tt-1) + beta .* (u(:,tt-1)-u(:,tt-2))./dt -...
                %                     alpha_1.*gamma_1.*sai_1(:,tt-1)-...
                %                     alpha_2.*gamma_2.*sai_2(:,tt-1));
                
                KU = K * (u(:,tt-1) + beta .* (u(:,tt-1)-u(:,tt-2))./dt -...
                    sum(alpha_1.*gamma_1.*sai_1)-...
                    sum(alpha_2.*gamma_2.*sai_2));
                
                CU_dot = C*(u(:,tt-1)-u(:,tt-2))/dt;
                
                u(:,tt) = M_inv * (F  - KU - CU_dot) * dt^2 +2*u(:,tt-1) - u(:,tt-2);
                
                
            end
            
%             xrt = 1;
            
            
            
            
        case 'BKT2'
            
            % generating alpha, beta, gamma for each element. (Q= 1 / (2*damping))
            
            
            alpha_1 = output.damping_param.bkt2.alpha_1;
            alpha_2 = output.damping_param.bkt2.alpha_2;
            
            gamma_1 = output.damping_param.bkt2.gamma_1;
            gamma_2 = output.damping_param.bkt2.gamma_2;
            
            beta = output.damping_param.bkt2.beta;
            
            % Memory allocation
            
            %             sai_1=zeros(n_e+1,nt_step);
            %             sai_2=zeros(n_e+1,nt_step);
            
            sai_1=zeros(n_e+1,1);
            sai_2=zeros(n_e+1,1);
            
            
            for tt = 3 : nt_step
                
                %                 F  = M_mat*unit_vec*acc_vec(tt,2);
                % adding direct force
                
                Force_1=Force.Force_1;
                Force_2=Force.Force_2;
                
                F = zeros(n_e+1,1);
                F(n_e+1,1) = Force_1(tt,1);
                F(n_e,1) = -Force_2(tt,1);
                
                %                 Storing the whole matrix is unnecessary.
                
                %                 sai_1(:,tt-1) = (dt/2)*((1-dt*gamma_1)*u(:,tt-1)+u(:,tt-2))+exp(-gamma_1*dt)*(sai_1(:,tt-2));
                %                 sai_2(:,tt-1) = (dt/2)*((1-dt*gamma_2)*u(:,tt-1)+u(:,tt-2))+exp(-gamma_2*dt)*(sai_2(:,tt-2));
                
             
              sai_1 = (dt/2)*((1-dt*gamma_1).*u(:,tt-1)+u(:,tt-2))+exp(-gamma_1*dt).*(sai_1);
              sai_2 = (dt/2)*((1-dt*gamma_2).*u(:,tt-1)+u(:,tt-2))+exp(-gamma_2*dt).*(sai_2);
                

                
                
                
                %                 KU = K * (u(:,tt-1) + beta .* (u(:,tt-1)-u(:,tt-2))./dt -...
                %                     alpha_1.*gamma_1.*sai_1(:,tt-1)-...
                %                     alpha_2.*gamma_2.*sai_2(:,tt-1));
                
                KU = K * (u(:,tt-1) + beta .* (u(:,tt-1)-u(:,tt-2))./dt -...
                     (alpha_1.*gamma_1.*sai_1)-...
                     (alpha_2.*gamma_2.*sai_2));
                
                CU_dot = C*(u(:,tt-1)-u(:,tt-2))/dt;
                
                u(:,tt) = M_inv * (F  - KU - CU_dot) * dt^2 +2*u(:,tt-1) - u(:,tt-2);
                
                
            end
            
            xrt = 1;
        case 'BKT3'
            
            % generating alpha, beta, gamma for each element. (Q= 1 / (2*damping))
            
            damping_zeta = element_index(:,8);  % damping
            Q = 1./ (2*damping_zeta/100);
            
            f_max=10;
            
            gamma_1 = 0.0151 * 2 * 3.14 * f_max;
            gamma_2 = 0.1    * 2 * 3.14 * f_max;
            gamma_3 = 0.4814 * 2 * 3.14 * f_max;
            
            alpha_1 = (-2.723  * Q.^-0.8206 + 1.601)./Q;
            alpha_2 = (-1.439  * Q.^-0.9668 + 1.040)./Q;
            alpha_3 = (-0.3037 * Q.^-0.8911 + 1.032)./Q;
            
            
            beta = (0.1249*Q.^(-0.804)+0.4782)./(Q*2*3.14*f_max);
            
            % this part does not make sense, check it later
            
            alpha_1=[alpha_1;alpha_1(end,1)];
            alpha_2=[alpha_2;alpha_2(end,1)];
            alpha_3=[alpha_3;alpha_3(end,1)];
            beta=[beta;beta(end,1)];
            
            % Memory allocation
            
            %             sai_1=zeros(n_e+1,nt_step);
            %             sai_2=zeros(n_e+1,nt_step);
            %             sai_3=zeros(n_e+1,nt_step);
            
            sai_1=zeros(n_e+1,1);
            sai_2=zeros(n_e+1,1);
            sai_3=zeros(n_e+1,1);
            
            for tt = 3 : nt_step
                
                F  = M_mat*unit_vec*acc_vec(tt,2);
                
                %                 Storing the whole matrix is unnecessary.
                %                 sai_1(:,tt-1) = (dt/2)*((1-dt*gamma_1)*u(:,tt-1)+u(:,tt-2))+exp(-gamma_1*dt)*(sai_1(:,tt-2));
                %                 sai_2(:,tt-1) = (dt/2)*((1-dt*gamma_2)*u(:,tt-1)+u(:,tt-2))+exp(-gamma_2*dt)*(sai_2(:,tt-2));
                %                 sai_3(:,tt-1) = (dt/2)*((1-dt*gamma_3)*u(:,tt-1)+u(:,tt-2))+exp(-gamma_3*dt)*(sai_3(:,tt-2));
                
                sai_1 = (dt/2)*((1-dt*gamma_1)*u(:,tt-1)+u(:,tt-2))+exp(-gamma_1*dt)*(sai_1);
                sai_2 = (dt/2)*((1-dt*gamma_2)*u(:,tt-1)+u(:,tt-2))+exp(-gamma_2*dt)*(sai_2);
                sai_3 = (dt/2)*((1-dt*gamma_3)*u(:,tt-1)+u(:,tt-2))+exp(-gamma_3*dt)*(sai_3);
                
                
                
                %                 KU = K * (u(:,tt-1) + beta .* (u(:,tt-1)-u(:,tt-2))./dt -...
                %                     alpha_1.*gamma_1.*sai_1(:,tt-1)-...
                %                     alpha_2.*gamma_2.*sai_2(:,tt-1)-...
                %                     alpha_3.*gamma_3.*sai_3(:,tt-1));
                
                
                KU = K * (u(:,tt-1) + beta .* (u(:,tt-1)-u(:,tt-2))./dt -...
                    alpha_1.*gamma_1.*sai_1-...
                    alpha_2.*gamma_2.*sai_2-...
                    alpha_3.*gamma_3.*sai_3);
                
                CU_dot = C*(u(:,tt-1)-u(:,tt-2))/dt;
                
                u(:,tt) = M_inv * (F  - KU - CU_dot) * dt^2 +2*u(:,tt-1) - u(:,tt-2);
            end
            
        case 'BKT3F'
            
            
            % generating alpha, beta, gamma for each element. (Q= 1 / (2*damping))
            
            damping_zeta = element_index(:,8);  % damping
            Q = 1./ (2*damping_zeta/100);
            
            f_max=10;
            
            gamma_1 = 0.002 * 2 * 3.14 * f_max;
            gamma_2 = 0.0116 * 2 * 3.14 * f_max;
            gamma_3 = 0.0798 * 2 * 3.14 * f_max;
            
            alpha_1 = (-2.809  * Q.^-0.7919 + 1.512)./Q;
            alpha_2 = (-1.748  * Q.^-0.8820 + 1.064)./Q;
            alpha_3 = (-2.358 * Q.^-1.725 + 1.581)./Q;
            
            
            beta = (0.09232*Q.^(-0.8876)+0.006941)./(Q*2*3.14*f_max);
            
            % this part does not make sense, check it later
            
            alpha_1=[alpha_1;alpha_1(end,1)];
            alpha_2=[alpha_2;alpha_2(end,1)];
            alpha_3=[alpha_3;alpha_3(end,1)];
            beta=[beta;beta(end,1)];
            
            % Memory allocation
            
            %             sai_1=zeros(n_e+1,nt_step);
            %             sai_2=zeros(n_e+1,nt_step);
            %             sai_3=zeros(n_e+1,nt_step);
            
            sai_1=zeros(n_e+1,1);
            sai_2=zeros(n_e+1,1);
            sai_3=zeros(n_e+1,1);
            
            for tt = 3 : nt_step
                
                F  = M_mat*unit_vec*acc_vec(tt,2);
                
                %                 Storing the whole matrix is unnecessary.
                %                 sai_1(:,tt-1) = (dt/2)*((1-dt*gamma_1)*u(:,tt-1)+u(:,tt-2))+exp(-gamma_1*dt)*(sai_1(:,tt-2));
                %                 sai_2(:,tt-1) = (dt/2)*((1-dt*gamma_2)*u(:,tt-1)+u(:,tt-2))+exp(-gamma_2*dt)*(sai_2(:,tt-2));
                %                 sai_3(:,tt-1) = (dt/2)*((1-dt*gamma_3)*u(:,tt-1)+u(:,tt-2))+exp(-gamma_3*dt)*(sai_3(:,tt-2));
                
                sai_1 = (dt/2)*((1-dt*gamma_1)*u(:,tt-1)+u(:,tt-2))+exp(-gamma_1*dt)*(sai_1);
                sai_2 = (dt/2)*((1-dt*gamma_2)*u(:,tt-1)+u(:,tt-2))+exp(-gamma_2*dt)*(sai_2);
                sai_3 = (dt/2)*((1-dt*gamma_3)*u(:,tt-1)+u(:,tt-2))+exp(-gamma_3*dt)*(sai_3);
                
                
                
                %                 KU = K * (u(:,tt-1) + beta .* (u(:,tt-1)-u(:,tt-2))./dt -...
                %                     alpha_1.*gamma_1.*sai_1(:,tt-1)-...
                %                     alpha_2.*gamma_2.*sai_2(:,tt-1)-...
                %                     alpha_3.*gamma_3.*sai_3(:,tt-1));
                
                
                KU = K * (u(:,tt-1) + beta .* (u(:,tt-1)-u(:,tt-2))./dt -...
                    alpha_1.*gamma_1.*sai_1-...
                    alpha_2.*gamma_2.*sai_2-...
                    alpha_3.*gamma_3.*sai_3);
                
                CU_dot = C*(u(:,tt-1)-u(:,tt-2))/dt;
                
                u(:,tt) = M_inv * (F  - KU - CU_dot) * dt^2 +2*u(:,tt-1) - u(:,tt-2);
            end
            
    end

elseif strcmp(use_damping,'Nonlinear') == 1
    
    
                % F  = M_mat*unit_vec*acc_vec(tt,2);
                % adding direct force
    
                Force_1=Force.Force_1;
                Force_2=Force.Force_2;
                
%                     Force_1= Force_1(433:end,1);
%                     Force_2= Force_2(433:end,1);
            
                F = zeros(n_e+1,1);
                Fnl = F * 0;         % force from nonlinear contribution

                
                % strain matrix
                
                %strain_matrix = zeros(n_e,1);
                strain_matrix_0 = zeros(n_e,1);
                

% %%%%%%%%% with k matrix involved  (strat)              
%             for tt = 3 : nt_step
%                 
% 
%                 F(n_e+1,1) = Force_1(tt,1);
%                 F(n_e,1)   = -Force_2(tt,1);
%                 
%                 KU = K * (u(:,tt-1));
%                 
%                 CU_dot = (C*(u(:,tt-1)-u(:,tt-2))/dt);
%                 
%                 u(:,tt) = M_inv * (F + Fnl  - KU - CU_dot) * dt^2 +2*u(:,tt-1) - u(:,tt-2);
%                 
%                 % extract total strain
%                 strain_matrix = extract_strain_stress_nonlin(u(:,tt),output);
%                 
%                 % placticity assumptions
%                 
%                 ep = zeros(3,3);
%                 ep_barn = 0;
%                 alpha_n = zeros(3,3);
%                 sigma0 = zeros(3,3);
%                 soo = zeros(3,3);
% 
%                 
%                 for il = 1: nl_e_n % only generate force for nonlinear elements
%                     k =  100;
%                     psi = 0.15;
%                     mu  = element_index(il,7); % N/m2 (kgm/s2)
%                     H_kin = psi*mu;
%                     Su = 10000; % (underained strength(N)
%                     nu = 0.3; % poisson
%                     E= mu*(2*(1+nu));
%                     lambda = E*nu/((1+nu)*(1-2*nu));
%                     kappa  = lambda + 2*mu/3;
% 
%                                  
%                     G = mu; %KN/m2 
%                     Ez  = strain_matrix(il,1);
%                     Ez0 = strain_matrix_0(il,1); % strain of t-1 
%                     e_n = [0 Ez 0; Ez 0 0; 0 0 0];
%                     e_n0= [0 Ez0 0; Ez0 0 0; 0 0 0];
%                     
%                     [fs,epl,sigma,alpha, ep_bar, dl, load_unl] = ...
%                         vonMises_ArmstrongFrederick_KinematicHrd_D2(k,H_kin, Su, kappa, G, e_n, ep, ep_barn, alpha_n, sigma0, soo, e_n0);
%                     
%                     strain_matrix_0 = strain_matrix;
%                     
%                     % update Fnl
%                     
%                     w_g = 1;
%                     el_size = element_index(il,5) - element_index(il,4);
%                     B = [-1;1];%/el_size;
%                     Fnl_local = w_g .* (B*sigma(2,1)); 
%                     Fnl(il:il+1,1) = Fnl(il:il+1,1)+Fnl_local;
%                     
%                   
%                     
%                 end
%             end
            
%%%%%%%%% with k matrix involved  (end)  
    
%%%%%%%%% without k matrix involved  (strat)      


                % placticity assumptions
                
                ep = zeros(3,3);
                ep_barn = 0;
                alpha_n = zeros(3,3);
                sigma0 = zeros(3,3);
                soo = zeros(3,3);



            for tt = 3 : nt_step
                 

                F(n_e+1,1) = Force_1(tt,1);
                F(n_e,1)   = -Force_2(tt,1);
                
                 KU = K * (u(:,tt-1));
                
%                  CU_dot = (C*(u(:,tt-1)-u(:,tt-2))/dt);
                
%                 u(:,tt) = M_inv * (F + Fnl  - KU - CU_dot) * dt^2 +2*u(:,tt-1) - u(:,tt-2);
                 u(:,tt) = M_inv * (F - KU) * dt^2 +2*u(:,tt-1) - u(:,tt-2);
                
                % extract total strain
                strain_matrix = extract_strain_stress_nonlin(u(:,tt-1),output);
                
                

                
                for il = 1: n_e % only generate force for nonlinear elements
%                     k =  1000;
%                     psi = 0.15;
                      mu  = element_index(il,7); % N/m2 (kgm/s2)
%                     H_kin = psi*mu;
%                     Su = 10000; % (underained strength(N)
%                     nu = 0.3; % poisson
%                     E= mu*(2*(1+nu));
%                     lambda = E*nu/((1+nu)*(1-2*nu));
%                     kappa  = lambda + 2*mu/3;
% 
%                                  
%                     G = mu; %KN/m2 
                    Ez  = strain_matrix(il,1);
%                     Ez0 = strain_matrix_0(il,1); % strain of t-1 
%                     e_n = [0 Ez 0; Ez 0 0; 0 0 0];
%                     e_n0= [0 Ez0 0; Ez0 0 0; 0 0 0];
                    
%                     [fs,epl,sigma,alpha, ep_bar, dl, load_unl] = ...
%                         vonMises_ArmstrongFrederick_KinematicHrd_D2(k,H_kin, Su, kappa, G, e_n, ep, ep_barn, alpha_n, sigma0, soo, e_n0);

                    sigma = simple_rheology(Ez,mu);

%                     strain_matrix_0 = strain_matrix;
                    
                    
%                     ep      = epl;
%                     alpha_n = alpha;
%                     sigma0  = zeros(3,3);
%                     ep_barn = ep_bar;
                    
                    % update Fnl
                    
                    % w_g = 1;
                    % B = [-1;1];
                    % Fnl_local = sigma(1,1).*B; 
                    Fnl_local = [-sigma(1,1);sigma(1,1)]; 
                    Fnl(il:il+1,1) = Fnl(il:il+1,1)+Fnl_local;

                    
                end
% %                 
%                 figure(1)
%                 plot(KU);
%                 hold on 
%                 plot(Fnl,'r');
%                 pause;
%                 close all 
            end
            
%%%%%%%%% without k matrix involved  (end)  

    
    
else
    
    fem_sol='implicit';
    %       fem_sol='explicit';
    
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



function q_index = Search_Quality_Table(Qtable,Q)
if ( Q > 525 )
    q_index = -1;
end

if ( ( Q > 475 ) && ( Q <= 525 ) )
    q_index = size(Qtable,1);
end

if (Q <= 475)
    
%     range = int(Q / 5);
    min = 1000;
    
    for i = 2:size(Qtable,1)
        
        diff = abs(Q - QTABLE(i,1));
        
        
        if(diff < min)
            min = diff;
            
        else
            
            q_index = i-1;
        end
    end
    
end


end
