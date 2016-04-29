function [C,output] = damping_mat_gen(output,element_index)

% Last update: 15 February 2016
n_e = size(element_index,1);

C_mat = zeros(n_e+1,n_e+1);

material_mat = output.simulationparams.material_mat;
use_damping = output.simulationparams.damping.use_damping;


switch use_damping
    
    case 1
        
        damping_model='Simplified Rayleigh';
        
        % For more details see Hashash_2001 and Hashash_2002
        
        H   = material_mat(end,5);
        v_s = sqrt(material_mat(:,3)./material_mat(:,2));
        h   = material_mat(:,5)-material_mat(:,4);
        one_over_vs = (1/H)*(sum(h./v_s));
        V_e = 1/one_over_vs;
        f_1 = (V_e)/(4*H);
        w_1 = 2*pi*f_1;
        
        
        for i=1:n_e
            h  = element_index(i,5)-element_index(i,4);
            Mu = element_index(i,7);
            element_damping = element_index(i,8)/100;
            
            c_local(1,1)=1;
            c_local(1,2)=-1;
            c_local(2,1)=-1;
            c_local(2,2)=1;
            
            c_local_m=c_local*(1/h)*Mu*element_damping*2/w_1;
            
            C_mat(i,i)     = C_mat(i,i)     +  c_local_m(1,1);
            C_mat(i,i+1)   = C_mat(i,i+1)   +  c_local_m(1,2);
            C_mat(i+1,i)   = C_mat(i+1,i)   +  c_local_m(2,1);
            C_mat(i+1,i+1) = C_mat(i+1,i+1) +  c_local_m(2,2);
            
            
        end
        
        C = C_mat;
        
    case 2
        
        damping_model='Freq-Independent Rayleigh';
        
        w_n = 2*pi*0.75;  % Assuming first mode is 1 Hz
        w_m = 2*pi*3.75; % Assuming second mode is 10 Hz
        
        alpha_1 = 2*(w_n*w_m)/(w_n+w_m);
        beta_1  = 2*(1/(w_n+w_m));
        n_e = size(element_index,1);
        
        % Stiffness matrix (Part Beta)
        
        for i=1:n_e
            h  = element_index(i,5)-element_index(i,4); % defining h here will help the to define elements with different size.
            Mu = element_index(i,7);
            element_damping = element_index(i,8)/100;
            
            c_local(1,1)=1;
            c_local(1,2)=-1;
            c_local(2,1)=-1;
            c_local(2,2)=1;
            
            c_local_m=c_local*(1/h)*Mu*element_damping;
            
            C_mat(i,i)     = C_mat(i,i)     +  c_local_m(1,1);                                     
            C_mat(i,i+1)   = C_mat(i,i+1)   +  c_local_m(1,2);
            C_mat(i+1,i)   = C_mat(i+1,i)   +  c_local_m(2,1);
            C_mat(i+1,i+1) = C_mat(i+1,i+1) +  c_local_m(2,2);
            
            
        end
        
        C_mat = beta_1*C_mat;
        
        
        % Mass Matrix (Part Alpha)
        M_mat_c = zeros(n_e+1,n_e+1);
        
        for i=1:n_e
            h  = element_index(i,5)-element_index(i,4); % defining h here will help the to define elements with different size.
            rho = element_index(i,6);
            element_damping = element_index(i,8)/100;
            
            m_local(1,1)=0.3333;
            m_local(1,2)=0.1666;
            m_local(2,1)=0.1666;
            m_local(2,2)=0.3333;
                
            
            m_local_m=m_local*(h)*rho*element_damping;
            
            
            M_mat_c(i,i)     = M_mat_c(i,i)     + m_local_m(1,1);
            M_mat_c(i,i+1)   = M_mat_c(i,i+1)   + m_local_m(1,2);
            M_mat_c(i+1,i)   = M_mat_c(i+1,i)   + m_local_m(2,1);
            M_mat_c(i+1,i+1) = M_mat_c(i+1,i+1) + m_local_m(2,2);
            
        end
        
        M_mat_c = M_mat_c * alpha_1;
        
        C = C_mat + M_mat_c;
        
             
    case 3 
        
        damping_model='Freq-Independent Rayleigh';
        
        
    case 4
        
        damping_model='BKT';
        
    case 5
        
        damping_model='None';
        
        C=C_mat;
end

F1 = sprintf('%s%s%s','-------> ',damping_model,' damping model is selected. Generating the damping matrix ...' );
disp(F1);




end