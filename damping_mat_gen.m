function [C,output] = damping_mat_gen(output,element_index)

% Last update: 15 February 2016
n_e = size(element_index,1);

C_mat = zeros(n_e+1,n_e+1);

material_mat = output.simulationparams.material_mat;
use_damping = output.simulationparams.damping.use_damping;


switch use_damping
    
    case 'SRD'
        
       
        
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
        
    case 'RD2'
        
        
        
        f_1 = 2.0;
        f_2 = 6.25;
        
        w_n = 2*pi*f_1;  % Assuming first mode is 1 Hz
        w_m = 2*pi*f_2; % Assuming second mode is 10 Hz
        
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
        
        
        
        % temporarily the value for plotting the damping is calculated for 
        % the first element.
        
        f = 0.01:0.01:2.4;
        w = 2*pi*f;
        sai_mass = alpha_1./(2*w)*element_index(1,8)/100;
        sai_stif = beta_1*w/2*element_index(1,8)/100;
        sai_tot  = sai_mass + sai_stif; 
        sai_mat  = [w' sai_mass' sai_stif' sai_tot'];
        output.simulationparams.damping.FIRayleigh = sai_mat;
        output.simulationparams.damping.w12=[w_n w_m];
        
        
        
    case 'FDRD'    
        
  
    case 'BKT'    
         C=C_mat;
         
    case 'BKT2'
        
         C=C_mat;
         
         % -------
            damping_zeta = element_index(:,8);  % damping
            Q = 1./ (2*damping_zeta/100);
            
            f_max=2.4;
            
            gamma_1 = 0.0373 * 2 * 3.14 * f_max;
            gamma_2 = 0.3082 * 2 * 3.14 * f_max;
            
            alpha_1 = (-2.656  * Q.^-0.8788 + 1.677)./Q;
            alpha_2 = (-0.5623 * Q.^-1.03   + 1.262)./Q;
            
            beta = (0.1876*Q.^(-0.9196)+0.6137)./(Q*2*3.14*f_max);
            
            % this part does not make sense, check it later
            
            alpha_1=[alpha_1;alpha_1(end,1)];
            alpha_2=[alpha_2;alpha_2(end,1)];
            beta=[beta;beta(end,1)];
            
            
            % shear velocity adjustment
            % --------------------------
            reference_f = 0;
            
            w = reference_f ./ f_max;
		    w2 = w .* w;
            
            shear_vel_corr_factor = sqrt( 1. - (alpha_1 * gamma_1 * gamma_1 / (gamma_1 * gamma_1 + w2) + alpha_2 * gamma_2 * gamma_2 / (gamma_2 * gamma_2 + w2)));
            shear_mod_corr_factor = shear_vel_corr_factor.^2;
            
%              element_index(:,7) = shear_mod_corr_factor(1:end-1,:) .* element_index(:,7);

            element_index(:,7) = 1*element_index(:,7);

            bkt2.alpha_1 = alpha_1;
            bkt2.alpha_2 = alpha_2;
            bkt2.gamma_1 = gamma_1;
            bkt2.gamma_2 = gamma_2;
            bkt2.beta    = beta;
            
            output.element_index = element_index;
            output.damping_param.bkt2 = bkt2;
        
    case 'BKT3'
        C=C_mat;
        
    case 'BKT3F'    
         C=C_mat;
        
    case 'None'
        
        
        
        C=C_mat;
end

F1 = sprintf('%s%s%s','-------> ',use_damping,' damping model is selected. Generating the damping matrix ...' );
disp(F1);



end