function [ fs, epl, sigma, alpha, ep_bar, dl, load_unl  ] = vonMises_ArmstrongFrederick_KinematicHrd_D2( k, H_kin, Su, K, G, e_n, ep,ep_barn, alpha_n, sigma0, soo, e_n0 )
%   Returning mapping algorithm for a kinematic hardening vonMises model.
%   Armstrong-Fredericks hardening function defined as:
%   d/dt(alpha) = Hkin * d/dt(epl) - Hnl * alpha * d/dt(ep_bar)
%   Implementation after :
%   (1995) F. AURICCHIOand R.L. TAYLOR.
%   TWO MATERIAL MODELS FOR CYCLIC PLASTICITY: NONLINEAR KINEMATIC HARDENING AND GENERALIZED PLASTICITY

% 	/* INPUTS:
%    * k             : Yield shear stress. It defines the initial elastic zone
%    * H_kin, H_nlin : Parameters of the backstress evolution law.
% 	 * G, K          : Material constants.
% 	 * e_n           : Total strain tensor.
% 	 * ep            : Plastic strain tensor at t-1. Used to compute the predictor state
% 	 * sigma0        : Self-weight tensor.
%    * alpha_n       : Backstress tensor at t-1.
% 	 *
% 	 * OUTPUTS:
% 	 * fs            : Updated yield function value
% 	 * epl           : Updated plastic strain
% 	 * sigma         : Updated stress tensor
% 	 * ep_bar        : Updated equivalent hardening variable
%    * alpha         : Updated backstress tensor
% 	 */
%
%    Notes:
%    From 1D analysis
%    1) H_kin = E
%    2) H_nlin = H_kin/( sqrt(2)*(Su - k) )
%       where Su is the unconfined compresion test value

    
%   tol=1E-10;
%% Initialization of parameters and state variables. Taken from David's implementation
% (i) Given delta e and the state variables at tn, evalaute the elastic trial state
ee_trial     = e_n - ep;
ee_vol_trial = trace(ee_trial);
ee_dev_trial = ee_trial - (1/3)*ee_vol_trial*eye(3,3);

sigma_trial  = K*ee_vol_trial*eye(3,3) + 2*G*ee_dev_trial + sigma0;


p_pr         = trace(sigma_trial)/3*eye(3,3);  % Volumetric predictor
S_pr         = sigma_trial - p_pr;             % Deviatoric predictor


n_pr = S_pr/sqrt(1/2*sum(sum(S_pr.*S_pr)));
toto=sum(sum(n_pr.*(e_n-e_n0)));

    
%% get deviatoric of soo
soo_I = trace(soo);
soo_dev = soo - (1/3)*soo_I*eye(3,3);

%%

% Kinematic Modulous
Sy     = k;
G1     = G + H_kin/2;

% (ii) Compute Predictor state
eta_pr  = S_pr - alpha_n;
% % q_pr    = sqrt(1*sum(sum(eta_pr.*eta_pr)));
q_pr    = sqrt(1/2*sum(sum(eta_pr.*eta_pr))); % checking in terms of sqrt(J2)

% % fs_pr = q_pr-Sy;
tol=1E-10;
fs_pr = q_pr-Sy;

if (fs_pr<=tol) %Elastic or unload state
    
    fs     = q_pr;
    epl    = ep;
    sigma  = sigma_trial - sigma0;
    alpha  = alpha_n;
    ep_bar = ep_barn;
    dl = 0;
    load_unl = toto;

else
    %%H_nlin = sqrt(1/2)*H_kin/(Su-k/sqrt(2));
    H_nlin = sqrt(1/2)*H_kin/(Su-Sy);
    
    Sy=sqrt(2)*Sy;
    
    % compute coeficients of quartic function
    S_ss    = sum(sum(S_pr.*S_pr));
    S_aa    = sum(sum(alpha_n.*alpha_n));
    S_sa    = sum(sum(S_pr.*alpha_n));
   
    C1 = (2*G*H_nlin)^2;
    
    C2 = (4*Sy*G*H_nlin + 8*G*G1)*H_nlin;
    
    C3 = (H_nlin^2)*(Sy^2 - S_ss) + 4*G1^2 + 4*H_nlin*Sy*(G+G1);
    
    C4 = 2*H_nlin*(Sy^2 + S_sa - S_ss) + 4*Sy*G1;
    
    C5 = Sy^2 - S_aa + 2*S_sa - S_ss;
    
    DL = roots([C1 C2 C3 C4 C5]);
    pos=find(DL>=0);
    dl=min(DL(pos));
    
    T_lambda = 1/(1+H_nlin*dl);
    Za = S_pr - T_lambda*alpha_n;
    n  = Za/sqrt((sum(sum(Za.*Za))));
           
    epl    = ep + dl*n;
    ep_bar = ep_barn + dl;
    alpha  = T_lambda*alpha_n + H_kin*T_lambda*dl*n;
    
    Sdev   = S_pr - 2*G*dl*n;
    
    ds    =  Sdev -soo_dev;
    dalpha =  alpha -alpha_n;
    
     load_unl= sum(sum(n.*ds));
    %load_unl= sum(sum(n.*(e_n-e_n0)));
    
    
     if load_unl < 0
         po=90;
     end
    %%[sum(sum(n.*ds))/(H_kin-H_nlin*(sum(sum(n.*alpha)))) dl]
   
    
    Z      = Sdev - alpha;    
    fs     = sqrt(sum(sum(Z.*Z)));
    sigma  = p_pr + Sdev - sigma0;
    
end



end

% 
% %%
% function [ fs, epl, sigma, alpha, ep_bar  ] = cutting_plane(Ko, nu, b, patm, eo, Gam, Lambda, pref, Mc, Me, kb, kb_ex, kd, kd_ex, Cm, ho, Ao, Cf, Fmax, Det, S_o, p_o, m_o, alpha_o,F_o, ev_o, fs_Error)
% %   Cutting plane algorithm for the 1997 Manzari & Dafalias constitutive model for sands.
% 
% %   Implementation after :
% %   2001, Manzari M., and Prachathananukit R.
% %   On integration of a cyclic soil plasticity model
% 
% % 	/* INPUTS:
% % 	 * eo, Ko, nu    : Material constants. K(bulk modulus) is asumed hypo-elastic 
% %                      according to: K=Ko(p/patm)^b, where:
% %                      pref : Reference pressure. Ussually assumed as the atmospheric pressure (100kPa), 
% %                      p    : mean pressure 
% %                      b    : material parameter 0.5<b<1.0 
% %                      Poisson's ratio (nu). Assumed constant
% %                      eo   : Initial void ratio
% 
% % 	 * Det           : Total strain tensor increment.
% %    * S_o           : Deviatoric effective stress tensor at t-1. 
% %    * p_o           : mean pressure at  at t-1. Compression positive
% %    * m_o           : Kinematic variable  at t-1.
% %    * alpha_o       : Backstress tensor   at t-1.
% %    * z_o           : Fabric tensor   at t-1.
% %    * ev_o          : Volumentric strain  at t-1. Compression positive
% 
% 
% % 	 *
% % 	 * OUTPUTS:
% % 	 * fs            : Updated yield function value
% % 	 * epl           : Updated plastic strain
% % 	 * sigma         : Updated stress tensor
% % 	 * ep_bar        : Updated equivalent hardening variable
% %    * alpha         : Updated backstress tensor
% % 	 */
% %
% 
% %% 
% Tol=1e100;
% 
% %% Total volumetric and deviatoric strain increments
% Dvol_T = -trace(Det);
% Ddev_T = Det+Dvol_T/3*eye(3,3);
% 
% %% Update void ratio
% evol_n = ev_o + Dvol_T;
% e      = eo -(1+eo)*evol_n;
% 
% %% Predictor state
% if Dvol_T ~= 0     % If Drained
%     p_k     = ( p_o^(1-b) + Ko*(1-b)/(patm^b)*(Dvol_T) )^(1/(1-b));
%     K_k     = (p_k - p_o)/(Dvol_T);
%     G_k     = 3*(1-2*nu)/(2*(1+nu))*K_k;
% else              % If Undrained
%     p_k     = p_o;
%     G_k     = 3*(1-2*nu)/(2*(1+nu))*Ko*(p_k/patm)^b;
% end
% S_k     = S_o + 2*G_k*(Ddev_T);
% alpha_k = alpha_o;
% 
% m_k     = m_o;
% F_k     = F_o;
% 
% f_k     = get_yielfnc( S_k, p_k, alpha_k, m_k );
% 
% %% Corrector state
% if f_k>0
%     
%     % Initialize plastic variables    
%     D_PlM  = 0;  % plastic multiplier increment   
%     DLam   = 0;
%     Tk     = zeros(3,3);
%     
%     while Tol> fs_Error
%        
%         % get yield function
%         f_k     = get_yielfnc( S_k, p_k, alpha_k, m_k );
% 
%         % get elastic constants
%         K_k = Ko*(p_k/patm)^b;
%         G_k = 3*(1-2*nu)/(2*(1+nu))*K_k;
%         
%         % update void ratio
%         psi    = eo - (Gam - Lambda*log(p_k/pref));
%         
%         % Get variable Nk
%         n_k  = get_unitTens( S_k - alpha_k );
%         N_k  = get_TensProd(n_k,alpha_k)+sqrt(2/3)*m_k;
% 
%         % Compute mapping rules
%         [ alpha_Theta_Bk, alpha_Theta_Dk, hb_k ] = mapping_rules( S_k,alpha_k, p_k, m_k, Mc, Me, kb, kb_ex, kd, kd_ex, psi, ho );
%         
%         Fn_k    = get_TensProd(F_k,n_k);
% 
%         if ( Fn_k >0)
%             A_k     = Ao*(1+Fn_k);
%         else
%             A_k     = Ao;
%         end
%         
%         % Get variable Dk
%         D_k = A_k*get_TensProd((alpha_Theta_Dk-alpha_k),n_k);
%         
%         % Get isotropic and kinematic hardening variable
%         mm_k  = Cm*(1+eo)*D_k;
%         Kp1_k = get_TensProd(hb_k,n_k);
%                 
%         % Get plastic  multiplier
%         D_PlM = f_k / (-N_k*D_k*K_k + 2*G_k + p_k*(Kp1_k+sqrt(2/3)*mm_k));
%         
%         % Update variables
%         S_k = S_k - 2*G_k*D_PlM*n_k;
%         p_k = p_k - K_k*D_PlM*D_k;
%         
%         alpha_k = alpha_k + D_PlM*hb_k;
%         m_k     = m_k + D_PlM*mm_k;
%         
%         if (D_k<0)
%             F_k = F_k - D_PlM*Cf*abs(D_k)*(Fmax*n_k+F_k);
%         end
%         
%         Tol = f_k;
% 
%         %% To update
%         % S_k, p_k, alpha_k, m_k, F_k
%         
% % %         Sk     = S_o + 2*Gk*(Ddev_T-DVol_p);
% % %         alphak = alpha_o + DLam*Tk;
% % %         rk     = get_dev(Sk-alphak);
%         
%               
%     end
%     
% end
% 
% 
% 
% 
% 
% 
% end
% 

