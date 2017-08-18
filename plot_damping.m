function plot_damping(output,element_no,plot_Q_D)

damping_model = output.simulationparams.damping.use_damping;


if strcmp(damping_model,'RD2')==1

    dp = output.simulationparams.damping.FIRayleigh;
    w12 = output.simulationparams.damping.w12;
    
    dp(:,1)=dp(:,1)/(2*pi);
    w12 = w12 /(2*pi);
    
    figure(1);
    
    plot(dp(:,1),dp(:,2),'b')
    hold on
    plot(dp(:,1),dp(:,3),'r')
    plot(dp(:,1),dp(:,4),'k')
    legend('Mass term','Stiffness term','Mass and Stiffness')
    xlabel('Frequency (Hz)')
    ylabel('Viscous Damping(%)')
    xlim([0 dp(end,1)])
    ylim([0 10])
    plot([w12(1) w12(1)],[0 15],'Color',[0.5 0.5 0.5],'LineStyle',':')
    plot([w12(2) w12(2)],[0 15],'Color',[0.5 0.5 0.5],'LineStyle',':')
   

elseif strcmp(damping_model,'BKT2')==1
    
    figure(1);
    
     node_damping = output.element_index(element_no,8)/100;
     Q_0 = 1/(2*node_damping);
          
     f_max=10;
     w_max = 2*pi*f_max;
     
     f =  0 : 0.01 : f_max;
     w = f * 2 * pi;

     gamma_1 = 0.0373 * 2 * pi * f_max;
     gamma_2 = 0.3082 * 2 * pi * f_max;

     alpha_1 = (-2.656  * Q_0.^-0.8788 + 1.677)./Q_0;
     alpha_2 = (-0.5623 * Q_0.^-1.03   + 1.262)./Q_0;

     beta = (0.1876*Q_0.^(-0.9196)+0.6137)./(Q_0*2*pi*f_max);
    
     alpha_1_b = alpha_1 * Q_0;
     alpha_2_b = alpha_2 * Q_0;
     
     gamma_1_b = gamma_1 / w_max;
     gamma_2_b = gamma_2 / w_max;
     
     w_b = w ./ w_max ;
     
     beta_b = beta * w_max * Q_0;
     
     alpha_gamma_1 = alpha_1_b*gamma_1_b./(gamma_1_b.^2+w_b.^2) +...
                     alpha_2_b*gamma_2_b./(gamma_2_b.^2+w_b.^2) ;
               
     alpha_gamma_2 = alpha_1_b*(gamma_1_b.^2)./(gamma_1_b.^2+w_b.^2) +...
                     alpha_2_b*(gamma_2_b.^2)./(gamma_2_b.^2+w_b.^2) ;            
     
     Q_inv = Q_0 * w_b.*((beta_b + alpha_gamma_1)./(1-(Q_0.^-1).*(alpha_gamma_2)));
     
       
     if plot_Q_D == 'Q'
        
         plot(f,Q_0,'r')
         hold on
         plot(f,Q_inv);
         xlabel('Frequency (Hz)')
         ylabel('Q')
         
         
     elseif plot_Q_D == 'D'
         
         plot(f,(1./(2*Q_0))*100,'r')
         hold on
         plot(f,(0.5*Q_inv)*100);
         xlabel('Frequency (Hz)')
         ylabel('Damping (%)')
         
     end
    
end

end