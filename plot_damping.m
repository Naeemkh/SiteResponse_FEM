function plot_damping(output)

damping_model = output.simulationparams.damping_model;


if strcmp(damping_model,'Freq-Independent Rayleigh')==1

    dp = output.simulationparams.damping.FIRayleigh;
    w12 = output.simulationparams.damping.w12;
    
    figure;
    
    semilogx(dp(:,1),dp(:,2),'b')
    hold on
    semilogx(dp(:,1),dp(:,3),'r')
    semilogx(dp(:,1),dp(:,4),'k')
    legend('Mass term','Stiffness term','Mass and Stiffness')
    xlabel('Angular frequency')
    ylabel('Viscous Damping(%)')
    xlim([0 dp(end,1)])
    ylim([0 10])
    semilogx([w12(1) w12(1)],[0 15],'Color',[0.5 0.5 0.5],'LineStyle',':')
    semilogx([w12(2) w12(2)],[0 15],'Color',[0.5 0.5 0.5],'LineStyle',':')
   

end

end