function wave_propagation_animation(output,responsetype)

figure

dt = output.simulationparams.dt;
u  = output.nodetime;
sim_time = output.simulationparams.sim_time;
nt_step = sim_time/dt;
element_index = output.element_index;
time_factor = floor(1/dt * 0.01);


if strcmp(responsetype,'disp')==1
    
    for ij =1: nt_step/time_factor;
        plot(u(:,ij*time_factor),[element_index(:,4);element_index(end,5)] ,'linewidth',1)
        set(gca,'Ydir','reverse')
        time_step = sprintf('%2.2f%s',(ij*dt)*time_factor,' s');
        text(1, 10, time_step)
        xlim([-0.5 0.5])
        pause(0.01)
        
        % file name
        %     file_name=sprintf('%s%04d%s','wave_snapshot/wave_',ij,'.png');
        %
        %     saveas(gcf,file_name,'png')
        
        %     close all
    end
    
elseif strcmp(responsetype,'vel')==1
    
    disp_ab = u; disp1 = [zeros(size(u,1))  disp_ab]; % Initial displacement is zero.
    vel_ab  = diff(disp1,1,2)/dt; 
    
    
    for ij =1: nt_step/time_factor;
        plot(vel_ab(:,ij*time_factor),[element_index(:,4);element_index(end,5)] ,'linewidth',1)
        set(gca,'Ydir','reverse')
        time_step = sprintf('%2.2f%s',(ij*dt)*time_factor,' s');
        text(1, 10, time_step)
        xlim([-10 10])
        pause(0.01)
        
        % file name
        %     file_name=sprintf('%s%04d%s','wave_snapshot/wave_',ij,'.png');
        %
        %     saveas(gcf,file_name,'png')
        
        %     close all
    end
    
    

    
    
elseif strcmp(responsetype,'acc')==1
    disp_ab = u; disp1 = [zeros(size(u,1))  disp_ab]; % Initial displacement is zero.
    vel_ab  = diff(disp1,1,2)/dt; vel1=[zeros(size(u,1)) vel_ab];
    acc_ab  = diff(vel1,1,2)/dt;
    
    
    
    

    
    for ij =1: nt_step/time_factor;
        plot(acc_ab(:,ij*time_factor),[element_index(:,4);element_index(end,5)] ,'linewidth',1)
        set(gca,'Ydir','reverse')
        time_step = sprintf('%2.2f%s',(ij*dt)*time_factor,' s');
        text(0.05, 10, time_step)
        xlim([-0.001 0.001])
         pause(0.01)
        

        
        
        
        
        
        
        % file name
        %     file_name=sprintf('%s%04d%s','wave_snapshot/wave_',ij,'.png');
        %
        %     saveas(gcf,file_name,'png')
        
        %     close all
    end
    
    
    
else
    warning('Please pick acc, vel, disp as a particle motion.')
end

end