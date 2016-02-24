function wave_snapshot(output,snaptime,dx)

dt = output.simulationparams.dt;
sim_time = output.simulationparams.sim_time;
sim_name = output.simulationparams.sim_name;

snaptime1=snaptime(snaptime<=sim_time & 0<=snaptime);

if numel(snaptime1) ~= numel(snaptime)
    
    warning('Out of range requested times have been deleted.')
    
end


snaptime=sort(snaptime1);



u  = output.nodetime;

element_index = output.element_index;




  FigHandle = figure;
  set(FigHandle, 'Position', [100, 100, 800, 1800]);

for i=1:numel(snaptime)
    if i > 4
        
        warning('Only first 4 requested times are plotted out.');
        break;
        
    else
        subplot(2,2,i)
        plot(u(:,floor(snaptime(i)/dt)),[element_index(:,4);element_index(end,5)] ,'linewidth',1)
        set(gca,'Ydir','reverse')
        time_step = sprintf('%s%s',num2str(snaptime(i)),' s');
        text(0.8*dx, 0.1*element_index(end,5) , time_step)
        xlim([-dx dx])
        set(gca,'XMinorGrid','on')
        
        if i==1
            sim_name_1 = strrep(sim_name, '_', '-');
            title(sim_name_1);
            xlabel('Amplitude')
            ylabel('Depth (m)')
            
        end
    end
    
    
end

sim_name_1


end