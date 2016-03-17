function plot_strain(filename_list,elementnumber_list,color_name,lw)


path1 = pwd;
close all


num_files = size(filename_list,2);


for i=1:size(filename_list,2)
    

output = filename_list{i};   
elementnumber = elementnumber_list{i};
    
    
% Extracting the time series values. 
% if strcmp(type,'absolute')==1
%     
%     F1=sprintf('%s%s','noderesults=output.node_',num2str(nodenumber),'.absolute.TDVA;');
%     eval(F1);
%     title1=sprintf('%s%s','Node Number : ',num2str(nodenumber),' - Absolute movement');
% elseif strcmp(type,'relative')==1
%     
%     F1=sprintf('%s%s','noderesults=output.node_',num2str(nodenumber),'.relative.TDVA;');
%     eval(F1);
%      title1=sprintf('%s%s','Node Number : ',num2str(nodenumber),' - Relative movement');
% else
%     
%     warning('Type is not defined. Absolute movement is presented.')
%     F1=sprintf('%s%s','noderesults=output.node_',num2str(nodenumber),'.absolute.TDVA;');
%     eval(F1);
%      title1=sprintf('%s%s','Node Number : ',num2str(nodenumber),' - Absolute movement');
%     
% end

strain_vec = output.results.strain_matrix(elementnumber,:);
dt         = output.simulationparams.dt;
time_vec   = (0:1:size(strain_vec,2)-1)*dt;

% Extracting the legend name
legend_1 = output.simulationparams.sim_name;
leg_1 = strrep(legend_1, '_', '-');
leg_2 = sprintf('%s%s%s',leg_1,'- Node : ', num2str(elementnumber));
legend_list{i} = leg_2;

color  = color_name{i};




%% Plots
figure(1);
subplot(2,1,1)
plot(time_vec,strain_vec,'color',color,'Linewidth',lw)
hold on
ylabel('Strain')
xlabel('Time (s)')
legend(legend_list)
subplot(2,1,2)
% Plot stress here.


end

end
    
