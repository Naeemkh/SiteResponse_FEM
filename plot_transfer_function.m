function plot_transfer_function(filename_list,iteration_list,color_name,lw,flimit,smoothing_degree)

close(figure(1));


for i=1:size(filename_list,2)
    

output = filename_list{i};   
it_number = iteration_list{i};
color  = color_name{i};    
    
% Extracting the time series values. 

    
    F1=sprintf('%s%s','TF =output.iteration.it_',num2str(it_number),'.TF;');
    eval(F1);
    title1=sprintf('%s%s','Iteration : ',num2str(it_number),' Transfer function');



% Extracting the legend name
legend_1 = output.simulationparams.sim_name;
leg_1 = strrep(legend_1, '_', '-');
leg_2 = sprintf('%s%s%s',leg_1,' - Iteration: ',num2str(it_number));
legend_list{i} = leg_2;


% smoothed_data = smooth(TF(:,2));

% color_1 = color * 2;
% color_1(color_1>1)=1;
% 
% lw_1 =lw./2;


%% Plots
h=figure(1);
set(h,'position',[100 100 500 500])
plot(TF(:,1),TF(:,2),'color',color,'Linewidth',lw);
hold on
xlim(flimit);
%title('Displacement(cm)')
ylabel('Spectral Ratio')
xlabel('Frequency(Hz)')
% plot(TF(:,1),smoothed_data,'color',color,'Linewidth',lw);
hold on
legend(legend_list)

xlim(flimit)

end

end
    
