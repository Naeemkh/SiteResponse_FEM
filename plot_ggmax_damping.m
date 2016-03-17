function plot_ggmax_damping(output,element_no)


element_index = output.element_index;

element_type = element_index(element_no,9);

f1=sprintf('%s%s%s%s','ggmx = output.simulationparams.soil_pro.s_r_type_',num2str(element_type),'.ggmax;');
f2=sprintf('%s%s%s%s','damp = output.simulationparams.soil_pro.s_r_type_',num2str(element_type),'.damping;');

eval(f1);
eval(f2);

f3=sprintf('%s%s%s','damp_ggmx_his = output.results.element_his.element_num_',num2str(element_no),';');
eval(f3)


figure;
subplot(2,1,1)
semilogx(ggmx(:,1),ggmx(:,2))
hold on
scatter(damp_ggmx_his(:,2),damp_ggmx_his(:,3))
set(gca,'xscale','log')
ylabel('G/Gmax')
xlabel('Shear Strain')
title1=sprintf('%s%s','Element Number:',num2str(element_no),' Soil type: ',num2str(element_type));
title(title1)

subplot(2,1,2)
semilogx(damp(:,1),damp(:,2))
hold on
scatter(damp_ggmx_his(:,2),damp_ggmx_his(:,4)*100)
set(gca,'xscale','log')
ylabel('Damping(%)')
xlabel('Shear Strain')




end