function plot_compare(output1,output2,nodenumber)


F1=sprintf('%s%s','noderesults_1=output1.node_',num2str(nodenumber),'.TDVA;');
eval(F1);

F2=sprintf('%s%s','noderesults_2=output2.node_',num2str(nodenumber),'.TDVA;');
eval(F2);

% simulation name 

S1=sprintf('%s%s','leg_1=output1.simulationparams.sim_name;');
eval(S1);

S2=sprintf('%s%s','leg_2=output2.simulationparams.sim_name;');
eval(S2);


title1=sprintf('%s%s','Node Number : ',num2str(nodenumber));

figure;
subplot(311)
plot(noderesults_1(:,1),noderesults_1(:,2))
hold on 
plot(noderesults_2(:,1),noderesults_2(:,2),'r')
ylabel('Displacement')
title(title1);
legend(leg_1,leg_2)

subplot(312)
plot(noderesults_1(:,1),noderesults_1(:,3))
hold on
plot(noderesults_2(:,1),noderesults_2(:,3),'r')
ylabel('Velocity')

subplot(313)
plot(noderesults_1(:,1),noderesults_1(:,4))
hold on
plot(noderesults_2(:,1),noderesults_2(:,4),'r')
ylabel('Acceleration')
xlabel('Time (s)')

end