function plot_wave(output,nodenumber)


F1=sprintf('%s%s','noderesults=output.node_',num2str(nodenumber),'.TDVA;');
eval(F1);

title1=sprintf('%s%s','Node Number : ',num2str(nodenumber));

figure;
subplot(311)
plot(noderesults(:,1),noderesults(:,2))
ylabel('Displacement')
title(title1);

subplot(312)
plot(noderesults(:,1),noderesults(:,3))
ylabel('Velocity')

subplot(313)
plot(noderesults(:,1),noderesults(:,4))
ylabel('Acceleration')
xlabel('Time (s)')


end