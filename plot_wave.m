function plot_wave(output,nodenumber,type)

if strcmp(type,'absolute')==1
    
    F1=sprintf('%s%s','noderesults=output.node_',num2str(nodenumber),'.absolute.TDVA;');
    eval(F1);
    title1=sprintf('%s%s','Node Number : ',num2str(nodenumber),' - Absolute movement');
elseif strcmp(type,'relative')==1
    
    F1=sprintf('%s%s','noderesults=output.node_',num2str(nodenumber),'.relative.TDVA;');
    eval(F1);
     title1=sprintf('%s%s','Node Number : ',num2str(nodenumber),' - Relative movement');
else
    
    warning('Type is not defined. Absolute movement is presented.')
    F1=sprintf('%s%s','noderesults=output.node_',num2str(nodenumber),'.absolute.TDVA;');
    eval(F1);
     title1=sprintf('%s%s','Node Number : ',num2str(nodenumber),' - Absolute movement');
    
end



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