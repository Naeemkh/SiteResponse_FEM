function plot_wave(filename_list,iteration_list,nodenumber_list,type,color_name,lw,tlimit,flimit,filter_it,filter_type,corner_freq,cpdf)

close(figure(1));


for i=1:size(filename_list,2)
    

output = filename_list{i};   
nodenumber = nodenumber_list{i};
iteration = iteration_list{i};    
    
% Extracting the time series values. 
if strcmp(type,'absolute')==1
    
    F1=sprintf('%s%s%s%s%s','noderesults=output.iteration.it_',num2str(iteration),'.node_',num2str(nodenumber),'.absolute.TDVA;');
    eval(F1);
    title1=sprintf('%s%s','Node Number : ',num2str(nodenumber),' - Absolute movement');
elseif strcmp(type,'relative')==1
    
    F1=sprintf('%s%s%s%s%s','noderesults=output.iteration.it_',num2str(iteration),'.node_',num2str(nodenumber),'.relative.TDVA;');
    eval(F1);
     title1=sprintf('%s%s','Node Number : ',num2str(nodenumber),' - Relative movement');
else
    
    warning('Type is not defined. Absolute movement is presented.')
    F1=sprintf('%s%s%s%s%s','noderesults=output.iteration.it_',num2str(iteration),'.node_',num2str(nodenumber),'.absolute.TDVA;');
    eval(F1);
     title1=sprintf('%s%s','Node Number : ',num2str(nodenumber),' - Absolute movement');
    
end

sim_time=output.simulationparams.sim_time;
 


% Extracting the legend name
legend_1 = output.simulationparams.sim_name;
leg_1 = strrep(legend_1, '_', '-');
leg_2 = sprintf('%s%s%s',leg_1,'- Node : ', num2str(nodenumber));
legend_list{i} = leg_2;



t_wave = noderesults;
color  = color_name{i};

dt = t_wave(2,1) - t_wave(1,1);

% pick 75 s of data
NPick = floor(sim_time / dt);
t_wave = t_wave(1:NPick,:);

N2          = size(t_wave,1);
NP2         = 2^(nextpow2(N2));
% NP2     = 8192;
FS      = 1/dt;
f       = (0:1/(NP2-1):1)*FS;
fc  = 1; % converting to cm 
ftype = filter_type;

% Filter defination
if strcmp(filter_it,'yes')==1
    
    
if strcmp(ftype,'bandpass')==1
        
    F_F_2 = corner_freq(2);
    F_F_1 = corner_freq(1);
    
    % low pass (elliptic)
    [zl,pl,kl] = ellip(70,0.0001,35,F_F_2/(FS/2),'low');
    [sosl,gl] = zp2sos(zl,pl,kl);	               % Convert to SOS form
    Hdl = dfilt.df2tsos(sosl,gl);
    % High pass (elliptic)
    [zh,ph,kh] = ellip(15,0.001,40,F_F_1/(FS/2),'high');
    [sosh,gh] = zp2sos(zh,ph,kh);	               % Convert to SOS form
    Hdh = dfilt.df2tsos(sosh,gh);                  % Create a dfilt object
    
    
    % filter the waves
    for iif=2:10
    F_wave_l=filtfilt(sosl,gl,t_wave(:,iif)); 
    F_wave_hl=filtfilt(sosh,gh,F_wave_l);
    t_wave(:,iif) = F_wave_hl;
    end
   
elseif   strcmp(ftype,'lowpass')==1
    
    F_F_1 = corner_freq(1);
    
    % low pass (elliptic)
%     [zl,pl,kl] = ellip(70,0.0001,35,F_F_1/(FS/2),'low');
%     [sosl,gl] = zp2sos(zl,pl,kl);	               % Convert to SOS form
%     Hdl = dfilt.df2tsos(sosl,gl);

 % low pass (elliptic)
    [zl,pl,kl] = butter(4,F_F_1/(FS/2),'low');
    [sosl,gl] = zp2sos(zl,pl,kl);	               % Convert to SOS form
    Hdl = dfilt.df2tsos(sosl,gl);

    
    for iif=2:4
    F_wave_l = filtfilt(sosl,gl,t_wave(:,iif));
    t_wave(:,iif) = F_wave_l;
    end
    
elseif   strcmp(ftype,'highpass')==1
    
    F_F_1 = corner_freq(1);
    
    % High pass (elliptic)
    [zh,ph,kh] = ellip(15,0.001,40,F_F_1/(FS/2),'high');
    [sosh,gh] = zp2sos(zh,ph,kh);	               % Convert to SOS form
    Hdh = dfilt.df2tsos(sosh,gh);                  % Create a dfilt object
    
    for iif=2:10
    F_wave_h = filtfilt(sosh,gh,t_wave(:,iif));
    t_wave(:,iif) = F_wave_h;
    end
    
    
end
end

f_wave_d  = fft(t_wave(:,2)*fc,NP2)/(N2); % Displacement
f_wave_v  = fft(t_wave(:,3)*fc,NP2)/(N2); % Velocity
f_wave_a  = fft(t_wave(:,4)*fc,NP2)/(N2); % Acceleration



%% Plots
h=figure(1);
set(h,'position',[100 100 1600 1200])
subplot(3,3,[1 2])
plot(t_wave(:,1),t_wave(:,2)*fc,'color',color,'Linewidth',lw)
hold on
xlim(tlimit);
%title('Displacement(cm)')
ylabel('Disp')
legend(legend_list)
subplot(3,3,3)
plot(f,abs(f_wave_d),'color',color,'Linewidth',lw)
hold on
xlim(flimit)
title('FFT - Amplitude')


subplot(3,3,[4 5])
plot(t_wave(:,1),t_wave(:,3)*fc,'color',color,'Linewidth',lw)
hold on
xlim(tlimit);
ylabel('Vel')
subplot(3,3,6)
plot(f,abs(f_wave_v),'color',color,'Linewidth',lw)
hold on
xlim(flimit)


subplot(3,3,[7 8])
plot(t_wave(:,1),t_wave(:,4)*fc,'color',color,'Linewidth',lw)
hold on
xlim(tlimit);
ylabel('Acc')
xlabel('Time (s)')
subplot(3,3,9)
plot(f,abs(f_wave_a),'color',color,'Linewidth',lw)
hold on
xlabel('Frequency (Hz)')
xlim(flimit)

% if cpdf == 1
% savename1=sprintf('%s%s',path1,'/Disp.pdf');
% saveas(gcf,savename1)
% end

end

end
    
%plot_wave({output,output,output,output},{1,2,3,4},{1,1,1,1},'absolute',{'b','r','g','c'},1,[0 10],[0 3],'no','bandpass',[0.1 3],'no')
