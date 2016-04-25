

FEM_rel = output.node_1.relative.TDVA(:,[1,4]);
FEM_abs = output.node_1.absolute.TDVA(:,[1,4]);

load simulation_results/deepsoil_ricker_fd.mat
load simulation_results/deepsoil_ricker_td.mat
% input_m = load('input_acc/sin_2_5_hz_8sec.txt');
seismosoil_TD_incident = load('simulation_results/ricker_5Hz_accel_on_surface_TD_incident.txt');
seismosoil_TD_outcrop  = load('simulation_results/ricker_5Hz_accel_on_surface_TD_outcrop.txt');
seismosoil_FD_incident = load('simulation_results/ricker_5Hz_accel_on_surface_FD_incident.txt');
seismosoil_FD_outcrop  = load('simulation_results/ricker_5Hz_accel_on_surface_FD_outcrop.txt');

time_shift=0.34;
dt=output.simulationparams.dt;
x=0:dt:FEM_rel(end,1)+time_shift;
y_rel = [zeros(time_shift/dt,1);FEM_rel(:,2)];
y_abs = [zeros(time_shift/dt,1);FEM_abs(:,2)];

FEM_t_rel = [x' y_rel];
FEM_t_abs = [x' y_abs];


figure;
hold on
% plot(FEM_rel(:,1),FEM_rel(:,2)/9.81)
% plot(FEM_abs(:,1),FEM_abs(:,2)/9.81,'k')
% plot(FEM_abs(:,1),(FEM_abs(:,2)+FEM_rel(:,2))/(2*9.81),'m')
% plot(Deep(:,1),Deep(:,2),'g')
% plot(seismosoil_TD_incident(:,1),seismosoil_TD_incident(:,2)/9.81,'b')
% plot(seismosoil_FD_incident(:,1),seismosoil_FD_incident(:,2)/9.81,':b')
plot(seismosoil_TD_outcrop(:,1),seismosoil_TD_outcrop(:,2)/9.81,'r')
plot(seismosoil_FD_outcrop(:,1),seismosoil_FD_outcrop(:,2)/9.81,':r')
plot(FEM_rel(:,1),-FEM_rel(:,2)/9.81,'k')
plot(Deepsoil_ricker_fd(:,1),Deepsoil_ricker_fd(:,2),'b')
plot(Deepsoil_ricker_td(:,1),Deepsoil_ricker_td(:,2),'g')

legend('Seismosoil-TD','Seismosoil-FD','FEM-Naeem-TD','DeepSoil-FD','DeepSoil-TD')
xlabel('Time(s)')
ylabel('Acceleration(m/s^2)')
% 
% figure;
% hold on
% % plot(seismosoil_TD_incident(:,1),seismosoil_TD_incident(:,2)/9.81,'b')
% % plot(seismosoil_FD_incident(:,1),seismosoil_FD_incident(:,2)/9.81,':b')
% plot(seismosoil_TD_outcrop(:,1),seismosoil_TD_outcrop(:,2)/9.81,'r')
% plot(seismosoil_FD_outcrop(:,1),seismosoil_FD_outcrop(:,2)/9.81,':r')
% plot(FEM_t_rel(:,1),FEM_t_rel(:,2)/9.81,'k')
% % plot(FEM_t_abs(:,1),FEM_t_abs(:,2)/9.81,'m')
% plot(Deepsoil_ricker_fd(:,1),Deepsoil_ricker_fd(:,2),'b')
