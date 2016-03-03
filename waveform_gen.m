function output = waveform_gen(sim_p,u,element_index,depth,t1)

% Last update : Feb 17 2016
% Query the elements according to the depth.
% Element numbers starts from the surface.
% in all elements I'm reporting the upper node value.
% The code is not interpolating, it is reporting the closest node.

dt=sim_p.dt;
t =sim_p.sim_time;

element_det = zeros(size(depth,2),3);

for i=1:size(element_det)
    
    element_det(i,1) = depth(1,i);
    element_det(i,2) = element_index(depth(1,i)-element_index(:,4)>=0 & element_index(:,5)-depth(1,i)>0,2);
    element_det(i,3) = element_det(i,2)+1;
    
end


node_g1 = [1;element_index(end,3);element_index(end-19,2)];
node_g2 = element_det(:,2);


time = (0:dt:t-dt)';
output.defualtnodes = node_g1;
output.userdepth = depth;

node_g = sort(unique([node_g1;node_g2]));

output.simulationparams=sim_p;
output.nodetime = u;
output.element_index = element_index;
output.simulationparams.n_element=size(element_index,1);
output.simulationparams.n_nodes=size(element_index,1)+1;
output.simulationparams.n_timest=t/dt;
% Displacement, Velocity, and Acceleration of the base


disp_base = u(end,:)'; disp1 = [0; disp_base]; % Initial displacement is zero.
vel_base  = diff(disp1)/dt; vel1=[0; vel_base];
acc_base  = diff(vel1)/dt;




for ij=1:size(node_g,1)
    
    disp_ab = u(node_g(ij,1),:)'; disp1 = [0; disp_ab]; % Initial displacement is zero.
    vel_ab  = diff(disp1)/dt; vel1=[0; vel_ab];
    acc_ab  = diff(vel1)/dt;
    
    disp_rel = disp_ab - disp_base;
    vel_rel  = vel_ab  - vel_base;
    acc_rel  = acc_ab  - acc_base;
    
    F1=sprintf('%s%s%s','output.node_',num2str(node_g(ij,1)),'.absolute.TDVA=[time disp_ab vel_ab acc_ab];');
    F2=sprintf('%s%s%s','output.node_',num2str(node_g(ij,1)),'.relative.TDVA=[time disp_rel vel_rel acc_rel];');
    eval(F1)
    eval(F2)
end

output.simulationparams.totaltime=toc(t1);  % Recording simulation time (sec)

end