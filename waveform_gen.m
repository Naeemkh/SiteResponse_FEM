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


% Displacement matrix
% The force is implemented at node 10 from bottom. (end - 9)
% The results for force point node and 10 node ahead of that (end - 19) and
% The results at the surface and 10 points below the surface are reported. 
%
%      . 1
%      .
%      . 
%       
%      |    
%
%      . 20
%      .
%      . 10d
%      .
%      .

node_g1 = [1;element_index(end-9,2);element_index(end-19,2)];
node_g2 = element_det(:,2);


time = (0:dt:t-dt)';
output.defualtnodes = node_g1;
output.userdepth = depth;

node_g = sort(unique([node_g1;node_g2]));

output.simulationparams=sim_p;
output.nodetime = u;
output.element_index = element_index;

for ij=1:size(node_g,1)
    
   disp = u(node_g(ij,1),:)'; disp1 = [0; disp]; % Initial displacement is zero.
   vel  = diff(disp1)/dt; vel1=[0; vel];
   acc  = diff(vel1)/dt;
   
   F1=sprintf('%s%s%s','output.node_',num2str(node_g(ij,1)),'.TDVA=[time disp vel acc];');
   eval(F1)

end

output.simulationparams.totaltime=toc(t1);  % Recording simulation time (sec)

end