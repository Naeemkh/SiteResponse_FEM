function output = strain_damping_ggmax(output,eq_it)

userdepth = output.userdepth;
userdepth = userdepth';
element_index = output.element_index;

for i=1:size(userdepth,1)

    if userdepth(i,1) == 0
        temp_depth = 0;
    else
        temp_depth = max(element_index(element_index(:,4) < userdepth(i,1),4));
    end
    
    
     element_num = element_index(element_index(:,4)==temp_depth,1);
     
     
     f1=sprintf('%s%s%s','output.results.element_his.element_num_',num2str(element_num),'(eq_it,1)=eq_it;'); % Iteration
     f2=sprintf('%s%s%s','output.results.element_his.element_num_',num2str(element_num),'(eq_it,2)=output.results.effective_strain(element_num,1);'); % Effective strain
     f3=sprintf('%s%s%s','output.results.element_his.element_num_',num2str(element_num),'(eq_it,3)=element_index(element_num,7)./element_index(element_num,11);'); % G/Gmax
     f4=sprintf('%s%s%s','output.results.element_his.element_num_',num2str(element_num),'(eq_it,4)=element_index(element_num,8);'); % Damping
     
     eval(f1);
     eval(f2);
     eval(f3);
     eval(f4);

end

end