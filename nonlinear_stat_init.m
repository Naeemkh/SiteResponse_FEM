function [C,K_mat,nl_e_n] = nonlinear_stat_init(C,K_mat,element_index)

% make the stiffness of the nonlinear elements as zero.

   nonlinear_elements_stat = element_index(:,12);
   j = 0;
   for i=1:size(nonlinear_elements_stat,1)
       
       if(nonlinear_elements_stat(i,1)==1)
            
               K_mat(i,i)     = 0;  
               K_mat(i,i+1)   = 0;
               K_mat(i+1,i)   = 0;  
               
               C(i,i)     = 0;  
               C(i,i+1)   = 0;
               C(i+1,i)   = 0;  
               
               j = j + 1;
           
       end
       
  
       
       
   end
   
   if (j~=0)
     K_mat(j,j)   = K_mat(j,j)./2;
   end
   
   nl_e_n = sum(nonlinear_elements_stat(:,1)==1);

end