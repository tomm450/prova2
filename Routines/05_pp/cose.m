close all
it = [1:size(fvals_f,2)]; 

for i = 1:size(fvals_f,2)
    
    figure(i)
    plot((fvals_f{i})-(fvals_c{i}),'ro--');
    hold on
    plot((fvals_f{i})-(fvals_copt{i}),'bo--');
    grid on
    err_m(i,:) = [abs(mean(fvals_f{i}-fvals_c{i})), abs(mean(fvals_f{i}-fvals_copt{i}))]; 
    
    x_mat(i,:) = copt{i}{1};
    
    figure(10*i);
    plot(fvals_f{i},'rs');
    hold on
    plot(fvals_c{i},'bo');
    plot(fvals_copt{i},'g<');
    grid on
end


figure; plot(err_m(:,1),'ro-'); hold on; plot(err_m(:,2),'bo-');grid on

figure; plot(it,x_mat(:,1),'ro--',it,x_mat(:,2),'bs--',it,x_mat(:,3),'g<--');grid on