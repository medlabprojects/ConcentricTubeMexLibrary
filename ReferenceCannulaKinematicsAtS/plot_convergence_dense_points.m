err_p = dlmread('err_p_over_dense.txt','\n');
err_q = dlmread('err_q_over_dense.txt','\n');
pts = [0,3,5,10,20,30,50,100,150,175]+2;

loglog(pts,err_p(1:end-1)*1000);
hold on;
loglog(pts,err_q(1:end-1));
legend('p [mm]','q [rad]')
xlabel('Density [pts/mm]');
ylabel('Error [units in legend]');

figure
err_J_pa = dlmread('err_J_pa_over_dense.txt','\n');
loglog(pts,err_J_pa(1:end-1)*1000); hold on;

err_J_pt = dlmread('err_J_pt_over_dense.txt','\n');
loglog(pts,err_J_pt(1:end-1));

err_J_wa = dlmread('err_J_wa_over_dense.txt','\n');
loglog(pts,err_J_wa(1:end-1));

err_J_wt = dlmread('err_J_wt_over_dense.txt','\n');
loglog(pts,err_J_wt(1:end-1)/1000);

legend('J_{pa} [mm/rad]','J_{pt} [mm/mm]','J_{wa} [rad/rad]','J_{wt} [rad/mm]')
xlabel('Density [pts/mm]');
ylabel('Error [units in legend]');