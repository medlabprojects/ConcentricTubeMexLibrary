data = dlmread('tmech_design_convergence_data.txt');

%%
%n_normal_points = data(:,1);
%dense_density = data(:,2);
L1 = data(:,1);
L2 = data(:,2);
L3 = data(:,3);
Lt1 = data(:,4);
Lt2 = data(:,5);
Lt3 = data(:,6);
k1 = data(:,7);
k2 = data(:,8);
k3 = data(:,9);
psi1 = data(:,10);
psi2 = data(:,11);
psi3 = data(:,12);
beta1 = data(:,13);
beta2 = data(:,14);
beta3 = data(:,15);
p_err = data(:,16);
q_err = data(:,17);
psi_err = data(:,18);
mz_err = data(:,19);
J_va_err = data(:,20);
J_vt_err = data(:,21);
J_wa_err = data(:,22);
J_wt_err = data(:,23);
kb1 = data(:,24);
kb2 = data(:,25);
kb3 = data(:,26);
ct1 = data(:,27);
ct2 = data(:,28);
ct3 = data(:,29);
stability = data(:,30);
q_angle_err = data(:,31);
q_mag_err = data(:,32);
psi_diff_mag = data(:,33);
mz_mag = data(:,34);
%q_err_x = data(:,37);
%q_err_y = data(:,38);
%q_err_z = data(:,39);
normal_density = data(:,35);
normal_length = data(:,36);
Fx = data(:,37);
Fy = data(:,38);
Fz = data(:,39);
tip_deflection = data(:,40);
normal_points_used = data(:,41);

Fmag = sqrt(Fx.^2 + Fy.^2 + Fz.^2);

%% compute the power
P = zeros(size(data, 1), 1);
for i=1:length(P)
   PP = get_section_info( [beta1(i),beta2(i),beta3(i)]',[Lt1(i),Lt2(i),Lt3(i)]',[L1(i),L2(i),L3(i)]',[k1(i),k2(i),k3(i)]');
   P(i) = PP;
end

%%
%plot convergence
L = length(p_err);
jj = 1:100:L;
err_up_to = @(k)(mean(q_err(1:k)));
p_err_convergence = arrayfun(err_up_to, jj);
loglog(jj,p_err_convergence)



%%
i = find(n_normal_points >= 0 & dense_density/1000 >= 0 );
loglog(q_angle_err(i), p_err(i)*1000, 'r.', 'MarkerSize',0.001)
%xlim([10^-4,10]);
%ylim([10^-25,2*10^-4]);

%%
theta1 = psi2(1001:2000) - psi1(1001:2000);
theta2 = psi3(1001:2000) - psi1(1001:2000);

p2 = p_err(1001:2000);
n2 = n_normal_points(1001:2000);
dd2 = dense_density(1001:2000);
i = find(abs(n2-50) < 10 );
loglog(dense_density(i), p2(i),'r.')

%%
i = find(Fmag <10);
ExposedLength = L1 + beta1;
semilogy(tip_deflection(i)./ExposedLength(i),p_err(i)./ExposedLength(i),'b.')
xlim([0,2])
ylim([10^-18,10^0])