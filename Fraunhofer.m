clear;
%% Parameters 

NR = 3;    % number of antennas in rx
NF = 64; %rows in RIS
NC = 64; %columns in RIS
f = 28e9;   % frequency of operation
frx = 28e9; % ?
c = 299792458; %speed of light
lambda = c/f;   %wavelength
kl = 2*pi/lambda; %wave number

N = NF*NC; %elements in RIS
d = 2.5*lambda; % patch to patch distance in RIS
dRx = 7.5*lambda;
d_RIS_RX = logspace(1,4,600);
eigen_near = zeros(size(d_RIS_RX,2),NR);
eigen_far = zeros(size(d_RIS_RX,2),NR);

%% Coordinates
for i = 1:size(d_RIS_RX,2)
    [H_far,H_near, H_near_sim, H_far_sim] = get_H(4, d_RIS_RX(i));
    S_near = svd(H_near_sim);
    %S_near=S_near/max(S_near);
    eigen_near(i,:) = S_near';
    S_far = svd(H_far_sim);
    %S_near=S_near/max(S_near);
    eigen_far(i,:) = S_far';
end
%%
FD = 2*(sqrt(2)*63*d)^2/lambda
figure
%subplot(2,1,1)
%plot(d_RIS/lambda, plot_FA)
%xlabel('Element to Element distance in the RIS (\lambda)')
%ylabel('Fraunhofer distance')
%title('Fraunhofer distance for a ' + string(NF) + 'x' + string(NC) + ' RIS and ' + string(NR)+ ' receiving elements separated d_{RX} = 7.5\lambda')
%subplot(2,1,2)
loglog(d_RIS_RX, eigen_near(:,1), 'b', 'LineWidth',2)
hold on
plot(d_RIS_RX, eigen_near(:,2), 'r', 'LineWidth',2)
hold on
plot(d_RIS_RX, eigen_near(:,3), 'g', 'LineWidth',2)
grid on
xlabel('Distance from RIS to RX (m)')
ylabel('Channel eigenvalues')
%legend('\lambda_0', '\lambda_1', '\lambda_2')
%figure
loglog(d_RIS_RX, eigen_far(:,1), 'b--', 'LineWidth',2)
% hold on
% loglog(d_RIS_RX, eigen_far(:,2), 'r--', 'LineWidth',2)
% hold on
% loglog(d_RIS_RX, eigen_far(:,3), 'g--', 'LineWidth',2)
% hold on
grid on
xlabel('Distance from RIS to RX (m)')
ylabel('Channel eigenvalues')
set(gca,'fontsize', 20)
%title('H eigenvalues at D_{RIS-RX} = '+ string(D)+' m')
legend('\lambda_{0, NF}', '\lambda_{1, NF}', '\lambda_{2, NF}', '\lambda_{0, FF}')
%ylim([0,2])