clear
diary  on;
diary log.txt;
%% Parameters
Nr = 3; % Receptores
Nt = 64*64; % Transmisores

n_sector = 4; % numero efectivo de simbolos transmitidos
SNR = 1e3; % SNR en lineal
%d_RIS_RX = 20;
max_epochs_near = 300;
max_epochs_far = 300;
%% Iterate distances
d_RIS_RX = [50];
%d_RIS_RX = logspace(1,4,150);
nn_MSE_near = zeros(size(d_RIS_RX));
nn_MSE_far = zeros(size(d_RIS_RX));
MMSE_near = zeros(size(d_RIS_RX));
MMSE_far = zeros(size(d_RIS_RX));
BER_ML_near = zeros(size(d_RIS_RX));
BER_ML_far = zeros(size(d_RIS_RX));
num_promedios = 10;
nn_MSE_near_p = zeros(num_promedios, size(d_RIS_RX,2));
BER_ML_near_p = zeros(num_promedios, size(d_RIS_RX,2));
BER_ML_far_p = zeros(num_promedios, size(d_RIS_RX,2));
%%
for p = 1:num_promedios
    for k = 1:size(d_RIS_RX, 2)
        %% Channel and Training/Test Dataset 
        disp(d_RIS_RX(k))
        [H_far,H_near, H_near_sim, H_far_sim] = get_H(n_sector, d_RIS_RX(k)); % debe depender de n_sector, SNR
        M = 5e3;  % numero de accesos al canal para entrenamiento, 
                    % el numero de simbolos enviados en total es MxNt
        r = binornd(1,0.5, [n_sector M])*2 -1;
        i = binornd(1,0.5, [n_sector M])*2 -1;
        s = r +1i*i;
        n = 1/sqrt(SNR)*(1/sqrt(2)*(rand(Nr, M) +1i*rand(Nr,M)));
        y_near = H_near_sim*s + n;
        y_far = H_far_sim*s + n;

        %% Channel estimation
        H_est_near = y_near*s'*inv(s*s');
        H_est_far = y_far*s'*inv(s*s');

        %% Preprocessing for NN (including CSI to received samples)
        y_p_near = H_est_near'*inv(H_est_near*H_est_near' + 1/SNR*eye(size(H_est_near*H_est_near')))*y_near;
        y_p_far = H_est_far'*inv(H_est_far*H_est_far' + 1/SNR*eye(size(H_est_far*H_est_far')))*y_far;

        %% NN
        %net_near = fitnet([2*n_sector], 'trainbr');
         net_near = fitnet([8,8], 'trainbr');
%         net_near.trainParam.epochs = max_epochs_near;
%         net_near.divideParam.trainRatio = 80/100;
%         net_near.divideParam.valRatio = 0/100;
%         net_near.divideParam.testRatio = 20/100;
         net_near = train(net_near,[real(y_p_near); imag(y_p_near)],[real(s); imag(s)]);
        view(net_near)

    %     net_near = fitnet([2*n_sector], 'trainbr');
    %     net_near.trainParam.epochs = 800;
    %     net_near = train(net_near,[real(y_near); imag(y_near)],[real(s); imag(s)]);

        %net_far = fitnet([2*n_sector], 'trainbr');
%         net_far = fitnet([8 8], 'trainbr');
%         net_far.trainParam.epochs = max_epochs_far;
%         net_far.divideParam.trainRatio = 80/100;
%         net_far.divideParam.valRatio = 0/100;
%         net_far.divideParam.testRatio = 20/100;
%         net_far = train(net_far,[real(y_p_far); imag(y_p_far)],[real(s); imag(s)]);

        %% Generate new samples for detection
        M = 2e2;  % numero de accesos al canal,
        r = binornd(1,0.5, [n_sector M])*2 -1;
        i = binornd(1,0.5, [n_sector M])*2 -1;
        s = r +1i*i;
        n = 1/sqrt(SNR)*(1/sqrt(2)*(rand(Nr, M) +1i*rand(Nr,M)));

        y_near = H_near_sim*s + n;
        y_p_near = H_est_near'*inv(H_est_near*H_est_near' + 1/SNR*eye(size(H_est_near*H_est_near')))*y_near;

        y_far = H_far_sim*s + n;
        y_p_far = H_est_far'*inv(H_est_far*H_est_far' + 1/SNR*eye(size(H_est_far*H_est_far')))*y_far;
        %% ML Detection
        tic
        BER_ML_near(k) = ML_detection(Nr,n_sector, y_near, H_est_near, s);
        toc
        BER_ML_far(k) = ML_detection(Nr,n_sector, y_far, H_est_far, s);

        %% MMSE Detection
    %     MMSE_near(k) = sum(abs((y_p_near - s)).^2, 'all')/(M*n_sector);
    %     MMSE_far(k) = sum(abs((y_p_far - s)).^2, 'all')/(M*n_sector);
        %% NN Detection
%         tic
%         output_nn_near = net_near([real(y_p_near); imag(y_p_near)]);
%         toc
%         nn_MSE_near(k) = sum(abs((output_nn_near - [real(s); imag(s)])).^2, 'all')/(M*n_sector);
%         
%         
%     %     
%         output_nn_far = net_far([real(y_p_far); imag(y_p_far)]);
%         nn_MSE_far(k) = sum(abs((output_nn_far - [real(s); imag(s)])).^2, 'all')/(M*n_sector);
        %output_b_nn = zeros(size(output_nn));
        %output_b_nn(output_nn>=0) = 1;
        %output_b_nn(output_nn<0) = -1;

        %BER_NN = sum(sum(ne([real(s); imag(s)],output_b_nn)))/(M*n_sector)

        disp(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss'))
    end
%     nn_MSE_near_p(p, :) =nn_MSE_near;
%     nn_MSE_far_p(p, :) =nn_MSE_far;
    BER_ML_near_p(p,:) = BER_ML_near;
    BER_ML_far_p(p,:) = BER_ML_far;
end
%% NN PLOT
% figure
% loglog(d_RIS_RX, sum(nn_MSE_near_p,1)/num_promedios)
% hold on
% loglog(d_RIS_RX, sum(nn_MSE_far_p,1)/num_promedios)
% %title('NN MSE for Near field')
% xlabel('distance [m]')
% ylabel('MSE')
% grid on
% legend('MSE_{NF}', 'MSE_{FF}')
% %set(gca, 'YScale', 'log')
% %set(gca, 'XScale', 'log')
% savefig('nn_MSE.fig');

%% MMSE PLOT
% figure
% loglog(d_RIS_RX, MMSE_near, 'LineWidth',2)
% hold on
% loglog(d_RIS_RX, MMSE_far, 'LineWidth',2)
% %title('NN MSE for Near field')
% xlabel('distance [m]')
% ylabel('MSE')
% %ylim([1e-2 1e1])
% set(gca,'fontsize', 20)
% grid on
% legend('MSE_{MMSE, NF}', 'MSE_{MMSE, FF}')

%% ML PLOT
figure
loglog(d_RIS_RX, sum(BER_ML_near_p,1)/num_promedios)
hold on
loglog(d_RIS_RX, sum(BER_ML_far_p,1)/num_promedios)
%title('NN MSE for Near field')
xlabel('distance [m]')
ylabel('BER')
ylim([1e-2 1.1e0])
grid on
legend('BER_{ML, NF}', 'BER_{ML, FF}')
set(gca,'fontsize', 20)
savefig('BER_ML.fig');
%%
% figure
% data = output_nn_near(1:4,:) + 1i*output_nn_near(5:8,:);
% unos = data(real(data)>0 & imag(data)>0);
% plot(real(unos), imag(unos), 'gx')
% hold on 
% unos = data(real(data)>0 & imag(data)<0);
% plot(real(unos), imag(unos), 'yx')
% hold on 
% unos = data(real(data)<0 & imag(data)>0);
% plot(real(unos), imag(unos), 'rx')
% hold on 
% unos = data(real(data)<0 & imag(data)<0);
% plot(real(unos), imag(unos), 'bx')
% yline(0)
% xline(0)
%%
% loglog(d_RIS_RX, nn_MSE_near_p(1,:))
% hold on
% loglog(d_RIS_RX, nn_MSE_near_p(2,:))
% hold on
% loglog(d_RIS_RX, nn_MSE_near_p(3,:))