clear;
%Detectores ML y MMSE
NR = 4;             %number of receptors
N = 64*64;          %number of transmitters
% agrupacio = 1*1;
% agrupacio_lateral = sqrt(agrupacio);
% sectors = N/agrupacio;
% sectors_lateral = sqrt(sectors);
P = [20,30];        %power
repeticiones = 1e1; 
% Mat_R = zeros(NR,2^sectors);
% seq_recibida_ML = zeros(sectors, repeticiones);
% Errors_ML = zeros(sectors,repeticiones,200);
% seq_recibida_MSE = zeros(sectors, repeticiones);
% Errors_MSE = zeros(sectors, repeticiones,200);
% BER_MSE = zeros(length(P),1);
% BER_ML = zeros(length(P),1);
% M = zeros(N, sectors); cont = 0; cont2 = 0;

% for j = 1:sectors_lateral   %definicio M
%     for i = 1: agrupacio_lateral
%         for m = 1: sectors_lateral
%             for n = 1:agrupacio_lateral
%                 M(n+cont*agrupacio_lateral,m+cont2*sectors_lateral)=1;
%             end
%         cont = cont+1;
%         end
%     end
%     cont2 = cont2+1;
% end

%sequencia_candidatas = zeros(sectors,1);
sequencia_candidatas = [zeros(NR,N);ones(NR,N)];
%Mat_seq_cand = zeros(sectors, 2^sectors);
 
% cont =1;
% for n = 1:2^sectors %generacio candidates
%     Mat_seq_cand(:,n) = sequencia_candidatas(:);
%     if cont < 2^sectors
%         for m = 1:log2(cont)+1
%             sequencia_candidatas(m) = bitget(cont,m);
%         end
%     end
%     cont = cont+1;
% end
% 
% cont = 0;
%VecDis = 100*1.259.^[0:1:13];
%for Dist = 1:length(VecDis)
    Dist = 5;
    %valor_nom_corba = VecDis(Dist);
    %nom_corba = string(valor_nom_corba);
    %nom_corba = replace(nom_corba,'.',',');
    %loadname = strcat("",nom_corba,".mat"); %introduir nom del arxiu
    loadname = strcat("H_TX64x64_RX4_D5.mat"); %introduir nom del arxiu
    %Mat_Hs = load(loadname);
    %H = Mat_Hs.H;
    H = cell2mat(struct2cell(load(loadname)));
%%    
    for idx =1: length(P)
        %Variacio Potencia
        amp = 10^(P(idx)/20);
        %HM = H*M*amp;
        HM = H*amp;
        Mat_seq_rx = HM.*(2*sequencia_candidatas-1);
        cont = 1;
        %for rep =1:repeticiones
            %seq_enviada = randi([0 1], sectors,1);
            seq_enviada = randi([0 1], N,1); %generate a Nx1 uniform distributed sequence of 0 and 1
            sequencia = seq_enviada*2-1; %generate BPSK
            w = (randn(NR,1)+1i*randn(NR,1))/sqrt(2); %sqrt...?
            r = HM*sequencia+w; %received sequence
            %ML
            %Mat_R = repmat(r,1,2^sectors);
            Mat_R = repmat(r,1,N);
            ML = vecnorm(Mat_R - Mat_seq_rx).^2;
            [MLmod, MLrx] = min(ML);
            MLrx = min();
            seq_recibida_ML(:,cont) = (de2bi(MLrx-1,sectors)); %,false int2bit
            Errors_ML(:,cont,idx) = ne(seq_recibida_ML(:,cont),seq_enviada);
            %MSE
%             A = ((HM*HM')+eye(NR))\(HM);
%             MSE = A'*r;
%             seq_recibida_MSE(:,cont) = (sign(real(MSE))+1)/2;
%             Errors_MSE(:,cont,idx) = ne(seq_recibida_MSE(:,cont),seq_enviada);
%             cont = cont+1;
        %end
%     Total_Errors_ML = sum(Errors_ML(:,:,idx), 'all');
%     Total_Errors_MSE = sum(Errors_MSE(:,:,idx), 'all');
%     BER_ML(idx) = Total_Errors_ML/(sectors*repeticiones);
%     BER_MSE(idx) = Total_Errors_MSE/(sectors*repeticiones);
    end
%resname = strcat('Cruzados_A16_txd2_rxdRx12_P20,30_D',nom_corba);
%save(resname,"BER_ML","BER_MSE"); 
%end