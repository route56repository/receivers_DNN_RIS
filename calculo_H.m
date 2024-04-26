clear;
%% Parameters 

NR = 12;    % number of antennas in rx
NF = 64; %rows in RIS
NC = 64; %columns in RIS
f = 28e9;   % frequency of operation
c = 299792458; %speed of light
lambda = c/f;   %wavelength
kl = 2*pi/lambda; %wave number

N = NF*NC; %elements in RIS
d = 2.5*lambda; % patch to patch distance in RIS
dRx = 2*lambda; % antenna distance in RX
%% Coordinates
D = 30;
% RX coordinates matrix (RX elements)
    RX_center = [0 0 D];
    AntenasRx = zeros(NR,1,3); %NRx3 
    %place elements at (x,0,D)
    AntenasRx(:,1,3) = D; % same z coordinate
    AntenasRx(:,1,1) = (-(NR-1)/2*dRx:dRx:(NR-1)/2*dRx).'; 
    scatter3(AntenasRx(:,1,1), AntenasRx(:,1,2), AntenasRx(:,1,3))
    hold on

% RIS coordinates matrix (N elements)
    RIS_center = [0 0 0];
    RIS = zeros(NF, NC, 3); % NFxNCx3
    RIS(:,:,1) = repmat(-(NC-1)/2*d:d:(NC-1)/2*d,NF,1);
    RIS(:,:,2) = repmat((-(NF-1)/2*d:d:(NF-1)/2*d).',1,NC);
    scatter3(reshape(RIS(:,:,1),1,[]), reshape(RIS(:,:,2),1,[]), reshape(RIS(:,:,3),1,[])) % function rectangle to see full RIS
    xlabel('x') 
    ylabel('y')
    zlabel('z')
    hold off
%% Fields
    %Calcul matriu Distancies, Azimut i elevació de element i-essim respecte
    %recetptor i-essim de AntenasRx. 
    r = zeros(NR,N);
    azimuth = zeros(NR,N);
    elevation = zeros(NR,N);
    %faseAntena = zeros(NR, N);eTx = zeros(N*NR,3);
    %DirectivitatRx = zeros(NR, N);
    %faseAntenaRx = zeros(NR, N);
    %eRx = zeros(N*NR,3);
    Antenna = patchMicrostrip('Height' ,9.9931e-05,...
                            'Length',0.0048, ...
                            'Width',0.0062,...
                            'GroundPlaneLength', 0.0100, ...
                            'GroundPlaneWidth', 0.0100,... 
                            'PatchCenterOffset', [0 0],...
                            'FeedOffset', [0.0010 0],...
                            'Load', lumpedElement, ...
                            'Tilt',[90],...
                            'TiltAxis',[0 0 1] );
    antRx = dipole('Length', 0.0046967 , ...
                    'Width', 9.9931e-05, ...
                    'Tilt',90,...
                    'TiltAxis',[1 0 0]);
    %H = zeros(NR,N);
    H_near = zeros(NR,N);
    H_far = ones(NR,N);
    %cte = (4*pi)/lambda; 
%%
    RIS_vector = zeros(N,1,3);
    RIS_vector(:,1,1) = reshape(RIS(:,:,1),N,1);
    RIS_vector(:,1,2) = reshape(RIS(:,:,2),N,1);
    RIS_vector(:,1,3) = reshape(RIS(:,:,3),N,1);
    for m = 1:NR
        for n=1:N
            %X = ['Step --> ', 'RX: ', num2str(m) , ' RIS ', num2str(n)];
            %disp(X)
            vec_dif = reshape(AntenasRx(m,1,:)-RIS_vector(n,1,:), 3,1);
            [azimuth(m,n),elevation(m,n),r(m,n)] = cart2sph(vec_dif(1),vec_dif(2),vec_dif(3));
            %[eT, hTx] = EHfields(Antenna,f,vec_dif);
            %[eR, hRx] = EHfields(antRx,frx,(-vec_dif));
            fase2 = exp(1i*kl.*r(m,n));%en contrafase per treure una de les dos de la polaritzacio
            %H(m,n) = cte.*r(m,n).* fase2.*(eT.'*eR);
            H_near(m,n) = 1/r(m,n)*conj(fase2);
        end
    end
    H_far = 1/D*exp(-1i*kl.*D)*H_far;
    %filename = 'H_' + string(N) + 'x' + string(NR)+'_D' + string(D);
    %nom_corba = replace(nom_corba,'.',',');
    %save(filename,'H_0');
%end
 %S = svd(H);
% S=S/max(S);
% 
S_near = svd(H_near);
S_near=S_near/max(S_near);
S_far = svd(H_far);
S_far = S_far/max(S_far);
% 
[S_near S_far]

