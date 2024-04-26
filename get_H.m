function [H_far,H_near, H_near_sim, H_far_sim] = get_H(n_sector, D)
    %% Parameters 
    f = 28e9;   % frequency of operation
    c = 299792458; %speed of light
    lambda = c/f;   %wavelength
    kl = 2*pi/lambda; %wave number
    
    %RX
    NF_rx = 1;
    NC_rx = 3;
    NR = NF_rx*NC_rx;    % number of antennas in rx
    dRx = 7.5*lambda; % antenna distance in RX
    
    %RIS
    NF = 64; %rows in RIS
    NC = 64; %columns in RIS
    N = NF*NC; %elements in RIS
    d = 2.5*lambda; % patch to patch distance in RIS
    
    %% Coordinates
    % RX coordinates matrix (RX elements)
        RX_center = [0 0 D];
        AntenasRx = zeros(NF_rx, NC_rx, 3); % NFxNCx3 

        %place elements
        AntenasRx(:,:,3) = D; % same z coordinate
        AntenasRx(:,:,1) = repmat(-(NC_rx-1)/2*dRx:dRx:(NC_rx-1)/2*dRx,NF_rx,1); %x coordinate
        AntenasRx(:,:,2) = repmat((-(NF_rx-1)/2*dRx:dRx:(NF_rx-1)/2*dRx).',1,NC_rx); 
        
        if n_sector == 16
            AntenasRx(:,:,3) = D + AntenasRx(:,:,2).^2 + AntenasRx(:,:,1).^2;
        else 
            % Exception
            AntenasRx(1,1,2) = AntenasRx(1,1,2) - 3*lambda;
            AntenasRx(1,2,2) = AntenasRx(1,2,2);
            AntenasRx(1,3,2) = AntenasRx(1,3,2) + 3*lambda;
        end
        
        figure(1)
        scatter3(reshape(AntenasRx(:,:,1),1,[]), reshape(AntenasRx(:,:,2),1,[]), reshape(AntenasRx(:,:,3),1,[]))
        xlabel('x (m)')
        ylabel('y (m)')
        set(gca,'fontsize', 20)
%         hold on

    % RIS coordinates matrix (N elements)
        RIS_center = [0 0 0];
        RIS = zeros(NF, NC, 3); % NFxNCx3
        RIS(:,:,1) = repmat(-(NC-1)/2*d:d:(NC-1)/2*d,NF,1);
        RIS(:,:,2) = repmat((-(NF-1)/2*d:d:(NF-1)/2*d).',1,NC);
%         figure(2)
%         scatter3(reshape(RIS(:,:,1),1,[]), reshape(RIS(:,:,2),1,[]), reshape(RIS(:,:,3),1,[])) % function rectangle to see full RIS
%         xlabel('x') 
%         ylabel('y')
%         zlabel('z')
%         hold off
    %% Fields
        %Calcul matriu Distancies, Azimut i elevació de element i-essim respecte
        %recetptor i-essim de AntenasRx. 
        r = zeros(NR,N);
        azimuth = zeros(NR,N);
        elevation = zeros(NR,N);
        H_near = zeros(NR,N);
        H_far = ones(NR,N);
        %cte = (4*pi)/lambda; 
    %%
        RIS_vector = zeros(N,1,3);
        AntenasRx_vector = zeros(NR,1,3);
        RIS_vector(:,1,1) = reshape(RIS(:,:,1),N,1);
        RIS_vector(:,1,2) = reshape(RIS(:,:,2),N,1);
        RIS_vector(:,1,3) = reshape(RIS(:,:,3),N,1);
        AntenasRx_vector(:,1,1) = reshape(AntenasRx(:,:,1),NR,1);
        AntenasRx_vector(:,1,2) = reshape(AntenasRx(:,:,2),NR,1);
        AntenasRx_vector(:,1,3) = reshape(AntenasRx(:,:,3),NR,1);
        AntenasRx_vector_farfield = AntenasRx_vector;
        AntenasRx_vector_farfield(:,1,3) = 100000;
        
        for m = 1:NR
            for n=1:N
                %X = ['Step --> ', 'RX: ', num2str(m) , ' RIS ', num2str(n)];
                %disp(X)
                vec_dif = reshape(AntenasRx_vector(m,1,:)-RIS_vector(n,1,:), 3,1);
                vec_dif_farfield = reshape(AntenasRx_vector_farfield(m,1,:)-RIS_vector(n,1,:), 3,1);
                %near field
                [azimuth(m,n),elevation(m,n),r(m,n)] = cart2sph(vec_dif(1),vec_dif(2),vec_dif(3));
                fase2 = exp(1i*kl.*r(m,n));
                H_near(m,n) = 1/r(m,n)*conj(fase2);
                
                %far field
                [azimuth(m,n),elevation(m,n),r(m,n)] = cart2sph(vec_dif_farfield(1),vec_dif_farfield(2),vec_dif_farfield(3));
                fase2 = exp(1i*kl.*r(m,n));
                H_far(m,n) = 1/r(m,n)*conj(fase2);
            end
        end
        H_far = H_far*100000/D;
        
    % simplificacion de la matriz
    div_f = sqrt(n_sector);
    div_c = sqrt(n_sector);
    
    H_near_sim = zeros(NR, n_sector);
    H_far_sim = zeros(NR, n_sector);
    count = 1;
    figure
    for j = 1:div_c
        for i = 1:div_f
            selection_matrix = zeros(64,64);
            selection_matrix((i-1)*NF/div_f+1:i*NF/div_f,(j-1)*NC/div_c+1:j*NC/div_c) = 1;
            selection_matrix = reshape(selection_matrix, [1 N]);
            
            %Plot sectors
            x = RIS_vector(:,1,1);
            y = RIS_vector(:,1,2);
            x = x(selection_matrix'==1);
            y = y(selection_matrix'==1);
            if j ==1 && i ==1
                scatter(x,y, [],'red')
                hold on
            elseif j ==1 && i ==2
                scatter(x,y, [],'blue')
                hold on
            elseif j ==2 && i ==1
                scatter(x,y, [],'green')
                hold on
            elseif j ==2 && i ==2
                scatter(x,y, [],'yellow')
                hold on
            end
            xlabel('x (m)')
            ylabel('y (m)')
            set(gca,'fontsize', 20) 
            hold on
            
            
            selection_matrix = repelem(selection_matrix,NR,1);
            H_near_sim(:,count) = sum(H_near.*selection_matrix, 2);
            H_far_sim(:,count) = sum(H_far.*selection_matrix, 2);
            count = count + 1;
        end
    end

    S_near_sim = svd(H_near_sim);
    S_near_sim=S_near_sim/max(S_near_sim);
    disp(['Near field eigenvalues: ' string(S_near_sim')])
end

