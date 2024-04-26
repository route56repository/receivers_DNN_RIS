function [BER] = ML_detection(Nr,Nt, y_received, H, s_sent)
    n_classes = 4^Nt;  %PSK modulation
    x = [1+1i, 1-1i, -1+1i, -1-1i];
    codebook = permn(x, Nt).';
    n_errors = 0;
    tic
    %H = H + 0.1*(randn(size(H)) + randn(size(H))*1i);
    for i = 1:size(y_received,2)
        d = sum(abs(y_received(:,i) - H*codebook).^2);
        [M,I] = min(d);
        %disp(codebook(:,I))
        b_sent = ([real(s_sent(:,i)),imag(s_sent(:,i))]+1)/2;
        b_detected = ([real(codebook(:,I)),imag(codebook(:,I))]+1)/2;
        n_errors = n_errors + sum(sum(ne(b_sent,b_detected)));
    %     if i ==10000
    %         disp(b_sent)
    %         disp(b_detected)
    %     end
    end
    toc
    BER = n_errors/(Nt*2*size(y_received,2));
end

