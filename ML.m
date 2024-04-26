clear
Nr = 2;
Nt = 8;
% y_received_test = [-0.9045+1.4983j;  ...
%                     0.0040-1.7052j]; %load file from colab
y_received_test = cell2mat( struct2cell( load('y_received_test.mat') ) );

% s_sent_test = [-1.-1.j; ...
%                 -1.-1.j; ...
%                 -1.-1.j; ...
%                 -1.-1.j; ...
%                 -1.+1.j]; %load file from colab
s_sent_test = cell2mat( struct2cell( load('s_sent_test.mat') ) );

% H = [ 1.9113+0.7957j, -0.8540-0.5098j, -0.2446+0.2601j, -0.7002-0.9410j, 0.8046-0.4157j;...
%       0.7886+0.0291j,  0.3001+0.1989j,  0.2862+0.1658j,  0.1568+0.1546j, -0.3060-0.6859j]; %load file from colab
H = cell2mat( struct2cell( load('H.mat') ) );
n_classes = 4^Nt;  %PSK modulation
x = [1+1i, 1-1i, -1+1i, -1-1i];
codebook = permn(x, Nt).';
%%
n_errors = 0;
tic
H = H + 0.1*(randn(size(H)) + randn(size(H))*1i);
for i = 1:size(y_received_test,2)
    d = sum(abs(y_received_test(:,i) - H*codebook).^2);
    [M,I] = min(d);
    %disp(codebook(:,I))
    b_sent = ([real(s_sent_test(:,i)),imag(s_sent_test(:,i))]+1)/2;
    b_detected = ([real(codebook(:,I)),imag(codebook(:,I))]+1)/2;
    n_errors = n_errors + sum(sum(ne(b_sent,b_detected)));
    if i ==10000
        disp(b_sent)
        disp(b_detected)
    end
end
toc
BER = n_errors/(Nt*2*size(y_received_test,2))*100


