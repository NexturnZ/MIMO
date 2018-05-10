clear all;  close all;
dbstop if error;
%% QPKS
datasize = 300000;
EbNo = 0:2:20;
M = 4;
x = randsrc(2,datasize/log2(M),0:M-1);
x1 = qammod(x,M);
h = randn(4,datasize/log2(M))+1i*randn(4,datasize/log2(M));
h = h./sqrt(2);

for idx = 1:length(EbNo)
   sigma1 = sqrt(1/(4*10.^(EbNo(idx)/10))); 
   n = sigma1*(randn(2,datasize/log2(M))+1i*randn(2,datasize/log2(M)));
   y = x1+n;
   y1 = x1+n./h(1:2,:);
   x2 = qamdemod(y,M);
   x3 = qamdemod(y1,M);
   sigma2 = sqrt(1/(2*10.^(EbNo(idx)/10)));
   n = sigma2*(randn(4,datasize/log2(M))+1i*randn(4,datasize/log2(M)));
   n1(1,:) = (conj(h(1,:)).*n(1,:)+h(2,:).*conj(n(2,:)))./(sum(abs(h(1:2,:)).^2));
   n1(2,:) = (conj(h(2,:)).*n(1,:)-h(1,:).*conj(n(2,:)))./(sum(abs(h(1:2,:)).^2));
   y = x1+n1;
   x4 = qamdemod(y,M);
   n2(1,:) = (conj(h(1,:)).*n(1,:)+h(2,:).*conj(n(2,:))+...
             conj(h(3,:)).*n(3,:)+h(4,:).*conj(n(4,:)))./(sum(abs(h).^2));
   n2(2,:) = (conj(h(2,:)).*n(1,:)-h(1,:).*conj(n(2,:))+...
             conj(h(4,:)).*n(3,:)-h(3,:).*conj(n(4,:)))./(sum(abs(h).^2));
   
    y1 = x1+n2;
    x5 = qamdemod(y1,M);
   
   [~, ber1(idx)] = biterr(x,x2,log2(M));
   [~, ber2(idx)] = biterr(x,x3,log2(M));
   [~, ber3(idx)] = biterr(x,x4,log2(M));
   [~, ber4(idx)] = biterr(x,x5,log2(M));
end
figure;
semilogy(EbNo, ber1,'-r*',EbNo, ber2, '-go', EbNo, ber3, '-bd', EbNo, ber4, '-k.');
title('QPSK');

%% 16-QAM
datasize = 300000;
EbNo = 0:2:20;
M = 16;
x = randsrc(2,datasize/log2(M),0:M-1);
x1 = qammod(x,M);
h = randn(4,datasize/log2(M))+1i*randn(4,datasize/log2(M));
h = h./sqrt(2);
ber = zeros(1,length(EbNo));

for idx = 1:length(EbNo)
   sigma1 = sqrt(1/(4*10.^(EbNo(idx)/10))); 
   n = sigma1*(randn(2,datasize/log2(M))+1i*randn(2,datasize/log2(M)));
   y = x1+n;
   y1 = x1+n./h(1:2,:);
   x2 = qamdemod(y,M);
   x3 = qamdemod(y1,M);
   sigma2 = sqrt(1/(2*10.^(EbNo(idx)/10)));
   n = sigma2*(randn(4,datasize/log2(M))+1i*randn(4,datasize/log2(M)));
   n3(1,:) = (conj(h(1,:)).*n(1,:)+h(2,:).*conj(n(2,:)))./(sum(abs(h(1:2,:)).^2));
   n3(2,:) = (conj(h(2,:)).*n(1,:)-h(1,:).*conj(n(2,:)))./(sum(abs(h(1:2,:)).^2));
   y = x1+n3;
   x4 = qamdemod(y,M);
   n4(1,:) = (conj(h(1,:)).*n(1,:)+h(2,:).*conj(n(2,:))+...
             conj(h(3,:)).*n(3,:)+h(4,:).*conj(n(4,:)))./(sum(abs(h).^2));
   n4(2,:) = (conj(h(2,:)).*n(1,:)-h(1,:).*conj(n(2,:))+...
             conj(h(4,:)).*n(3,:)-h(3,:).*conj(n(4,:)))./(sum(abs(h).^2));
   
    y1 = x1+n4;
    x5 = qamdemod(y1,M);
   
   [~, ber1(idx)] = biterr(x,x2,log2(M));
   [~, ber2(idx)] = biterr(x,x3,log2(M));
   [~, ber3(idx)] = biterr(x,x4,log2(M));
   [~, ber4(idx)] = biterr(x,x5,log2(M));
end
figure; 
semilogy(EbNo, ber1,'-r*',EbNo, ber2, '-go', EbNo, ber3, '-bd', EbNo, ber4, '-k.');
title('16-QAM');

