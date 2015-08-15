i=(10:10:200);

iH33=zeros(1,length(i));
iH43=zeros(1,length(i));
iH32=zeros(1,length(i));
iH42=zeros(1,length(i));
iH31=zeros(1,length(i));
iH41=zeros(1,length(i));
for j=1:length(i)

load('3_very small areaebsd2.mat')
load('3_very small areaN6H4.mat');
load('3_very small areaT6.mat');
H3 = T(1:3:end,1)+T(2:3:end,2)+T(3:3:end,3);
[n, x] = hist(log10(abs(H4)+0.1),i(j));

iH33(1,j)=10^x(n == max(n))-0.1;
[n, x] = hist(log10(abs(H3)+0.1),i(j));

iH43(1,j)=10^x(n == max(n))-0.1;

load('2_very small areaebsd2.mat')
load('2_very small areaN6H4.mat');
load('2_very small areaT6.mat');

H3 = T(1:3:end,1)+T(2:3:end,2)+T(3:3:end,3);
[n, x] = hist(log10(abs(H4)+0.1),i(j));

iH32(1,j)=10^x(n == max(n))-0.1;
[n, x] = hist(log10(abs(H3)+0.1),i(j));

iH42(1,j)=10^x(n == max(n))-0.1;

load('1_small_areaebsd2.mat')
load('N6H41_small_area.mat');
load('T61_small_area.mat');

H3 = T(1:3:end,1)+T(2:3:end,2)+T(3:3:end,3);
[n, x] = hist(log10(abs(H4)+0.1),i(j));

iH31(1,j)=10^x(n == max(n))-0.1;
[n, x] = hist(log10(abs(H3)+0.1),i(j));

iH41(1,j)=10^x(n == max(n))-0.1;
end
disp('1_small area');
disp(iH31);
figure;plot(i,iH31);hold on;plot(i,iH41);grid on;hold off;
disp(iH41);
disp('2_very small area');
disp(iH32);
figure;plot(i,iH32);hold on;plot(i,iH42);grid on;hold off;
disp(iH42);
disp('3_very small area');
disp(iH33);
figure;plot(i,iH33);hold on;plot(i,iH43);grid on;hold off;
disp(iH43);
