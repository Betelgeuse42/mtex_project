function figuration(name1,name2,name3,name4)
close all
%%Figuration
ebsd = loadEBSD(name1);
% ebsd = loadEBSD([name1 'ebsd2.mat']);
ebsd0 = ebsd;
% Horizontal mirror 
ebsd0 = flipud(ebsd0); 

% X-axis direct to right
plotx2east;

% Plot EBSD
% plot(ebsd0);

%%
reg = [0 0 400 400];
ebsd1 = ebsd0(inpolygon(ebsd0, reg));
figure;plot(ebsd1);
title('');
xlabel('���', 'FontSize', 18);
ylabel('���','FontSize', 18);
ci = get(ebsd1,'ci');
figure;
plot(ebsd1(ci > 0.1));
xlabel('���', 'FontSize', 18);
ylabel('���','FontSize', 18);

%%
%reconstruction
load(name2);
figure();
plot(ebsd2);
xlabel('���', 'FontSize', 18);
ylabel('���','FontSize', 18);
%%
%%grad6
load(name3);
figure();
hist(abs(H4),256);
z=input('');
figure();
plot(ebsd2, 'property',abs(H4));
xlabel('���', 'FontSize', 18);
ylabel('���','FontSize', 18);
colormap(1-gray);
H4m = abs(H4);
H4m(H4m > z) =  z;
figure();
plot(ebsd2, 'property',H4m);
xlabel('���', 'FontSize', 18);
ylabel('���','FontSize', 18);
colormap(1-gray);��������
%%grad12
load(name4);
figure();
hist(abs(H4),256);
z=input('');
figure();
plot(ebsd2, 'property',abs(H4));
xlabel('���', 'FontSize', 18);
ylabel('���','FontSize', 18);
colormap(1-gray);
H4m = abs(H4);
H4m(H4m > z) =  z;
figure();
plot(ebsd2, 'property',H4m);
xlabel('���', 'FontSize', 18);
ylabel('���','FontSize', 18);
colormap(1-gray);
H4m = abs(H4);
H4m(H4m > z) =  5;
figure();
plot(ebsd2, 'property',H4m);
xlabel('���', 'FontSize', 18);
ylabel('���','FontSize', 18);
colormap(1-gray);
%%
%%KAM
grains0 = calcGrains(ebsd2, 'threshold', 5*degree);
grains = calcGrains(grains0(grainSize(grains0) > 25), 'threshold', 2*degree);

figure;plot(grains0);
xlabel('���');
ylabel('���');

figure;plot(grains);
%title('������ ����� ������ 25 ��������');
xlabel('���');
ylabel('���');

KAM = calcKAM(grains);
%max(KAM);
p=max(KAM)/degree;
disp('max(KAM)/degree=');
disp(p);
KAM0 = KAM;
KAM0(KAM>5*degree) = 5*degree;
figure;plot(get(grains,'ebsd'),'property',KAM0);
colorbar;
xlabel('���');
ylabel('���');

colormap(1-gray);
figure;plot(get(grains,'ebsd'),'property',KAM0);
colorbar;
xlabel('���');
ylabel('���');
grains0
grains
end