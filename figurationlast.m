function figurationlast(name)
close all
%%Figuration
ebsd = loadEBSD([name '.ang']);
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
figure(1);plot(ebsd1);
title('');
xlabel('לךל', 'FontSize', 14);
ylabel('לךל','FontSize', 14);
ci = get(ebsd1,'ci');
print('-dpng', [name 'ebsd.png']);
figure;
plot(ebsd1(ci > 0.1));
xlabel('לךל', 'FontSize', 14);
ylabel('לךל','FontSize', 14);
print('-dpng', [name 'ci.png']);
%%
%reconstruction
 load([name 'ebsd2.mat']);
 figure();
plot(ebsd2);
xlabel('לךל', 'FontSize', 14);
ylabel('לךל','FontSize', 14);
print('-dpng', [name 'ebsd2.png']);
%%
%%grad6
load([name 'N6H4.mat']);
f=get(ebsd2,'unitCell');
figure(); plot(f(:,1),f(:,2));
grid on;
xlabel('X, לךל','FontSize', 14); 
ylabel('Y, לךל','FontSize', 14);
print('-dpng', [name 'unitcell.png'])
prompt = 'נאחלונû ÿקויךט? ';
d = input(prompt);
figure();
hist(abs(H4*d),256);xlabel('Êנטגטחםא', 'FontSize', 14);
ylabel('׳טסכמ חםאקוםטי','FontSize', 14);
print('-dpng', [name 'hist6.png']);
prompt = '־דנאםטקוםטו הכÿ דנאםטצ? ';
z = input(prompt);

% figure();
% plot(ebsd2, 'property',abs(H4));
% xlabel('לךל', 'FontSize', 18);
% ylabel('לךל','FontSize', 18);
% colormap(1-gray);
H4m = abs(H4*d);
H4m(H4m > z) =  z/10;
figure();
plot(ebsd2, 'property',H4m);
xlabel('לךל', 'FontSize', 14);
ylabel('לךל','FontSize', 14);
colormap(1-gray);
colorbar;
print('-dpng', [name 'grad6.png']);
 
% figure();
% plot(ebsd2, 'property',abs(H4));
% xlabel('לךל', 'FontSize', 18);
% ylabel('לךל','FontSize', 18);
% colormap(1-gray);
H4m = abs(H4);
figure();
plot(ebsd2, 'property',H4m);
xlabel('לךל', 'FontSize', 14);
ylabel('לךל','FontSize', 14);
colormap(1-gray);
colorbar;
print('-dpng', [name 'grad0.png']);
% %%grad12
% load([name 'N12H4.mat']);
% % figure();
% % hist(abs(H4),256);
% % print('-dpng', [name 'hist12.png']);
% % z=input('');
% % figure();
% % plot(ebsd2, 'property',abs(H4));
% % xlabel('לךל', 'FontSize', 18);
% % ylabel('לךל','FontSize', 18);
% % colormap(1-gray);
% H4m = abs(H4);
% H4m(H4m > z) =  z;
% figure();
% plot(ebsd2, 'property',H4m);
% xlabel('לךל', 'FontSize', 14);
% ylabel('לךל','FontSize', 14);
% colormap(1-gray);
% print('-dpng', [name 'grad12.png']);
% H4m = abs(H4);
% H4m(H4m > z) =  5;
% figure();
% plot(ebsd2, 'property',H4m);
% xlabel('לךל', 'FontSize', 18);
% ylabel('לךל','FontSize', 18);
% colormap(1-gray);
%%
%%KAM
grains0 = calcGrains(ebsd2, 'threshold', 5*degree);
grains = calcGrains(grains0(grainSize(grains0) > 25), 'threshold', 2*degree);

% figure;plot(grains0);
% xlabel('לךל','FontSize', 14);
% ylabel('לךל','FontSize', 14);
% print('-dpng', [name 'grains0.png']);
% figure;plot(grains);
% %title('׀אחלונ חונםא במכüרו 25 ןטךסוכוי');
% xlabel('לךל','FontSize', 14);
% ylabel('לךל','FontSize', 14);
% print('-dpng', [name 'grains.png']);
KAM = calcKAM(grains);
%max(KAM);
p=max(KAM)/degree;
disp('max(KAM)/degree=');
disp(p);
KAM0 = KAM;
KAM0(KAM>5*degree) = 5*degree;
figure;plotKAM(grains);
figure;plot(get(grains,'ebsd'),'property',KAM0/degree);
colorbar;
xlabel('לךל','FontSize', 14);
ylabel('לךל','FontSize', 14);
colormap(1-gray);
print('-dpng', [name 'KAM.png']);
% % figure;plot(get(grains,'ebsd'),'property',KAM0);
% % colorbar;
% % xlabel('לךל');
% % ylabel('לךל');
% % grains0;
% % grains;
figure; plot(ebsd0, 'property', get(ebsd0, 'iq'));
xlabel('לךל','FontSize', 14);
ylabel('לךל','FontSize', 14);
colormap(gray);
print('-dpng', [name 'iq.png']);
end