%%full program
clear all
close all

% ebsd - loaded data
% ebsd0 - full data
% ebsd1 - working data
% ebsd2 - modified working data
name = 'testEBSD001';
%% First loading
% % Load EBSD data
%ebsd = loadEBSD('testEBSD001.ang');
% save('data0.mat', 'ebsd');

%% Second loading
load('testEBSD001.txt');

%% EBSD representation properties
ebsd0 = ebsd;
% Horizontal mirror 
ebsd0 = flipud(ebsd0); 

% X-axis direct to right
plotx2east;

% Plot EBSD
% plot(ebsd0);

%% Cut data
reg = [0 0 400 400];
ebsd1 = ebsd0(inpolygon(ebsd0, reg));
figure;plot(ebsd1);
title('');
xlabel('мкм');
ylabel('мкм');
colorbar;

%% Get data
X = get(ebsd1, 'X');
Y = get(ebsd1, 'Y');
ot1 = get(ebsd1,'orientation');

%% Filter
ci = get(ebsd1,'ci');
figure;
plot(ebsd1(ci > 0.1));
xlabel('мкм');
ylabel('мкм');
colorbar;
%% Preparation

% Number of elements in circle
h = 6;

%% Reconstruct 4+
work = 1;

while work > 0
    [ot1, ci, work] = recon(ot1, ci, X, Y, h, 4, 0);
end

ebsd2 = set(ebsd1,'rotations',ot1);
% ebsd2 = set(ebsd2,'CI',ci);
figure; plot(ebsd2(ci > 0.1));
xlabel('мкм');
ylabel('мкм');


%% Reconstruct 3+
work = 1;

while work > 0
    [ot1, ci, work] = recon(ot1, ci, X, Y, h, 3, 1);
end
% ebsd2 = set(ebsd2,'CI',ci);
ebsd2 = set(ebsd1,'rotations',ot1);
figure; plot(ebsd2(ci > 0.1));
xlabel('мкм');
ylabel('мкм');


ebsd2 = set(ebsd1,'rotations',ot1);
figure;plot(ebsd2(ci > 0.1));
xlabel('мкм');
ylabel('мкм');

%%
%% Filter
ci = get(ebsd2,'ci');
% plot(ebsd(ci > 0.1));
save(['ebsd2' name '.mat'], 'ebsd2', 'reg');

%% Get data
X = get(ebsd2, 'X');
Y = get(ebsd2, 'Y');
ot = get(ebsd2,'orientation');
CS = symmetry('m-3m');
%%Array of radius calculation
%%Variables 
Radius=zeros(size(X));
dx=zeros(size(X));
dy=zeros(size(X));
k=length(X);
H1=zeros(size(X));
H2=zeros(size(X));
H3=zeros(size(X));
H4=zeros(size(X));
H5=zeros(size(X));
%number of elements in circle
h=6;
tic
%% I is under study, J -other
for i=1:k
    dx(1:k)=(X(1:k)-X(i));
    dy(1:k)=(Y(1:k)-Y(i));
    Radius=(dy.^2+dx.^2).^0.5; 
   
  
     %just in case
    Radius(i)=30;%больше чем область
  
    
    
    % search for the minimum
    z= min(Radius);
    %tolerance 20%
    if h==12
        eps=1.8;
    else if h==6
            eps=1.3;
        end;
    end;
  
    ind = find(Radius < z*eps);
    if numel(ind)<h
        continue;
    else
        if numel(ind)>h
            disp(ind);
            break;
        end;
        % any((i - offset) < 1) || any((i - offset) > length(X))
        %Top coordinate matrix
        Rtop=zeros(3,h+1);
        %Bottom coordinate matrix
        Rbottom=zeros(3,h+1);
        %Orientation Matrix Rodrigues
        U1=zeros(3,h+1);
        %Orientation Matrix Euler
        U2=zeros(3,h+1);
        %Orientation Matrix axes
        U3=zeros(3,h+1);
        U4=zeros(3,h+1);
          if norm(Rtop)==0 
                if norm(Rbottom)==0
                   for c=1:h
                    g=ind(c,1); 
                    %coordinate matrix
                    Rtop(1,(c+0))=dx(g);
                    Rtop(2,(c+0))=dy(g);
                    Rtop(3,(c+0))=1;
                    Rbottom(1,(c+0))=dx(g);
                    Rbottom(2,(c+0))=dy(g);
                    Rbottom(3,(c+0))=-1;
                   end;
                    Rtop(1,c+1)=0;
                    Rtop(2,c+1)=0;
                    Rtop(3,c+1)=1;
                    Rbottom(1,c+1)=0;
                    Rbottom(2,c+1)=0;
                    Rbottom(3,c+1)=-1;
               else
                continue
               end;
          end;
%           for 
              c = 1:h;
                    %disorientation 
                    m = inverse(ot(ind(c)))*ot(i);
%                     %rodrigues vector
%                     v= Rodrigues(m);
%                     U1(1,(c+0))=getx(v);
%                     U1(2,(c+0))=gety(v);
%                     U1(3,(c+0))=getz(v);
%                     %euler angle
%                     u= Euler(m)/degree;
%                     U2(1,(c+0))=u(1,1);
%                     U2(2,(c+0))=u(1,2);
%                     U2(3,(c+0))=u(1,3);
                    %axis and angle
                    [m1,a] = project2FundamentalRegion(m);
                    a = a/degree;
                    
                    b = axis(rotation(m1));
                    b = Miller(b,CS,CS);

%                     b = axis(m);
                    b0 = [getx(b),gety(b),getz(b)];
%                     b1 = b0/norm(b0);
                    b1 = b0;
                    
                    u1_b = repmat(a,1,3).*b1;
                    U3(:,c) = u1_b';
%                     U3(1,(c+0))=u1_b(1);
%                     U3(2,(c+0))=u1_b(2);
%                     U3(3,(c+0))=u1_b(3);            
%           end;
                
            U3(1,c(end)+1)=0;
            U3(2,c(end)+1)=0;
            U3(3,c(end)+1)=0;
             U4 = U3;

     end; 
    
   

%     [~,L1]=gradientmatrix([Rtop Rbottom],[U1 U1]);
%       H1(i)=trace(L1);
%     [~,L2]=gradientmatrix([Rtop Rbottom],[U2 U2]);
%    H2(i)=trace(L2);
%    [~,L3]=gradientmatrix([Rtop Rbottom],[U3 U3]);
   [~,L4]=gradientmatrix([Rtop Rbottom],[U4 U4]);
%     H3(i)=trace(L3);
    H4(i)=trace(L4);

if (mod(i,fix(k/25)) == 0)
    fprintf('|')
end
end;
fprintf('Done')

t = toc;
disp('t=');
disp(t);


%%
H4m = abs(H4);
H4m(H4m >  round(mean(H4m)/2)) =  round(mean(H4m)/2);
c=(find(H4m==14));
cc=size(c)/k;
disp('Процент значений кривизны,которые срезаются');
disp(cc(1, 1));
figure; plot(ebsd2, 'property',abs( H4m));
xlabel('мкм');
ylabel('мкм');
colormap(1-colormap(gray));

save(['6H4' name '.mat'], 'H4', 'reg');
%save(['12H4' name '.mat'], 'H4', 'reg');
%% grains & KAM
grains0 = calcGrains(ebsd2, 'threshold', 5*degree);
grains = calcGrains(grains0(grainSize(grains0) > 25), 'threshold', 2*degree);

figure;plot(grains0);
xlabel('мкм');
ylabel('мкм');

figure;plot(grains);
%title('Размер зерна больше 25 пикселей');
xlabel('мкм');
ylabel('мкм');

KAM = calcKAM(grains);
%max(KAM);
p=max(KAM)/degree;
disp('max(KAM)/degree=');
disp(p);
KAM0 = KAM;
KAM0(KAM>5*degree) = 5*degree;
figure;plot(get(grains,'ebsd'),'property',KAM0);
colorbar;
xlabel('мкм');
ylabel('мкм');

colormap(1-gray);
figure;plot(get(grains,'ebsd'),'property',KAM0);
colorbar;
xlabel('мкм');
ylabel('мкм');
grains0
grains

