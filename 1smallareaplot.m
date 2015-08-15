name='1_small_area';
close all;
ebsd = loadEBSD([name '.ang']);
figure('Position', [100 100 800 600]);
plot(ebsd);
xlabel('мкм','FontSize', 14);
ylabel('мкм','FontSize', 14);
hold on; 
rectangle( 'Position',  [0 70 20 20], 'lineWidth', 3,'EdgeColor','b'); 
hold off;
print('-dpng', [name 'ebsd.png']);
ebsd2 = ebsd(inpolygon(ebsd, [0 70 20 20]));
figure('Position', [100 100 800 600]);
plot(ebsd2);
xlabel('мкм','FontSize', 14);
ylabel('мкм','FontSize', 14);
print('-dpng', [name 'ebsd2.png']);
grains = calcGrains(ebsd2, 'threhold', 5*degree);
grains1 = calcGrains(grains(grainSize(grains) > 10), 'threhold', 5*degree);
figure('Position', [100 100 800 600]);
plot(ebsd2, 'property', 'iq');
xlabel('мкм','FontSize', 14);
ylabel('мкм','FontSize', 14);
colormap(gray);
print('-dpng', [name 'iq.png']);
figure('Position', [100 100 800 600]);
plotKAM(grains1, 'secondorder');
xlabel('мкм','FontSize', 14);
ylabel('мкм','FontSize', 14);
colormap(1-gray);
colorbar;
hold on; 
plotBoundary(grains1, 'property', [0 15],'lineWidth', 3, 'linecolor', 'b');
hold off;
xlabel('мкм','FontSize', 14);
ylabel('мкм','FontSize', 14);
pause(10);
print('-dpng', [name 'KAM_grains.png']);
%% Filter


%% Get data
X = get(ebsd2, 'X');
Y = get(ebsd2, 'Y');
ot = get(ebsd2,'orientation');
CS = symmetry('m-3m');
%%Array of radius calculation
%%Variables 
dx=zeros(size(X));
dy=zeros(size(X));
k=length(X);
H3=zeros(size(X));
H4=zeros(size(X));
at=k*3;
T=zeros(at,3);
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
    T(i*3-2,1)=L4(1,1);
    T(i*3-1,1)=L4(2,1);
    T(i*3,1)=L4(3,1);
    T(i*3-2,2)=L4(1,2);
    T(i*3-1,2)=L4(2,2);
    T(i*3,2)=L4(3,2);
    T(i*3-2,3)=L4(1,3);
    T(i*3-1,3)=L4(2,3);
    T(i*3,3)=L4(3,3);
    Vx=[L4(1,1) L4(1,2) L4(1,3)];
     Vy=[L4(2,1) L4(2,2) L4(2,3)];
      Vz=[L4(3,1) L4(3,2) L4(3,3)];
    H4(i)=max([sqrt(sum(Vx.*Vx,2)),sqrt(sum(Vy.*Vy,2)),sqrt(sum(Vz.*Vz,2))],[],2);
    H3(i)=L4(1,1)+L4(2,2)+L4(3,3);
if (mod(i,fix(k/25)) == 0)
    fprintf('|')
end
end;
fprintf('Done')

t = toc;
disp('t=');
disp(t);


%%
% H4m = abs(H4);
% H4m(H4m >  round(mean(H4m)/2)) =  round(mean(H4m)/2);
% c=(find(H4m==14));
% cc=size(c)/k;
% disp('Процент значений кривизны,которые срезаются');
% disp(cc(1, 1));
% figure; plot(ebsd2, 'property',abs( H4m));
% xlabel('мкм');
% ylabel('мкм');
% colormap(1-colormap(gray));
% figure; plot(ebsd2, 'property',abs( H4m));
% xlabel('мкм');
% ylabel('мкм');
% colormap(gray);
save(['N6H4' name '.mat'], 'H4', 'reg');
save(['T6' name '.mat'], 'T', 'reg');
%%



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
h=12;
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
T(i*3-2,1)=L4(1,1);
    T(i*3-1,1)=L4(2,1);
    T(i*3,1)=L4(3,1);
    T(i*3-2,2)=L4(1,2);
    T(i*3-1,2)=L4(2,2);
    T(i*3,2)=L4(3,2);
    T(i*3-2,3)=L4(1,3);
    T(i*3-1,3)=L4(2,3);
    T(i*3,3)=L4(3,3);
    Vx=[L4(1,1) L4(1,2) L4(1,3)];
     Vy=[L4(2,1) L4(2,2) L4(2,3)];
      Vz=[L4(3,1) L4(3,2) L4(3,3)];
        H4(i)=max([sqrt(sum(Vx.*Vx,2)),sqrt(sum(Vy.*Vy,2)),sqrt(sum(Vz.*Vz,2))],[],2);
    H3(i)=L4(1,1)+L4(2,2)+L4(3,3);
if (mod(i,fix(k/25)) == 0)
    fprintf('|')
end
end;
fprintf('Done')

t = toc;
disp('t=');
disp(t);


%%
% H4m = abs(H4);
% H4m(H4m >  round(mean(H4m)/2)) =  round(mean(H4m)/2);
% c=(find(H4m==14));
% cc=size(c)/k;
% disp('Процент значений кривизны,которые срезаются');
% disp(cc(1, 1));
% figure; plot(ebsd2, 'property',abs( H4m));
% xlabel('мкм');
% ylabel('мкм');
% colormap(1-colormap(gray));
% figure; plot(ebsd2, 'property',abs( H4m));
% xlabel('мкм');
% ylabel('мкм');
% colormap(gray);
save(['N12H4' name '.mat'], 'H4', 'reg');
save(['T12' name '.mat'], 'T', 'reg');
load(['N6H4'name'.mat']);
f=get(ebsd2,'unitCell');
figure(); plot(f(:,1),f(:,2));
grid on;
xlabel('X, мкм','FontSize', 14); 
ylabel('Y, мкм','FontSize', 14);
print('-dpng', [name 'unitcell.png'])
figure();
hist(log10(H4+0.1),128));xlabel('Логарифм Кривизны', 'FontSize', 14);
ylabel('Число значений','FontSize', 14);
print('-dpng', [name 'hist6.png']);
prompt = 'Ограничение для границ? ';
z = input(prompt);

% figure();
% plot(ebsd2, 'property',abs(H4));
% xlabel('мкм', 'FontSize', 18);
% ylabel('мкм','FontSize', 18);
% colormap(1-gray);
H4m = abs(H4);
H4m(H4m > z) =  z/10;
figure();
plot(ebsd2, 'property',H4m);
xlabel('мкм', 'FontSize', 14);
ylabel('мкм','FontSize', 14);

 colormap(1-gray)
 print('-dpng', [name 'gradmaxwithoutgrains.png']);
figure();
plot(ebsd2, 'property',H4);
xlabel('мкм', 'FontSize', 14);
ylabel('мкм','FontSize', 14);colormap(1-gray)
print('-dpng', [name 'grad6.png']);
figure();
plot(ebsd2, 'property',log10(H4+0.1));
xlabel('мкм', 'FontSize', 14);
ylabel('мкм','FontSize', 14);colormap(1-gray)
colormap(1-gray);
colorbar;
print('-dpng', [name 'gradmaxlog.png']);
 
