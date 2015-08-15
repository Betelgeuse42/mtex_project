clear all
close all
ebsd = createTestSample('Isl2A.txt', 100, 100, 'sin', 'Amplitude', 3, 'Period', 6, 'Noise', 0.5, 'Coeff', [1 1; 1 1; 1 1]);

% ebsd0=ebsd; 
% figure(); plot(ebsd0);
% xlabel('לךל','FontSize',14);
% ylabel('לךל','FontSize',14);

%% Get data
X = get(ebsd, 'X');
Y = get(ebsd, 'Y');
ot = get(ebsd,'orientation');
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
h=12;
tic
%% I is under study, J -other
for i=1:k
    dx(1:k)=(X(1:k)-X(i));
    dy(1:k)=(Y(1:k)-Y(i));
    Radius=(dy.^2+dx.^2).^0.5; 
   
  
     %just in case
    Radius(i)=30;%במכüרו קול מבכאסעü
  
    
    
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
    H3(i)=norm(L4);

if (mod(i,fix(k/25)) == 0)
    fprintf('|')
end
end;
fprintf('Done')

t = toc;
disp('t=');
disp(t);


%%

figure; plot(ebsd, 'property',abs( H4));
xlabel('לךל','FontSize',14);
ylabel('לךל','FontSize',14);
colormap(1-gray);

figure; plot(ebsd, 'property',abs( H3));
xlabel('לךל','FontSize',14);
ylabel('לךל','FontSize',14);
colormap(1-gray);

grains0 = calcGrains(ebsd, 'threshold', 5*degree);
grains = calcGrains(grains0(grainSize(grains0) > 5), 'threshold', 2*degree);
figure();
KAM = calcKAM(grains0);
plot(get(grains0,'ebsd'),'property',KAM);
xlabel('לךל','FontSize',14);
ylabel('לךל','FontSize',14);
colormap(1-gray);