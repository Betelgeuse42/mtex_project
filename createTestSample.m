function ebsd = createTestSample(fname, Xn, Yn, type, varargin)
%% Test data

[X,Y] = meshgrid(1:Xn,1:Yn);

gridType = get_option(varargin, 'gridType', 'Hex');

if strcmp(gridType,'Hex')
	X(2:2:end,:) = X(2:2:end,:) + 0.5;
    Y = Y*sqrt(3)/2;
end

switch (type)
    case 'random'
        ori = randomOri(X(:),Y(:), varargin{:});
    case 'lineGrad'
        ori = lineGradOri1(X(:),Y(:), varargin{:});
    case 'circGrad'
        ori = circGradOri(X(:),Y(:), varargin{:});
    case 'lineGrad2'
        ori = lineGrad2Ori1(X(:),Y(:), varargin{:});
    case 'sin'
        ori = sinOri1(X(:),Y(:), varargin{:});
    case 'tan'
        ori = tanOri1(X(:),Y(:), varargin{:});
    case 'island'
        ori = islandOri1(X(:),Y(:), varargin{:});
    otherwise
        error('Bad type');
end

CS = symmetry('m-3m');
SS = symmetry('m-3m');
ebsd = EBSD(ori, CS, SS);

ebsd = set(ebsd, 'x', X(:));
ebsd = set(ebsd, 'y', Y(:));
ebsd = set(ebsd,'phaseMap',1);
ebsd = set(ebsd,'phase',ones(length(ebsd),1));
ebsd = set(ebsd,'ci',ones(length(ebsd),1));

plotMisMap(ebsd, varargin{:});
xlabel('ìêì', 'FontSize', 14);
ylabel('ìêì','FontSize', 14);

export(ebsd, fname);
end

function plotMisMap(ebsd, varargin)

e1 = get_option(varargin, 'Euler1', [10 20 30]);

ori1 = orientation('Euler', e1(1)*degree, e1(2)*degree, e1(3)*degree, symmetry('m-3m'), symmetry('-1'));

ori = get(ebsd, 'orientation');
mis = ori.\ori1;

plot(ebsd, 'property', angle(mis)/degree);
xlabel('ìêì', 'FontSize', 12);
ylabel('ìêì','FontSize', 12);
end

function ori = lineGradOri1(X,Y, varargin)
% Set
%  Euler angles of point in degree
%  Axis and angle
%  Direction

e1 = get_option(varargin, 'Euler1', [10 20 30]);
ag = get_option(varargin, 'Amplitude', 30);
an = get_option(varargin, 'Noise', 0.0);

xmax = max(X);
ymax = max(Y);

% Variant 1
phi1 = (e1(1)+ag*X/xmax);
Phi  = (e1(2)-ag*(3*X+Y)/(3*xmax+ymax));
phi2 = (e1(3)+ag*Y/ymax);

phi1 = phi1 + an*randn(length(X),1);
Phi  = Phi  + an*randn(length(X),1);
phi2 = phi2 + an*randn(length(X),1);

ori = orientation('Euler', phi1*degree, Phi*degree, phi2*degree, symmetry('m-3m'), symmetry('-1'));
end

function ori = lineGrad2Ori1(X,Y, varargin)
% Set
%  Euler angles of point in degree
%  Axis and angle
%  Direction

e1 = get_option(varargin, 'Euler1', [10 20 30]);
e2 = get_option(varargin, 'Euler2', [ 5 79 16]);
ag = get_option(varargin, 'Amplitude', 30);
an = get_option(varargin, 'Noise', 0.0);

xmax = max(X);
ymax = max(Y);

ind = (X < xmax/2);

phi1 = zeros(1,length(X));
Phi  = zeros(1,length(X));
phi2 = zeros(1,length(X));

% Variant 1
phi1(ind) = (e1(1)+ag*X(ind)/xmax);
Phi(ind)  = (e1(2)-ag*(3*X(ind)+Y(ind))/(3*xmax+ymax));
phi2(ind) = (e1(3)+ag*Y(ind)/ymax);

ind = ~ind;
phi1(ind) = (e2(1)+ag*X(ind)/xmax);
Phi(ind)  = (e2(2)-ag*(3*X(ind)+Y(ind))/(3*xmax+ymax));
phi2(ind) = (e2(3)+ag*Y(ind)/ymax);

phi1 = phi1 + an*randn(length(X),1);
Phi  = Phi  + an*randn(length(X),1);
phi2 = phi2 + an*randn(length(X),1);

ori = orientation('Euler', phi1*degree, Phi*degree, phi2*degree, symmetry('m-3m'), symmetry('-1'));
end

function ori = sinOri1(X,Y, varargin)
% Sin gragient

% Get parameters
e1 = get_option(varargin, 'Euler1', [10 20 30]);
as = get_option(varargin, 'Amplitude', 30);
ns = get_option(varargin, 'Period', 5);
an = get_option(varargin, 'Noise', 0.0);
c  = get_option(varargin, 'Coeff', [0 0; 0 0; 1 0]);

% Get limits
xmax = max(X);
ymax = max(Y);

% Make gradient
phi1 = (e1(1)+c(1,1)*as*sin(ns*X/xmax*pi)+c(1,2)*as*sin(ns*Y/ymax*pi));
Phi  = (e1(2)+c(2,1)*as*sin(ns*X/xmax*pi)+c(2,2)*as*sin(ns*Y/ymax*pi));
phi2 = (e1(3)+c(3,1)*as*sin(ns*X/xmax*pi)+c(3,2)*as*sin(ns*Y/ymax*pi));

% Add noise
phi1 = phi1 + an*randn(length(X),1);
Phi  = Phi  + an*randn(length(X),1);
phi2 = phi2 + an*randn(length(X),1);

% Set orientation
ori = orientation('Euler', phi1*degree, Phi*degree, phi2*degree, symmetry('m-3m'), symmetry('-1'));
end

function ori = tanOri1(X,Y, varargin)
% Tan gragient

% Get parameters
e1 = get_option(varargin, 'Euler1', [10 20 30]);
as = get_option(varargin, 'Amplitude', 30);
ns = get_option(varargin, 'Period', 5);
an = get_option(varargin, 'Noise', 0.0);
c  = get_option(varargin, 'Coeff', [0 0; 0 0; 1 0]);

% Get limits
xmax = max(X);
ymax = max(Y);

% Make gradient
phi1 = (e1(1)+c(1,1)*as*tan(ns*X/xmax*pi)+c(1,2)*as*tan(ns*Y/ymax*pi));
Phi  = (e1(2)+c(2,1)*as*tan(ns*X/xmax*pi)+c(2,2)*as*tan(ns*Y/ymax*pi));
phi2 = (e1(3)+c(3,1)*as*tan(ns*X/xmax*pi)+c(3,2)*as*tan(ns*Y/ymax*pi));

% Add noise
phi1 = phi1 + an*randn(length(X),1);
Phi  = Phi  + an*randn(length(X),1);
phi2 = phi2 + an*randn(length(X),1);

% Set orientation
ori = orientation('Euler', phi1*degree, Phi*degree, phi2*degree, symmetry('m-3m'), symmetry('-1'));
end

function ori = islandOri1(X,Y, varargin)
% Tan gragient

% Get parameters
e1 = get_option(varargin, 'Euler1', [10 20 30]);
as = get_option(varargin, 'Amplitude', 30);
ns = get_option(varargin, 'Period', 5);
an = get_option(varargin, 'Noise', 0.0);
c  = get_option(varargin, 'Coeff', [0 0; 0 0; 1 0]);
h  = get_option(varargin, 'High', 0);

% Get limits
xmax = max(X);
ymax = max(Y);

% Make gradient
phi1 = (e1(1)+c(1,1)*as*(sin(ns*X/xmax*pi)>h)+c(1,2)*as*(sin(ns*Y/ymax*pi)>h));
Phi  = (e1(2)+c(2,1)*as*(sin(ns*X/xmax*pi)>h)+c(2,2)*as*(sin(ns*Y/ymax*pi)>h));
phi2 = (e1(3)+c(3,1)*as*(sin(ns*X/xmax*pi)>h)+c(3,2)*as*(sin(ns*Y/ymax*pi)>h));

% Add noise
phi1 = phi1 + an*randn(length(X),1);
Phi  = Phi  + an*randn(length(X),1);
phi2 = phi2 + an*randn(length(X),1);

% Set orientation
ori = orientation('Euler', phi1*degree, Phi*degree, phi2*degree, symmetry('m-3m'), symmetry('-1'));
end