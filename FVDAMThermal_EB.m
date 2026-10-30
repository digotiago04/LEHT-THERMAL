function FVDAMThermal_EB(ny1, ny2, k_m, k_i, frac, field,  SVy1_cut, SVy2_cut)
% =========================================================================
% EFFECTIVE THERMAL CONDUCTIVITY AND TEMPERATURE FIELDS
%               FVDAM - ENERGY BASED
% =========================================================================
clc;

% 1. DIMENSIONS AND PARAMETERS
L = 1; H = 1;             % Dimensions of the repeating unit cell (RUC)    
Hy1_plot = 1;             % Macroscopic temperature gradient in y1 for plots
Hy2_plot = 0;             % Macroscopic temperature gradient in y2 for plots
l = L/ny1; h = H/ny2;     % Subvolume dimensions
R = sqrt(frac*L*H/pi);    % Radius of the inclusion
points = 3;               % Stride for plot markers in cut profiles


% 2. FINITE-VOLUME THEORY
% Periodic boundary conditions and degrees of freedom
[i, j] = meshgrid(1:ny1, 1:ny2);
s = (i + (j - 1)*ny1)';   % Subvolume numbering
faces = [s(:), s(:)+(ny1*ny2)+1, s(:)+ny1, s(:)+(ny1*ny2)];
faces(end-ny1+1:end, 3) = faces(1:ny1, 1);
faces(ny1:ny1:end, 2)    = faces(1:ny1:end-ny1+1, 4);
ndof = max(faces(:));   % Total number of faces
dofIS = unique([faces(1:ny1, 1); faces(end-ny1+1:end, 3)]);
dofIS = [dofIS(1), dofIS(end)];
dofDE = unique([faces(ny1:ny1:end, 2); faces(1:ny1:end-ny1+1, 4)]);
dofDE = [dofDE(1), dofDE(end)];
fixed = [dofIS dofDE];          % Fixed degrees
free = setdiff(1:ndof, fixed);  % Free degrees

% Sparse mapping indices
iK = reshape(kron(faces, ones(4,1))', 16*(ny1*ny2), 1);
jK = reshape(kron(faces, ones(1,4))', 16*(ny1*ny2), 1);
iF = repmat(faces', 2, 1);
jF = [ones(4, (ny1*ny2)); 2*ones(4, (ny1*ny2))];
k0 = eye(2);

% Auxiliary matrices
a = ones(4,1);
N1 = [0,-1]; N2 = [1,0]; N3 = [0,1]; N4 = [-1,0];
N = [N1, zeros(1,6); zeros(1,2), N2, zeros(1,4); ...
    zeros(1,4), N3, zeros(1,2); zeros(1,6), N4];
A = [0 -h/2 0 h^2/4; l/2 0 l^2/4 0; ...
    0 h/2 0 h^2/4; -l/2 0 l^2/4 0];
E = [0 0 0 0; 0 -1 0 3*h/2; -1 0 -3*l/2 0; 0 0 0 0; ...
    0 0 0 0; 0 -1 0 -3*h/2; -1 0 3*l/2 0; 0 0 0 0];
B = N * eye(8) * E;
ab = (B * (A\a)) \ (B/A);
Ab = A \ (eye(4) - a*ab);
K0 = B * Ab;
K0 = [K0(1,:)*l; K0(2,:)*h; K0(3,:)*l; K0(4,:)*h]; % Local thermal conductivity matrix

% Local load vector
H0 = [N1*l; N2*h; N3*l; N4*h] * k0;

% Material Design
[I, J] = ndgrid(1:ny1, 1:ny2);
Ic = I(:);
Jc = J(:);
% T0_x: Matrix containing the Y1 coordinate of the center of each face.
T0_x = [(0.5 + Ic - 1)*l, (1 + Ic - 1)*l, (0.5 + Ic - 1)*l, (Ic - 1)*l];
% T0_y: Matrix containing the Y2 coordinate of the center of each face.
T0_y = [ (Jc - 1)*h, (Jc - 0.5)*h,   (Jc)*h, (Jc - 0.5)*h ];
inclusion = ((I - 0.5)*l - L/2).^2 + ((J - 0.5)*h - H/2).^2 < R^2;
x = k_m * ones(ny1, ny2);
x(inclusion) = k_i;

% Material Interpolation
sK = K0(:) * x(:)';
sF = H0(:) * x(:)';

% Global thermal conductivity matrix
K = sparse(iK, jK, sK, ndof, ndof);
K = (K + K') / 2;

% Assembly of heat flux vectors corresponding to two unit temperature gradient tests
Q0 = sparse(iF(:), jF(:), sF, ndof, 2);

% Compute fluctuating temperatures for two unit temperature gradient tests
Tf = zeros(ndof, 2);
Tf(free, :) = K(free, free) \ Q0(free, :); 
Tfx = Tf(:, 1);
Tfy = Tf(:, 2);

Tx = T0_x + Tfx(faces);
Ty = T0_y + Tfy(faces);
T = cat(3, Tx, Ty);

% HOMOGENIZATION BASED ON ENERGY THEORY
C  = zeros(2);
for dir_i = 1:2
    Ti = T(:, :, dir_i);
    for dir_j = 1:2
        Tj = T(:, :, dir_j);
        sumE = sum((Ti * K0) .* Tj, 2) / (L * H);
        C(dir_i, dir_j) = -sum(x(:) .* sumE);
    end
end


% PRINT OUTPUT (K*)
fprintf('====================================================\n');
fprintf('EFFECTIVE THERMAL CONDUCTIVITY MATRIX (K*)\n');
fprintf('====================================================\n');
fprintf('FVDAM - ENERGY BASED \n');
disp(C);

% 3. POST-PROCESSING FOR PLOTS (FIELDS & CUTS)
doVerticalCut   = (SVy1_cut >= 1 && SVy1_cut <= ny1);
doHorizontalCut = (SVy2_cut >= 1 && SVy2_cut <= ny2);

% FVDAM Total Temperature Faces
Ttot_faces = Hy1_plot * Tx + Hy2_plot * Ty;
T00 = (ab * Ttot_faces')';

% ----- Extract Cuts -----
px = (SVy1_cut*l) - l/2;
py = (SVy2_cut*h) - h/2;

% Vertical Cut
if doVerticalCut
    svV = SVy1_cut + (0:ny2-1)'*ny1;
    yV  = ((1:ny2)' - 0.5)*h;
    TFVDAM_Y = T00(svV);
end

% Horizontal Cut
if doHorizontalCut
    svH = (1:ny1)' + (SVy2_cut-1)*ny1;
    xH  = ((1:ny1)' - 0.5)*l;
    TFVDAM_X = T00(svH);
end

% ----- Nodal Temperature -----
invA = A \ eye(4);                        
Taux = (Ttot_faces' - a*T00');
% Tcell_tot = Ttot_faces;   
% invA = A \ eye(4);                        
% T00  = (ab*Tcell_tot')';                      
% Taux = (Tcell_tot' - a*T00');                 
TijM = (invA*Taux)';                          
vert = [-l/2,-h/2; l/2,-h/2; l/2, h/2; -l/2, h/2];
phi = zeros(4,4); 
for p = 1:4
    x1 = vert(p,1); x2 = vert(p,2);
    phi(:,p) = [x1; x2; 0.5*(3*(x1^2) - (l^2)/4); 0.5*(3*(x2^2) - (h^2)/4)];
end
Tnodes = T00(:,ones(1,4)) + TijM*phi;    

nNx = ny1+1;  nNy = ny2+1;                     
nidBL = Ic + (Jc-1)*nNx;
nidBR = (Ic+1) + (Jc-1)*nNx;
nidTR = (Ic+1) + (Jc)*nNx;
nidTL = Ic + (Jc)*nNx;
nodeIDs = [nidBL, nidBR, nidTR, nidTL];       
sumT = accumarray(nodeIDs(:), Tnodes(:), [nNx*nNy,1], @sum, 0);
cntT = accumarray(nodeIDs(:), 1,              [nNx*nNy,1], @sum, 0);
Tnod = sumT ./ max(cntT,1);
Z_FVDAM = reshape(Tnod, [nNx, nNy])';           
[Xg,Yg] = meshgrid((0:ny1)*l, (0:ny2)*h);

% 4. PLOTTING

% Total Temperature Field - FVDAM
if field == 1
figure();
surf(Xg, Yg, Z_FVDAM);
view(2); shading flat; axis equal tight;
colormap('jet'); colorbar;
xlabel('$y_1$','Interpreter','latex','FontSize',16,'FontWeight','bold');
ylabel('$y_2$','Interpreter','latex','FontSize',16,'FontWeight','bold');
title('Total temperature field - FVDAM', 'Interpreter','latex','FontSize',16,'FontWeight','bold');
set(gca,'FontSize',12);
end

% Micro-temperature fields (Cuts)
if doHorizontalCut
figure();
hold on;
plot(xH, TFVDAM_X, 'bo', 'MarkerIndices', 1:points:length(xH));
grid on; xlim([0 L]); box on;
legend({'FVDAM'}, 'Interpreter', 'latex', 'Location', 'best');
xlabel('$y_1$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Temperature (K)', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
title(sprintf('Temperature field at $y_2 = %.3f$', py), 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
end 

if doVerticalCut
figure();
hold on;
plot(yV, TFVDAM_Y, 'bo', 'MarkerIndices', 1:points:length(yV));
grid on; xlim([0 H]); box on;
ylim('padded');
legend({'FVDAM'}, 'Interpreter', 'latex', 'Location', 'best');
xlabel('$y_2$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Temperature (K)', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
title(sprintf('Temperature field at $y_1 = %.3f$', px), 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
end
end
