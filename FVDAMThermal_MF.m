function FVDAMThermal_MF(ny1, ny2, k_m, k_i, frac, field, SVy1_cut, SVy2_cut)
% =========================================================================
% CALCULATION OF THE EFFECTIVE THERMAL CONDUCTIVITY AND TEMPERATURE FIELDS
%                      FVDAM - MEAN FIELD
% =========================================================================
clc;

% 1. DIMENSIONS AND PARAMETERS
L = 1; H = 1;             % Dimensions of the repeating unit cell (RUC)
Hy1_plot = 1;             % Macroscopic temperature gradient in y1 for plots
Hy2_plot = 0;             % Macroscopic temperature gradient in y2 for plots
l = L/ny1; h = H/ny2;     % Subvolume dimensions     
R = sqrt(L*H*frac/pi);    % Radius of the inclusion
points = 3;               % Stride for plot markers in cut profiles


% 2. FINITE-VOLUME THEORY
% Periodic boundary conditions and degrees of freedom
[i,j] = meshgrid(1:ny1,1:ny2);
s = (i+(j-1)*ny1)';  % Subvolume numbering
faces = [s(:),s(:)+ny1*ny2+1,s(:)+ny1,s(:)+ny1*ny2];
faces(end-ny1+1:end,3) = faces(1:ny1,1);
faces(ny1:ny1:end,2) = faces(1:ny1:end-ny1+1,4);
ndof = max(faces(:)); % Total number of faces 
dofIS = unique([faces(1:ny1, 1); faces(end-ny1+1:end, 3)]);
dofIS = [dofIS(1),dofIS(end)];
dofDE = unique([faces(ny1:ny1:end, 2); faces(1:ny1:end-ny1+1, 4)]);
dofDE = [dofDE(1),dofDE(end)];
fixed = [dofIS dofDE];          % Fixed degrees 
free = setdiff(1:ndof,fixed);   % Free degrees

% Sparse mapping indices
iK = reshape(kron(faces,ones(4,1))',16*ny1*ny2,1);
jK = reshape(kron(faces,ones(1,4))',16*ny1*ny2,1);
iF = repmat((faces)',2,1);
jF = [ones(4,ny1*ny2); 2*ones(4,ny1*ny2)];

% Auxiliary matrices 
k0 = eye(2);
a_vec = ones(4,1);
N1 = [0,-1]; N2 = [1,0]; N3 = [0,1]; N4 = [-1,0];
N_mat = [N1,zeros(1,6); zeros(1,2),N2,zeros(1,4); zeros(1,4),N3,zeros(1,2); zeros(1,6),N4];
A = [0 -h/2 0 h^2/4; l/2 0 l^2/4 0; 0 h/2 0 h^2/4; -l/2 0 l^2/4 0];
E = [0 0 0 0; 0 -1 0 3*h/2; -1 0 -3*l/2 0; 0 0 0 0; 0 0 0 0; 0 -1 0 -3*h/2; -1 0 3*l/2 0; 0 0 0 0];
B = N_mat*eye(8)*E;
ab = (B*(A\a_vec))\(B/A);
Ab = A\(eye(4)-a_vec*ab);
K0 = B*Ab;
K0 = [K0(1,:)*l;K0(2,:)*h;K0(3,:)*l;K0(4,:)*h]; % Local thermal conductivity matrix
H0 = [N1*l;N2*h;N3*l;N4*h]*k0;
G0 = eye(2);

% Material Design
x = InitialMaterialDesign(ny1, ny2, l, h, L, H, R, k_m, k_i);

% Material Interpolation
sK = K0(:)*x(:)';  
sF = H0(:)*x(:)';

% Global thermal conductivity matrix
K = sparse(iK, jK, sK, ndof, ndof);
K = (K+K')/2;

% Assembly of heat flux vectors corresponding to two unit temperature gradient tests
Q0 = sparse(iF(:), jF(:), sF, ndof, 2);

% Compute fluctuating temperatures for two unit temperature gradient tests
T_til = zeros(ndof,2);
T_til(free,:) = K(free,free) \ Q0(free,:);

% Calculates effective thermal conductivity from local temperature gradients.
Tx = (T_til(faces(:,2), :) - T_til(faces(:,4), :)) / l;
Ty = (T_til(faces(:,3), :) - T_til(faces(:,1), :)) / h;
K_ast_FVDAM = G0 * mean(x(:)) + [x(:)' * Tx; x(:)' * Ty] / (ny1*ny2);



% PRINT OUTPUT (K*)
fprintf('====================================================\n');
fprintf('EFFECTIVE THERMAL CONDUCTIVITY MATRIX (K*)\n');
fprintf('====================================================\n');
fprintf('FVDAM - MEAN-FIELD \n');
disp(K_ast_FVDAM);

% 3. POST-PROCESSING FOR PLOTS (FIELDS & CUTS) e % 4. PLOTTING
doVerticalCut   = (SVy1_cut >= 1 && SVy1_cut <= ny1);
doHorizontalCut = (SVy2_cut >= 1 && SVy2_cut <= ny2);

% FVT Total Temperature Faces
kcase = 1; 
if Hy1_plot==0 && Hy2_plot==1, kcase = 2; end
[I_mesh,J_mesh] = ndgrid(1:ny1,1:ny2);
I_vec = I_mesh(:);
J_vec = J_mesh(:);
Xfaces = [(I_vec-0.5)*l, I_vec*l, (I_vec-0.5)*l, (I_vec-1)*l];
Yfaces = [(J_vec-1)*h, (J_vec-0.5)*h, J_vec*h, (J_vec-0.5)*h];
Ttil_case = reshape(T_til(faces(:),kcase), ny1*ny2, 4);
Ttot_faces = Xfaces*Hy1_plot + Yfaces*Hy2_plot + Ttil_case;


% ----- Nodal Temperature -----  
invA = A \ eye(4);                        
T00  = (ab*Ttot_faces')';                      
Taux = (Ttot_faces' - a_vec*T00');                 
TijM = (invA*Taux)';                          
vert = [-l/2,-h/2; l/2,-h/2; l/2, h/2; -l/2, h/2];
phi = zeros(4,4); 
for p = 1:4
    x1l = vert(p,1); x2l = vert(p,2);
    phi(:,p) = [ x1l; x2l; 0.5*(3*(x1l^2) - (l^2)/4); 0.5*(3*(x2l^2) - (h^2)/4) ];
end
Tnodes_cell = T00(:,ones(1,4)) + TijM*phi;    
nNx = ny1+1;  nNy = ny2+1;
Ic = I_mesh(:);  Jc = J_mesh(:);                        
nidBL = Ic + (Jc-1)*nNx;
nidBR = (Ic+1) + (Jc-1)*nNx;
nidTR = (Ic+1) + (Jc)*nNx;
nidTL = Ic + (Jc)*nNx;
nodeIDs = [nidBL, nidBR, nidTR, nidTL];       
sumT = accumarray(nodeIDs(:), Tnodes_cell(:), [nNx*nNy,1], @sum, 0);
cntT = accumarray(nodeIDs(:), 1, [nNx*nNy,1], @sum, 0);
Tnod = sumT ./ max(cntT,1);
T_FVDAM = reshape(Tnod, [nNx, nNy])';           
[Xg,Yg] = meshgrid((0:ny1)*l, (0:ny2)*h);


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

% 4. PLOTTING

% Total Temperature Field - FVDAM
if field == 1
figure();
surf(Xg, Yg, T_FVDAM);
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


% AUXILIARY FUNCTIONS 

% Material design - FVDAM
function x = InitialMaterialDesign(nx,ny,l,h,L,H,R,k_m,k_i)
    [I, J] = ndgrid(1:nx, 1:ny);
    x = k_m * ones(nx, ny);
    inclusion = ((I - 0.5)*l - L/2).^2 + ((J - 0.5)*h - H/2).^2 < R^2;
    x(inclusion) = k_i;
end
