function LEHT_Thermal(k_m, k_i, frac, field, y1_cut, y2_cut)

clc;

% 1. DIMENSIONS AND PARAMETERS
L = 1; H = 1;            % Dimensions of the Representative Unit Cell (RUC)
Hx = 1;                  % Macroscopic temperature gradient in y1 for plots
Hy = 0;                  % Macroscopic temperature gradient in y2 for plots
R = sqrt((L*H*frac)/pi); % Radius of the inclusion  

% 2. LEHT
K_ast_LEHT = zeros(2);
for col = 1:2
    if col==1, Hx=1; Hy=0; else, Hx=0; Hy=1; end
    coef = LEHT_isotropic_square(L, H, L/2, H/2, R, k_m, k_i, 30, 250, 400, Hx, Hy); 
    [qxbar, qybar] = mean_flux_LEHT_boundary(coef, L, H, k_m, Hx, Hy, 1000);
    K_ast_LEHT(:,col) = -[qxbar; qybar];
end

coef_plot = LEHT_isotropic_square(L, H, L/2, H/2, R, k_m, k_i, 30, 250, 400, Hx, Hy);
points_plot = 100; 
xH = linspace(0, L, points_plot)';
yV = linspace(0, H, points_plot)';
TLEHT_Y = evaluate_temperature(y1_cut*ones(size(yV)), yV, coef_plot, Hx, Hy);
TLEHT_X = evaluate_temperature(xH, y2_cut*ones(size(xH)), coef_plot, Hx, Hy);

% 3. EFFECTIVE THERMAL CONDUCTIVITY MATRICES
fprintf('====================================================\n');
fprintf('EFFECTIVE THERMAL CONDUCTIVITY MATRICES (K*)\n');
fprintf('====================================================\n');
disp(K_ast_LEHT);
             
% 4. PLOTTING

% Total Temperature Field
if field == 1
    [Xg, Yg] = meshgrid(linspace(0, L, points_plot), linspace(0, H, points_plot));  
    T_col = evaluate_temperature(Xg(:), Yg(:), coef_plot, Hx, Hy);
    T_LEHT = reshape(real(T_col), size(Xg))';
    figure(1);
    surf(Xg, Yg, T_LEHT); 
    view(2);            
    shading flat;     
    colormap('jet');    
    colorbar;              
    xlabel('$y_1$','Interpreter','latex','FontSize',16,'FontWeight','bold');
    ylabel('$y_2$','Interpreter','latex','FontSize',16,'FontWeight','bold');
    title('Total temperature field - LEHT', 'Interpreter','latex','FontSize',16,'FontWeight','bold');
    set(gca,'FontSize',12);
end

% Micro-temperature fields
if y2_cut > 0 && y2_cut <= 1
figure(2);
plot(xH, TLEHT_X, 'k-', 'LineWidth', 1); hold on;
grid on; xlim([0 L]); box on;
ylim([0 1]);
legend({'LEHT'}, 'Interpreter', 'latex', 'Location', 'best');
xlabel('$y_1$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Temperature ($^\circ$C)', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
title(sprintf('Temperature field at $y_2 = %.2f$', y2_cut), 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
end

if y1_cut > 0 && y1_cut <= 1
figure(3);
plot(yV, TLEHT_Y, 'k-', 'LineWidth', 1); hold on;
grid on; xlim([0 H]); box on;
ylim([0 1]);
legend({'LEHT'}, 'Interpreter', 'latex', 'Location', 'best');
xlabel('$y_2$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Temperature ($^\circ$C)', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
title(sprintf('Temperature field at $y_1 = %.2f$', y1_cut), 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
end
end


% LEHT functions
function coef = LEHT_isotropic_square(L,H,xc,yc,a,km,ki,N,Mbc,Mint,Hx,Hy)
    yy = (1:Mbc)'/(Mbc+1)*H;
    xx = (1:Mbc)'/(Mbc+1)*L;
    xL = zeros(size(yy));   yL = yy;
    xR = L + zeros(size(yy)); yR = yy;
    xB = xx;                yB = zeros(size(xx));
    xT = xx;                yT = H + zeros(size(xx));
    
    thI = (0:Mint-1)'*(2*pi/Mint);
    
    A_bcT = [ periodic_Ttilde_rows(xR,yR,xL,yL,xc,yc,a,N); ...
              periodic_Ttilde_rows(xT,yT,xB,yB,xc,yc,a,N) ];
    b_bcT = zeros(size(A_bcT,1),1);
    
    A_bcQ = [ periodic_flux_rows('x',xR,yR,xL,yL,xc,yc,a,km,N); ...
              periodic_flux_rows('y',xT,yT,xB,yB,xc,yc,a,km,N) ];
    b_bcQ = zeros(size(A_bcQ,1),1);
    
    [A_intT, b_intT] = interface_T_rows(thI,a,N);
    [A_intQ, b_intQ] = interface_qn_rows(thI,km,ki,N,Hx,Hy);
    
    A_mat = [A_bcT; A_bcQ; A_intT; A_intQ];
    b_mat = [b_bcT; b_bcQ; b_intT; b_intQ];
    
    normA = sqrt(sum(A_mat.^2, 1)); 
    normA(normA == 0) = 1; 
    A_scaled = A_mat ./ normA; 
    y = lsqminnorm(A_scaled, b_mat);
    u = y ./ normA';  
    coef.Af = u(1:N);       coef.Bf = u(N+1:2*N);
    coef.Am = u(2*N+1:3*N); coef.Bm = u(3*N+1:4*N);
    coef.Cm = u(4*N+1:5*N); coef.Dm = u(5*N+1:6*N);
    coef.N=N; coef.a=a; coef.km=km; coef.ki=ki; coef.xc=xc; coef.yc=yc;
end

function T_total = evaluate_temperature(x, y, coef, Hx, Hy)
    rx = x - coef.xc;
    ry = y - coef.yc;
    r = sqrt(rx.^2 + ry.^2);
    th = atan2(ry, rx);
    T_tilde = zeros(size(x)); 
    mask_f = (r <= coef.a);
    mask_m = ~mask_f;
    
    if any(mask_f, 'all')
        xi_f = r(mask_f) / coef.a; 
        th_f = th(mask_f);
        T_f = zeros(size(xi_f));
        for n = 1:coef.N
            T_f = T_f + coef.a * (coef.Af(n)*(xi_f.^n).*cos(n*th_f) + coef.Bf(n)*(xi_f.^n).*sin(n*th_f));
        end
        T_tilde(mask_f) = T_f;
    end
    
    if any(mask_m, 'all')
        xi_m = max(r(mask_m) / coef.a, 1e-14); 
        th_m = th(mask_m);
        T_m = zeros(size(xi_m));
        for n = 1:coef.N
            T_m = T_m + coef.a * (...
                coef.Am(n)*(xi_m.^n).*cos(n*th_m) + coef.Bm(n)*(xi_m.^n).*sin(n*th_m) + ...
                coef.Cm(n)*(xi_m.^(-n)).*cos(n*th_m) + coef.Dm(n)*(xi_m.^(-n)).*sin(n*th_m));
        end
        T_tilde(mask_m) = T_m;
    end
    T_total = Hx.*x + Hy.*y + T_tilde;
end

function [qxbar,qybar] = mean_flux_LEHT_boundary(coef,L,H,km,Hx,Hy,M)
    yy = (1:M)'/(M+1)*H;              
    xx = (1:M)'/(M+1)*L;              
    xR = L * ones(size(yy));  
    [dTdxR, ~] = eval_grad_Ttil_matrix(xR, yy, coef);
    qxbar = mean(-km * (Hx + dTdxR)); 
    yT = H * ones(size(xx));
    [~, dTdyT] = eval_grad_Ttil_matrix(xx, yT, coef);
    qybar = mean(-km * (Hy + dTdyT)); 
end

function [dTdx, dTdy] = eval_grad_Ttil_matrix(x,y,coef)  
    rx = x - coef.xc;
    ry = y - coef.yc;
    r  = sqrt(rx.^2 + ry.^2);
    th = atan2(ry,rx);
    n = 1:coef.N;
    C = cos(th*n);
    S = sin(th*n);
    Cn = C.*n;
    Sn = S.*n;
    xi = r/coef.a;
    xiP  = bsxfun(@power, xi,  n);
    xiM  = bsxfun(@power, xi, -n);
    xiP1 = bsxfun(@power, xi,  n-1);
    xiM1 = bsxfun(@power, xi, -n-1);
    dTdr  = (xiP1.*Cn)*coef.Am + (xiP1.*Sn)*coef.Bm - (xiM1.*Cn)*coef.Cm - (xiM1.*Sn)*coef.Dm;
    dTdth = coef.a * ( -(xiP.*Sn)*coef.Am + (xiP.*Cn)*coef.Bm - (xiM.*Sn)*coef.Cm + (xiM.*Cn)*coef.Dm );
    ct = cos(th);
    st = sin(th);
    invr = 1./r; 
    dTdx = dTdr.*ct - (dTdth.*st).*invr;
    dTdy = dTdr.*st + (dTdth.*ct).*invr;
end

function Arows = periodic_Ttilde_rows(x2,y2,x1,y1,xc,yc,a,N)
    [PhiAm2,PhiBm2,PhiCm2,PhiDm2] = matrix_Ttilde_basis(x2,y2,xc,yc,a,N);
    [PhiAm1,PhiBm1,PhiCm1,PhiDm1] = matrix_Ttilde_basis(x1,y1,xc,yc,a,N);
    dAm = PhiAm2 - PhiAm1;
    dBm = PhiBm2 - PhiBm1;
    dCm = PhiCm2 - PhiCm1;
    dDm = PhiDm2 - PhiDm1;
    Z = zeros(size(dAm));
    Arows = [Z Z dAm dBm dCm dDm];
end

function [PhiAm,PhiBm,PhiCm,PhiDm] = matrix_Ttilde_basis(x,y,xc,yc,a,N)
    rx = x - xc; ry = y - yc;
    r = sqrt(rx.^2+ry.^2); th = atan2(ry,rx);
    xi = max(r/a, 1e-14);
    n = 1:N;
    C = cos(th*n); S = sin(th*n);
    xiP = bsxfun(@power,xi,n);
    xiM = bsxfun(@power,xi,-n);
    PhiAm = a*(xiP.*C);
    PhiBm = a*(xiP.*S);
    PhiCm = a*(xiM.*C);
    PhiDm = a*(xiM.*S);
end

function Arows = periodic_flux_rows(dir,x2,y2,x1,y1,xc,yc,a,km,N)
    [QxAm2,QxBm2,QxCm2,QxDm2, QyAm2,QyBm2,QyCm2,QyDm2] = matrix_flux_basis(x2,y2,xc,yc,a,km,N);
    [QxAm1,QxBm1,QxCm1,QxDm1, QyAm1,QyBm1,QyCm1,QyDm1] = matrix_flux_basis(x1,y1,xc,yc,a,km,N);
    
    if dir == 'x'
        dAm = QxAm2 - QxAm1; dBm = QxBm2 - QxBm1;
        dCm = QxCm2 - QxCm1; dDm = QxDm2 - QxDm1;
    else
        dAm = QyAm2 - QyAm1; dBm = QyBm2 - QyBm1;
        dCm = QyCm2 - QyCm1; dDm = QyDm2 - QyDm1;
    end
    
    Z = zeros(size(dAm));
    Arows = [Z Z dAm dBm dCm dDm];
end

function [QxAm,QxBm,QxCm,QxDm, QyAm,QyBm,QyCm,QyDm] = matrix_flux_basis(x,y,xc,yc,a,km,N)
    rx = x - xc; ry = y - yc;
    r = sqrt(rx.^2+ry.^2); th = atan2(ry,rx);
    xi = max(r/a, 1e-14);
    n = 1:N;
    C = cos(th*n); S = sin(th*n);
    Cn = C.*n; Sn = S.*n;
    xiP  = bsxfun(@power,xi,n);
    xiM  = bsxfun(@power,xi,-n);
    xiP1 = bsxfun(@power,xi,n-1);
    xiM1 = bsxfun(@power,xi,-n-1);
    ct = cos(th); st = sin(th);
    invr = 1./max(r, 1e-14);
   
    dTdr_Am = (xiP1).*Cn;     dTdth_Am = -a*(xiP).*Sn;
    dTdr_Cm = -(xiM1).*Cn;    dTdth_Cm = -a*(xiM).*Sn;
    dTdr_Bm = (xiP1).*Sn;     dTdth_Bm =  a*(xiP).*Cn;
    dTdr_Dm = -(xiM1).*Sn;    dTdth_Dm =  a*(xiM).*Cn;
    
    dTdx_Am = dTdr_Am.*ct - (dTdth_Am.*st).*invr;
    dTdx_Bm = dTdr_Bm.*ct - (dTdth_Bm.*st).*invr;
    dTdx_Cm = dTdr_Cm.*ct - (dTdth_Cm.*st).*invr;
    dTdx_Dm = dTdr_Dm.*ct - (dTdth_Dm.*st).*invr;
    
    dTdy_Am = dTdr_Am.*st + (dTdth_Am.*ct).*invr;
    dTdy_Bm = dTdr_Bm.*st + (dTdth_Bm.*ct).*invr;
    dTdy_Cm = dTdr_Cm.*st + (dTdth_Cm.*ct).*invr;
    dTdy_Dm = dTdr_Dm.*st + (dTdth_Dm.*ct).*invr;
    
    QxAm = -km*dTdx_Am; QxBm = -km*dTdx_Bm; QxCm = -km*dTdx_Cm; QxDm = -km*dTdx_Dm;
    QyAm = -km*dTdy_Am; QyBm = -km*dTdy_Bm; QyCm = -km*dTdy_Cm; QyDm = -km*dTdy_Dm;
end

function [Arows, brows] = interface_T_rows(thI,a,N)
    n = 1:N;
    C = cos(thI*n);
    S = sin(thI*n);
    Arows = [ -a*C, -a*S,  a*C,  a*S,  a*C,  a*S ];
    brows = zeros(size(thI));
end

function [Arows, brows] = interface_qn_rows(thI,km,ki,N,Hx,Hy)
    n = 1:N;
    C = cos(thI*n);
    S = sin(thI*n);
    Cn = C.*n;
    Sn = S.*n;
    Arows = [ -ki*Cn, -ki*Sn,  km*Cn,  km*Sn, -km*Cn, -km*Sn ];
    dT0dr = Hx*cos(thI) + Hy*sin(thI);
    brows = -(km-ki)*dT0dr;
end