function LEHT_Thermal(k_m, k_i, frac, field, y1_cut, y2_cut)

% 1. DIMENSIONS AND NUMERICAL PARAMETERS
L = 1;           % CELL LENGTH
H = 1;           % CELL HEIGHT
N = 60;          % HARMONICS 
Mbc = 300;       % COLLOCATION POINTS PER FACE
a = sqrt(frac*L*H/pi);
if ~(frac > 0 && a < min(L,H)/2)
    error('LEHT:frac', 'frac must satisfy 0 < frac < pi/4 (= 0.785) for a centered inclusion in a square RUC.');
end

% 2. EFFECTIVE THERMAL CONDUCTIVITY MATRIX K* 
K_ast = zeros(2); coefs = cell(1,2);
for col = 1:2
    Hv = [0 0]; Hv(col) = 1;
    coefs{col} = LEHT_isotropic_square(L, H, L/2, H/2, a, k_m, k_i, N, Mbc, Hv(1), Hv(2));
    [qx, qy] = mean_flux_fibre(coefs{col}, L, H, Hv(1), Hv(2));
    K_ast(:,col) = -[qx; qy];
end

fprintf('EFFECTIVE THERMAL CONDUCTIVITY MATRIX (K*) \n');
disp(K_ast);

% 3. POST-PROCESSING AND VISUALIZATION 
if field == 1 || (y1_cut > 0 && y1_cut <= L) || (y2_cut > 0 && y2_cut <= H)
    plot_thermal_results(coefs{1}, L, H, 1, 0, field, y1_cut, y2_cut);
end
end

function plot_thermal_results(coef, L, H, Hx, Hy, field, y1_cut, y2_cut)
T_func = @(x, y) evaluate_temperature(x, y, coef, Hx, Hy);

if field == 1
    [Xg, Yg] = meshgrid(linspace(0, L, 150), linspace(0, H, 150));
    Z_LEHT = T_func(Xg, Yg);
    figure();
    surf(Xg, Yg, Z_LEHT, 'EdgeColor', 'none');
    view(2); axis equal tight;
    colormap('jet'); colorbar;
    title('Total temperature field - LEHT', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
    xlabel('$y_1$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('$y_2$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
    set(gca,'FontSize',12);
end

if y2_cut > 0 && y2_cut <= H
    xH = linspace(0, L, 500);
    TLEHT_X = T_func(xH, y2_cut * ones(size(xH)));
    figure();
    plot(xH, TLEHT_X, 'k-', 'LineWidth', 1.5);
    grid on; xlim([0 L]); box on;
    xlabel('$y_1$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('Temperature (K)', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
    title(sprintf('Temperature field at $y_2 = %.2f$', y2_cut), 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
end

if y1_cut > 0 && y1_cut <= L
    yV = linspace(0, H, 500);
    TLEHT_Y = T_func(y1_cut * ones(size(yV)), yV);
    figure();
    plot(yV, TLEHT_Y, 'k-', 'LineWidth', 1.5);
    grid on; xlim([0 H]); box on;
    xlabel('$y_2$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('Temperature (K)', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
    title(sprintf('Temperature field at $y_1 = %.2f$', y1_cut), 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
end
end

% LEHT FUNCTIONS
function T_total = evaluate_temperature(x, y, coef, Hx, Hy)
sz = size(x); x = x(:); y = y(:);
z  = complex(x - coef.xc, y - coef.yc) / coef.a;
n  = 1:coef.N;
T_tilde = zeros(size(x));
mf = abs(z) <= 1; mm = ~mf;
if any(mf)
    Zp = z(mf).^n;
    T_tilde(mf) = coef.a * (real(Zp)*coef.Af + imag(Zp)*coef.Bf);
end
if any(mm)
    Zp = z(mm).^n; Zn = z(mm).^(-n);
    T_tilde(mm) = coef.a * (real(Zp)*coef.Am + imag(Zp)*coef.Bm + real(Zn)*coef.Cm - imag(Zn)*coef.Dm);
end
T_total = reshape(Hx*x + Hy*y + T_tilde, sz);
end

function coef = LEHT_isotropic_square(L, H, xc, yc, a, km, ki, N, Mbc, Hx, Hy)
if a <= 0 || a >= min([xc, L-xc, yc, H-yc])
    error('LEHT:geometry', 'The inclusion must lie strictly inside the RUC.');
end
kap = ki/km;
e1  = double((1:N).' == 1);
cp  = (1+kap)/2; cq = (1-kap)/2;
I = eye(N); O = zeros(N);
P = [cp*I, O; O, cp*I; cq*I, O; O, cq*I];              
f = [-cq*Hx*e1; -cq*Hy*e1; cq*Hx*e1; cq*Hy*e1];

t  = ((1:Mbc).' - 0.5) / Mbc;
on = ones(Mbc,1);
[TR, GxR]     = basis_m(L*on, t*H,  xc, yc, a, N);
[TL, GxL]     = basis_m(0*on, t*H,  xc, yc, a, N);
[TT, ~, GyT]  = basis_m(t*L,  H*on, xc, yc, a, N);
[TB, ~, GyB]  = basis_m(t*L,  0*on, xc, yc, a, N);
lref = sqrt(L*H);
D = [ (TR - TL)/lref;  (TT - TB)/lref;                 
       GxR - GxL;       GyT - GyB ];                     
A = D*P; b = -D*f;
s = sqrt(sum(A.^2, 1)); s(s == 0) = 1;                   
g = ((A./s) \ b) ./ s.';                                 
c = P*g + f;
coef.Af = g(1:N);        coef.Bf = g(N+1:2*N);
coef.Am = c(1:N);        coef.Bm = c(N+1:2*N);
coef.Cm = c(2*N+1:3*N);  coef.Dm = c(3*N+1:4*N);
coef.N = N; coef.a = a; coef.km = km; coef.ki = ki; coef.xc = xc; coef.yc = yc;
end

function [PT, PGx, PGy] = basis_m(x, y, xc, yc, a, N)
z  = complex(x(:) - xc, y(:) - yc) / a;
n  = 1:N;
Zp = z.^n;          Zn = z.^(-n);
Dp = n.*z.^(n-1);   Dn = -n.*z.^(-n-1);                 
PT  = a*[real(Zp), imag(Zp), real(Zn), -imag(Zn)];
PGx =   [real(Dp), imag(Dp), real(Dn), -imag(Dn)];
PGy =   [-imag(Dp), real(Dp), -imag(Dn), -real(Dn)];
end

function [qx, qy] = mean_flux_fibre(coef, L, H, Hx, Hy)
vf = pi*coef.a^2 / (L*H);
qx = -coef.km*Hx - (coef.ki - coef.km)*vf*(Hx + coef.Af(1));
qy = -coef.km*Hy - (coef.ki - coef.km)*vf*(Hy + coef.Bf(1));
end
