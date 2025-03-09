% Solving Parabolic Wave equation with finite differences and
% Crank Nicolson, with PML

% Convergence rate
% Convergence Analysis
% clc;
% dz_v = [0.1 0.05 0.025];
% dx_v = [10 5 2.5];
% error = zeros(length(dz_v)-1,1);
% U_old = 0;
% 
% % Nyqvist criterion
% lambda = 0.1;
% alpha = atan(50/250);
% deltaz = lambda/(4*sin(alpha));
% deltax = lambda/(8*sin(0.5*alpha)^2);
% 
% for m = 1:length(dz_v)
%     dz = dz_v(m);
%     dx = dx_v(m);
%     tic
%     [U, x, z] = wavefunc(dx, dz);
%     toc
%     drawnow
%     if m > 1
%         U_coarse = U(1:2:end, 1:2:end);
%         matrix_diff = U_old - U_coarse;
%         error(m-1) = norm(matrix_diff);
%     end
% 
%     U_old = U;
%     xold = x;
%     zold = z;
% end
% 
% figure;
% loglog(dx_v(1:end-1), error);
% hold on;
% loglog(dx_v(1:end-1), dx_v(1:end-1).^2);
% legend('Method', 'dx^2');
% xlabel('N');
% ylabel('Error');
% 
% order = log(error(1) / error(2)) / log(2);
% disp(['Computed order of accuracy: ', num2str(order)]);


%%

dx = 0.5;
dz = 0.5;

z0 = 0;
zend = 100;

x0 = 0;
xend = 10000;

z = z0:dz:zend;
x = x0:dx:xend;

Nz = length(z);
Nx = length(x);

%tStart = tic;
[U_coarse] = wavefunc(dx,dz,xend,zend);
% figure
% surf(abs(U_coarse))
% shading interp
%toc
%fprintf('Time elapsed: %.6f seconds\n', toc(tStart));
%tStart = tic;
[U_fine] = wavefunc(dx/2,dz/2,xend,zend);
% figure
% surf(abs(U_fine))
% shading interp
%toc
%fprintf('Time elapsed: %.6f seconds\n', toc(tStart));
%tStart = tic;
[U_finest] = wavefunc(dx/4,dz/4,xend,zend);
% figure
% surf(abs(U_finest))
% shading interp
%toc
%fprintf('Time elapsed: %.6f seconds\n', toc(tStart));
%tic

U_coarse_interp = U_fine(1:2:end, 1:2:end);  % Take every second point
error_1 =  norm(U_coarse - U_coarse_interp)/sqrt(Nx*Nz);
max_e_1 = (max(max(abs(U_coarse - U_coarse_interp))));

U_fine_interp = U_finest(1:2:end, 1:2:end);
error_2 =  norm(U_fine - U_fine_interp)/sqrt(2*Nx*2*Nz);
max_e_2 = (max(max(abs(U_fine - U_fine_interp))));

order = log(error_1 / error_2) / log(2);
disp(['Computed order of accuracy: ', num2str(order)]);
%toc

save('U_coarse.mat', 'U_coarse');
save('U_fine.mat', 'U_fine');
save('U_finest.mat', 'U_finest');


%%
% Plot errors

% figure
% surf(abs(U_coarse - U_coarse_interp))
% shading interp
% 
% figure
% surf(abs(U_fine - U_fine_interp))
% shading interp



%% 


function [U] = wavefunc(dx, dz, xend, zend)
    
    z = 0:dz:zend;
    x = 0:dx:xend;

    Nz = length(z);
    Nx = length(z);
    
    lambda = 0.1;    % Radar wavelength
    wk = 2 * pi / lambda; % Wave number
    ra = 6360000;    % Radius of the earth
    
    x = linspace(0, xend, Nx);
    z = linspace(0, zend, Nz);
    
    % Initial condition
    beta = 3 * pi / 180;
    theta0 = 0;
    ha = 31;
    k = wk;
    u_fs = @(z) (k * beta) / (2 * sqrt(2 * pi * log(2))) .* exp(-1i * k * theta0 * z) .* exp(-(beta^2 / (8 * log(2))) * k^2 * (z - ha).^2);
    U = zeros(Nz, Nx);
    U(:, 1) = u_fs(z) - u_fs(-z);
    
    % Absorption function sigma(z)                                      
    d = 40;
    sigma_max = 20;
    z1 = zend - d;

% Absorption function sigma(z)
sigma = zeros(Nz, 1);
for j = 1:Nz
    if z(j) > z1
        sigma(j) = sigma_max * ((z(j) - z1) / (zend - z1))^3;
    end
end

% Derivative of sigma
sigma_prime = zeros(Nz, 1);
for j = 2:Nz-1
    sigma_prime(j) = (sigma(j+1) - sigma(j-1)) / (2 * dz);
end


% % Integral of sigma, Trapezoidal method
I_sigma = zeros(Nz, 1);
for j = 2:Nz
    if z(j) > z1  % Only accumulate in the PML region
        I_sigma(j) = I_sigma(j-1) + 0.5*(sigma(j) + sigma(j-1))* dz;
    else
        I_sigma(j) = 0;  % No stretching in the main computational domain
    end
end


% F(z) for PML region
F_z = z; 
for j = 1:Nz
    if z(j) > z1
        F_z(j) = z(j) + 1i * I_sigma(j);
    end
end
    
    % F(z) for PML region
    F_z = z + 1i * I_sigma;
    N_F = 315 * exp(-0.136 * (F_z) / 1000);
    n_F = 1 + N_F / 1e6;
    m_F = sqrt(n_F.^2 + 2 * (F_z) / ra);
    n_z = (wk^2 * (m_F.^2 - 1));
    
    % Construct matrices for Crank-Nicolson
    F = zeros(Nz, Nz);
   F = zeros(Nz, Nz);
for j = 2:Nz-1
    factor = 1 / (1 + 1i * sigma(j));
    F(j, j-1) = factor * 1/dz^2 + factor^2*1i*sigma_prime(j)/(2*dz);
    F(j, j)   = factor * -2/dz^2;
    F(j, j+1) = factor * 1/dz^2 - factor^2*1i*sigma_prime(j)/(2*dz);
end
    
    I = eye(Nz);
    A = sparse(I - (F + diag(n_z)) * dx * 1i / (2 * 2 * wk));
    B = sparse(I + (F + diag(n_z)) * dx * 1i / (2 * 2 * wk));

    % Change boundaries, dirichlet BC
    A(1,:) = 0; A(1,1) = 1;
    A(end,:) = 0; A(end,end) = 1;
    B(1,:) = 0; B(1,1) = 1;
    B(end,:) = 0; B(end,end) = 1;

    [L_mat, U_mat] = lu(A);
    
    for i = 1:Nx-1
        U(:, i+1) = U_mat \ (L_mat \ (B * U(:, i)));
    end
end




