% Solving Parabolic Wave equation with finite differences and
% Crank Nicolson, with PML

% Works well, replicates what is done in the report.
% Implementing varying refractive index and curved earth

clc
tic

z0 = 0;
d = 40;                % Thickness of PML
sigma_max = 20;        % Maximum damping
Nz = 1000;%14080;%17564;            % Number of points in z-direction
% Nyqvist criterion: Nz = 17 564
% Same as in report, dz = 1/32, Nz = 14080
zend = 140;
z1 = zend-d;           % Transition between domain and PML
z0 = 0;
dz = (zend-z0)/(Nz-1); % Compute spatial step

x0 = 0;
xend = 10000;
Nx = 1000;            % Number of steps in x-direction
% Nyqvist criterion: Nx = 3 745 300,
% Same as in report, dx = 50, Nx = 2 000,


dx = (xend-x0)/(Nx-1);

z = linspace(z0,zend,Nz);
x = linspace(x0,xend,Nx);

lambda = 0.1;        % Radar wavelength, 3 ghz = 0.1
wk = 2*pi/lambda;    % Wave number
n = 1.000293;        % varying refractive index
ra = 6360000;        % Radius of the earth
%n_z = wk^2*(n^2-1+2*z/ra)';    % n(z)

% Nyqvist criterion
alpha = atan(5000/319);
deltaz = lambda/(4*sin(alpha));
deltax = lambda/(8*sin(0.5*alpha)^2);
minNz = zend/deltaz;
minNx = xend/deltax;


a = 0;               % Dirichlet bundary condition at z0
b = 0;               % Dirichlet bundary condition at zN

% % Gaussian wave as initial condition
% z_c = (z0+z1)/2;                          % Center of the wave
% wave_par = 5;                             % Controls wave width, low -> wide
% u0 = 10*exp(-wave_par * (z - z_c).^2)';   % Gaussian wave as initial condition

% Initial condition from report
beta = 3*pi/180;    % Beam width parameter
theta0 = 0;         % Elevation angle in radians
ha = 31;            % Antenna height
k = wk;             % Wave number

% Compute initial field
u_fs = @(z) (k*beta)/(2*sqrt(2*pi*log(2))) .* exp(-1i*k*theta0*z) .* exp(-(beta^2/(8*log(2)))*k^2*(z-ha).^2);

u0 = u_fs(z) - u_fs(-z);

% % Plot of initial wave
% plot(z, abs(u0), 'r');
% title('Initial Wave');
% xlabel('z'); ylabel('Amplitude');



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


% Refractive index, standard atmosphere
N_F = 315 * exp(-0.136 * (F_z)/1000);
n_F = 1 + N_F / 10^6;
figure
plot(abs(N_F),z)
xlabel('M')
ylabel('z')
title('Standard atmosphere')

% % Refractive index, evaporation duct
% M0 = 315;
% delta = 40;
% z0 = 1.5e-4;
% N_F = M0 + 0.125*z - 0.125*delta*log((z+z0)/z0);
% n_F = 1 + N_F / 10^6;
% figure
% plot(N_F,z)
% xlabel('M')
% ylabel('z')
% title('Evaporation duct')


% Refractive index, elevated duct



% Refractive index, curved earth and n_F
m_F = sqrt(n_F.^2 + 2*(F_z)/ra);
n_z = (wk^2 * (m_F.^2 - 1));

% Matrix from finite differences
% F matrix with PML
F = zeros(Nz, Nz);
for j = 2:Nz-1
    factor = 1 / (1 + 1i * sigma(j));
    F(j, j-1) = factor * 1/dz^2 + factor^2*1i*sigma_prime(j)/(2*dz);
    F(j, j)   = factor * -2/dz^2;
    F(j, j+1) = factor * 1/dz^2 - factor^2*1i*sigma_prime(j)/(2*dz);
end


% Construct system of equations A u_n+1 = B u_n, from Cranc Nicolson
I = eye(Nz);

A = sparse(I - (F + diag(n_z))*dx*1i/(2*2*wk));
B = sparse(I + (F + diag(n_z))*dx*1i/(2*2*wk));

% Change boundaries, dirichlet BC
A(1,:) = 0; A(1,1) = 1;
A(end,:) = 0; A(end,end) = 1;
B(1,:) = 0; B(1,1) = 1;
B(end,:) = 0; B(end,end) = 1;


% Solving A u_n+1 = B u_n by first LU-factorizing, makes loop more
% efficient

[L_mat, U_mat] = lu(A);  % LU decomposition

U = sparse(zeros(Nz, Nx));
U(:,1) = u0;

toc

tic

% Iterate in x
for i = 1:Nx-1

    %U(:,i+1) = A \ (B * U(:,i));
    U(:,i+1) = U_mat \ (L_mat \ (B * U(:,i)));

end

toc

tic

% Propagation factor

F_dB = 20 * log10(abs(U) + eps) + 10 * log10(x + eps) + 10 * log10(lambda);

% Plot results
figure;
surf(x, z, F_dB, 'EdgeColor', 'none');
xlabel('x (Propagation direction)');
ylabel('z (Height)');
zlabel('F_{dB}');
title('Propagation Factor in dB');
colorbar;
clim([-50 10]);
view(2); % 2D view for better visualization


% Plot domain
% Find indices corresponding to z values from 0 to z_limit
z_limit = 350;
z_indices = z <= z_limit;

% F_dB withuot PML region
z_plot = z(z_indices);
F_dB_plot = F_dB(z_indices, :);

% Plot results without the PML region
figure;
surf(x, z_plot, F_dB_plot, 'EdgeColor', 'none');
xlabel('x (Propagation direction)');
ylabel('z (Height)');
zlabel('F_{dB}');
title('Propagation Factor in dB');
colorbar;
clim([-50 10]);
view(2); % 2D view


% % 3D plot, abs
% figure
% surf(x,z,abs(U));
% title('Wave Propagation |u|');
% xlabel('x'); ylabel('z'); zlabel('|u|');
% shading interp


% % Visualization of damping function
% figure
% plot(z, sigma);
% title('PML Damping Function \sigma(z)');
% xlabel('z'); ylabel('\sigma(z)');


% % Check energy evolution
% energy = zeros(1, Nx);
% for i = 1:Nx
%     energy(i) = sum(abs(U(:, i)).^2);
% end
% figure;
% plot(1:Nx, energy);
% title('Energy Over Time');
% xlabel('Time Step'); ylabel('Energy');

toc