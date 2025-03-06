% Solving Parabolic Wave equation with finite differences and 
% Crank Nicolson, with PML

% Works well, replicates what is done in the report.
% Flat earth, constant refraction index n(z) = k^2(n^2-1))

clear all, close all, clc
tic

z0 = 0; 
d = 40;                % Thickness of PML
sigma_max = 20;        % Maximum damping
Nz = 934;%5490;%1120;             % Number of points in z-direction
% Nyqvist stability: Ny = 5 490
zend = 140;
z1 = zend-d;           % Transition between domain and PML
z0 = 0;
dz = (zend-z0)/(Nz-1); % Compute spatial step

x0 = 0;                             
xend = 10000;
Nx = 66667;%12799;             % Number of steps in x-direction
% Nyqvist stability: Nx = 762 300, kommer ta 10 h
dx = (xend-x0)/(Nx-1);

lambda = 0.3;        % Radar wavelength
wk = 2*pi/lambda;    % Wave number
n = 1.000293;        % refractive index air

c = wk^2*(n^2-1);    % Constant n0 

a = 0;               % Dirichlet bundary condition at z0
b = 0;               % Dirichlet bundary condition at zN

z = linspace(z0,zend,Nz);
x = linspace(x0,xend,Nx);

% % Gaussian wave as initial condition
% z_c = (z0+z1)/2;                          % Center of the wave
% wave_par = 5;                             % Controls wave width, low -> wide
% u0 = 10*exp(-wave_par * (z - z_c).^2)';   % Gaussian wave as initial condition

% Initial condition from report
beta = 8*pi/180;    % Beam width parameter
theta0 = 0;         % Elevation angle in radians
ha = 50;            % Antenna height
k = wk;             % Wave number

% Compute initial field
u_fs = @(z) (k*beta)/(2*sqrt(2*pi*log(2))) .* exp(-1i*k*theta0*z) .* exp(-(beta^2/(8*log(2)))*k^2*(z-ha).^2);

u0 = u_fs(z) - u_fs(-z);

% Plot of initial wave
plot(z, abs(u0), 'r');
title('Initial Wave');
xlabel('z'); ylabel('Amplitude');



% Absorption function sigma(z)
sigma = zeros(Nz, 1);
% for j = 1:Nz
%     if z(j) > z1
%         sigma(j) = sigma_max * ((z(j) - z1) / (zend - z1))^3;
%     end
% end
sigma(z > z1) = sigma_max * ((z(z > z1) - z1) / (zend - z1)).^3;


% Derivative of sigma
sigma_prime = zeros(Nz, 1);
% for j = 2:Nz-1
%     sigma_prime(j) = (sigma(j+1) - sigma(j-1)) / (2 * dz);
% end
sigma_prime(2:Nz-1) = diff(sigma, 2) / (2 * dz);


%Matrix from finite differences
%Construct F matrix with PML
F = zeros(Nz, Nz);
for j = 2:Nz-1
    factor = 1 / (1 + 1i * sigma(j));
    F(j, j-1) = factor * 1/dz^2 + factor^2*1i*sigma_prime(j)/(2*dz);
    F(j, j)   = factor * -2/dz^2;
    F(j, j+1) = factor * 1/dz^2 - factor^2*1i*sigma_prime(j)/(2*dz);
end

%F = spdiags([factor(2:Nz) ./ dz^2, -2 * factor ./ dz^2, factor(1:Nz-1) ./ dz^2], [-1 0 1], Nz, Nz);


% Dirichlet BC, u0 = 0, uN = 0
F(1, :) = 0; F(1, 1) = 1;
F(Nz, :) = 0; F(Nz, Nz) = 1; 

% Construct system of equations A u_n+1 = B u_n, from Cranc Nicolson
I = eye(Nz);

A = sparse(I - (F + I*c)*dx*1i/(2*2*wk));
B = sparse(I + (F + I*c)*dx*1i/(2*2*wk));


[L, U_mat] = lu(A);  % Perform LU decomposition outside the loop

U = zeros(Nz, Nx); 
U(:,1) = u0;

% Iterate in x
for i = 1:Nx-1

    %U(:,i+1) = A \ (B * U(:,i));
    U(:,i+1) = U_mat \ (L \ (B * U(:,i)));

end

%%
% % Simulation
% figure
% for step = 1:10:Nx
%     plot(z, real(U(:,step)), 'b', z, abs(U(:,step)), 'r--');
%     title(['Step ' num2str(step)]);
%     xlabel('z'); ylabel('Amplitude');
%     legend('Re(u)', '|u|');
%     drawnow;
%     pause(0.5)
% end

%%


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


% Find indices corresponding to z values from 0 to 100
z_limit = 100;
z_indices = z <= z_limit;

% Extract the portion of F_dB and z that we want to plot
z_plot = z(z_indices);
F_dB_plot = F_dB(z_indices, :);

% Plot results excluding the PML region
figure;
surf(x, z_plot, F_dB_plot, 'EdgeColor', 'none');
xlabel('x (Propagation direction)');
ylabel('z (Height)');
zlabel('F_{dB}');
title('Propagation Factor in dB');
colorbar;
clim([-50 10]);
view(2); % 2D view for better visualization



% 3D plot, abs
figure
surf(x,z,abs(U));
title('Wave Propagation |u|');
xlabel('x'); ylabel('z'); zlabel('|u|');
shading interp


% Visualization of damping function
figure
plot(z, sigma);
title('PML Damping Function \sigma(z)');
xlabel('z'); ylabel('\sigma(z)');



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