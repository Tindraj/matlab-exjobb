% Solving Parabolic Wave equation with finite differences and 
% Crank Nicolson, without PML


tic

z0 = 0; 
zend = 10; 
Nz = 800;       % Number of points in z-direction
dz = (zend-z0)/(Nz-1);

Nx = 100;      % Number of steps in x-direction
x0 = 0;
xend = 10;
dx = (xend-x0)/(Nx-1);

lambda = 0.3;        % Radar wavelength
wk = 2*pi/lambda;    % Wave number
n = 1.000293;        % refractive index air

c = wk^2*(n^2-1);    % Constant n0


stability_limit = 1 / (2 * abs(wk));
dx_ratio = dx / dz^2;
disp(['dx / dz^2 = ', num2str(dx_ratio)]);
disp(['Stability limit: ', num2str(stability_limit)]);

a = 0;      % Dirichlet bundary condition at z0
b = 0;      % Dirichlet bundary condition at zN

z = linspace(z0,zend,Nz);
x = linspace(x0,xend,Nx);

% % Initial wave
% z_c = (z0+zend)/2;               % Center of the wave
% wave_par = 0.5;                           % Controls wave width
% u0 = 10*exp(-wave_par * (z - z_c).^2)';  % Gaussian wave as initial condition
% u0(1) = 0;
% u0(end) = 0;


% Initial condition from report
beta = 8*pi/180;    % Beam width parameter
theta0 = 0;         % Elevation angle in radians
ha = 5;            % Antenna height
k = wk;             % Wave number

% Compute initial field
u_fs = @(z) (k*beta)/(2*sqrt(2*pi*log(2))) .* exp(-1i*k*theta0*z) .* exp(-(beta^2/(8*log(2)))*k^2*(z-ha).^2);

u0 = u_fs(z) - u_fs(-z);

% Plot of initial wave
plot(z, abs(u0));
title('Initial Wave');
xlabel('z'); ylabel('Amplitude');

% Matrix from finite differences
F = diag(-2*ones(Nz,1)) + diag(ones(Nz-1,1),1) + diag(ones(Nz-1,1),-1);

F = F/dz^2;

% Change boundaries, dirichlet BC
F(1,:) = [1 zeros(1,Nz-1)];
F(Nz,:) = [zeros(1,Nz-1) 1];

I = eye(Nz);
% Construct system of equations A u_n+1 = B u_n, from Cranc Nicolson
A = I - (F + I*c)*dx*1i/(2*2*wk);
B = I + (F + I*c)*dx*1i/(2*2*wk);


U = [];
U(:,1) = u0;

% Iteration in x
for i = 1:Nx-1

    U(:,i+1) = A \ (B * U(:,i));
    U(1,1) = a;
    U(end,end) = b;

end


figure
surf(x,z,abs(U))
xlabel('x')
ylabel('z')
zlabel('u')
shading interp

%% Propagation factor

% Compute F_dB 
F_dB = 20 * log10(abs(U) + eps) + 10 * log10(x + eps) + 10 * log10(lambda);


figure;
surf(x, z, F_dB, 'EdgeColor', 'none');
xlabel('x (Propagation direction in km)');
ylabel('z (Height)');
zlabel('F_{dB}');
title('Propagation Factor in dB');
colorbar;
clim([-50 10]);
view(2); % 2D view for better visualization

toc
