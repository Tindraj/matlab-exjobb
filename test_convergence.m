% Solving Parabolic Wave equation with finite differences and
% Crank Nicolson, without PML
% test convergence

dx = 1;
dz = 0.1;

z0 = 0;
zend = 100;

x0 = 0;
xend = 10000;

z = z0:dz:zend;
x = x0:dx:xend;

Nz = length(z);
Nx = length(x);

tic
[U_coarse] = wavefunc(dx,dz,xend,zend);
% figure
% surf(abs(U_coarse))
% shading interp
toc
drawnow
tic
[U_fine] = wavefunc(dx/2,dz/2,xend,zend);
% figure
% surf(abs(U_fine))
% shading interp
toc
drawnow
tic
[U_finest] = wavefunc(dx/4,dz/4,xend,zend);
% figure
% surf(abs(U_finest))
% shading interp
toc
drawnow
tic


U_coarse_interp = U_fine(1:2:end, 1:2:end);  % Take every second point
error_1 =  norm(U_coarse - U_coarse_interp)/sqrt(Nx*Nz);

U_fine_interp = U_finest(1:2:end, 1:2:end);
error_2 =  norm(U_fine - U_fine_interp)/sqrt(2*Nx*2*Nz);

order = log(error_1 / error_2) / log(2);
disp(['Computed order of accuracy: ', num2str(order)]);
toc

function [U] = wavefunc(dx,dz,xend,zend)

z = 0:dz:zend;
x = 0:dx:xend;

Nz = length(z);
Nx = length(x);

a = 0;      % Dirichlet bundary condition at z0
b = 0;      % Dirichlet bundary condition at zN

lambda = 0.3;        % Radar wavelength,
wk = 2*pi/lambda;    % Wave number
n = 1.000293;        % refractive index

% Nyqvist criterion
alpha = atan(50/250);
deltaz = lambda/(4*sin(alpha));
deltax = lambda/(8*sin(0.5*alpha)^2);
minNz = zend/deltaz;
minNx = xend/deltax;

c = wk^2*(n^2-1);    % Constant n0

% % Initial wave gaussian
% z_c = (z0+zend)/2;               % Center of the wave
% wave_par = 0.5;                           % Controls wave width
% u0 = 10*exp(-wave_par * (z - z_c).^2)';  % Gaussian wave as initial condition
% u0(1) = 0;
% u0(end) = 0;

% Initial condition from report
beta = 8*pi/180;    % Beam width parameter
theta0 = 0;         % Elevation angle in radians
ha = 50;            % Antenna height
k = wk;             % Wave number

% Compute initial field
u_fs = @(z) (k*beta)/(2*sqrt(2*pi*log(2))) .* exp(-1i*k*theta0*z) .* exp(-(beta^2/(8*log(2)))*k^2*(z-ha).^2);

u0 = u_fs(z) - u_fs(-z);

% % Plot of initial wave
% plot(z, abs(u0));
% title('Initial Wave');
% xlabel('z'); ylabel('Amplitude');

% Matrix from finite differences
F = diag(-2*ones(Nz,1)) + diag(ones(Nz-1,1),1) + diag(ones(Nz-1,1),-1);

F = F/dz^2;

I = eye(Nz);
% Construct system of equations A u_n+1 = B u_n, from Cranc Nicolson
A = I - (F + I*c)*dx*1i/(4*wk);
B = I + (F + I*c)*dx*1i/(4*wk);

% Change boundaries, dirichlet BC
A(1,:) = 0; A(1,1) = 1;
A(end,:) = 0; A(end,end) = 1;
B(1,:) = 0; B(1,1) = 1;
B(end,:) = 0; B(end,end) = 1;

[L_mat, U_mat] = lu(A);  % LU decomposition

U = [];
U(:,1) = u0;

% Iteration in x
for i = 1:Nx - 1

    U(:,i+1) = U_mat \ (L_mat \ (B * U(:,i)));

end



end


