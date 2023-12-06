clc; 
clear; 

f = @(x,t) (pi^2 - 1)*exp(-t)*sin(pi.*x); % given function -- scalar
method = 0; % forward euler = 1; backward euler = 0; 

% user-inputs 
N=11; % nodes 
n=551; % time step -> dt = 1/n

%left and right boundary = [0, 1] 
xi = linspace(0, 1, N); % may need to be N + 1
ts = linspace(0, 1,n+1);% ctime 
dt = 1/n; % delta time 

xr = 0; 
xl = 1; 
% h = xi(2)-xi(1); % uniform step size
h = (xr - xl) / (N - 1); % --- xr = 0; xl = 1; 
x_eta = @(eta, xi) (eta + 1)*(h/2) + xi; % equation given in pseudocode

% phi functions -- mapping to [-1, 1]: scalar 
phi_1 = @(eta)(1 - eta) / 2;
phi_2 = @(eta)(1 + eta) / 2; 

phi_deta = [-1/2, 1/2];

% more derivatives 
eta_dx = 2/h;
dx_deta = h/2;

% taken from  gaussian quadrature table: 
quad_1 = -1/sqrt(3);
quad_2 = 1/sqrt(3);
% weights = 1 

%initialize global mass, stiffness, and force matrices 
M = zeros(N,N);
K = zeros(N,N);
f_global = zeros([N n+1]); % turn into array of n+1 elements?

%initializing stiffness matrix
k_local = zeros(2,2);
f_local = zeros([1 n+1]);

%local to global mapping matrix
iee = [linspace(0,N-2,N-1); linspace(1,N-1,N-1)]'; % 10 x 2 

phi_eta = [phi_1(quad_1), phi_1(quad_2); phi_2(quad_1), phi_2(quad_2)]; % solutions

for k = 1:N-1 % looping over global grid elements 
    for l = 1:2 % calculating local elements  
        % mapping to parent function [-1, 1]
        % calculating local element 
        y1 = x_eta(quad_1, xi(k+l-1));
        % N = 10 & l = 2 -> 12 - 1 = 11; 
        y2 = x_eta(quad_2, xi(k+l-1));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f_local = (f(y1,ts)*phi_1(quad_1) + f(y2, ts)*phi_1(quad_2))*dx_deta;
        
        for m = 1:2 
            k_local(l,m) = phi_deta(l)*eta_dx*phi_deta(m)*eta_dx*dx_deta; 
        end
   
    end
    
    for l = 1:2 % finite element assembly 
        global_node = iee(k,l);
        global_node = global_node + 1;

        for m = 1:2
            global_node2 = iee(k,m);
            global_node2 = global_node2 + 1;
            K(global_node,global_node2) = K(global_node,global_node2) + k_local(l,m);
        end

        f_global(:, k+l-1) = f_global(:, k+l-1) + f_local(1, l);
    end
    
end

K(1,:) = 0;
K(:,1) = 0; 
K(1,1) = 1;

K(N,:) = 0;
K(:,N) = 0;
K(N,N) = 1;

%% 

inv_M = inv(M);

B = (1/dt)*M + K;
inv_B = inv(B);

u = zeros(N, n+1); % 1 x 552  

%initial conditions
f_boundary = @(x) sin(pi*x); % u(x,0)
u(1, 1) = f_boundary(1); % only need to use f_boundary for this??

%dirichlet boundary conditions
boundaries = [0,0];
dbc = eye(N); % identity matrix: 11 x 11 
dbc(1,1) = boundaries(1);
dbc(N,N) = boundaries(2);

%forward euler method is type = 1
if method == 1
for j = 1:n
    u(:, j+1) = u(:, j) - dt.*inv_M.*K*u(:, j) + dt.*inv_M*f_global(:, j);  
    % 11 x 1 = (11x1) - (11x11)(11x1) + (11x11)(11x1)
    
    % need to apply the boundary condition????
    u(:, j+1) = dbc*u(:, j+1);
end

else % backward euler method = 0
for j = 1:n
    u(:, j+1) = (1/dt).*inv_B*M*u(:, j) + inv_B*f_global(:, j+1); 

    % applying the boundary condition
    u(:, j+1) = dbc*u(:, j+1);
end
end

%% plotting 
xn = linspace(0,1,1000);
exact_sol = exp(-1)*sin(pi*xn); % exact solution 
plot(xn, exact_sol); 
hold on
x = linspace(0,1,N);
plot(x,u);
title('1D Finite Element Method');
subtitle('BEM - nodes: 11'); 
legend({'exact solution'}, 'Location', 'northeast'); 
