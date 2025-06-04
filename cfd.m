clear; clc;

%% Input Parameters
nx = 50; ny = 50;
Lx = 1.0; Ly = 1.0;
rho = 1.0;
nu = 0.1;
dt = 0.001;
nsteps = 100;
u_top = 1.0; u_bot = 0.0;
v_left = 0.0; v_right = 0.0;

%% Index Extents and Mesh
imin=2; imax=imin+nx-1;
jmin=2; jmax=jmin+ny-1;
x(imin:imax+1)=linspace(0,Lx,nx+1);
y(jmin:jmax+1)=linspace(0,Ly,ny+1);
xm(imin:imax)=0.5*(x(imin:imax)+x(imin+1:imax+1));
ym(jmin:jmax)=0.5*(y(jmin:jmax)+y(jmin+1:jmax+1));
dx = x(imin+1)-x(imin);
dy = y(jmin+1)-y(jmin);
dxi = 1/dx; dyi = 1/dy;

%% Initialize Fields
u = zeros(imax+1,jmax+2);
v = zeros(imax+2,jmax+1);
us = u; vs = v;
p = zeros(imax,jmax);

%% Create Laplacian Operator
L = zeros(nx*ny,nx*ny);
for j = 1:ny
    for i = 1:nx
        row = i + (j-1)*nx;
        L(row,row) = 2*dxi^2 + 2*dyi^2;
        for ii = i-1:2:i+1
            if ii > 0 && ii <= nx
                col = ii + (j-1)*nx;
                L(row,col) = L(row,col) - dxi^2;
            else
                L(row,row) = L(row,row) - dxi^2;
            end
        end
        for jj = j-1:2:j+1
            if jj > 0 && jj <= ny
                col = i + (jj-1)*nx;
                L(row,col) = L(row,col) - dyi^2;
            else
                L(row,row) = L(row,row) - dyi^2;
            end
        end
    end
end
L(1,:) = 0; L(1,1) = 1;

%% Time Stepping
for step = 1:nsteps
    %% Predictor Step
    for j=jmin:jmax
        for i=imin+1:imax
            v_here = 0.25*(v(i-1,j)+v(i-1,j+1)+v(i,j)+v(i,j+1));
            us(i,j) = u(i,j) + dt * (...
                nu * ( (u(i-1,j)-2*u(i,j)+u(i+1,j)) * dxi^2 + ...
                       (u(i,j-1)-2*u(i,j)+u(i,j+1)) * dyi^2 ) - ...
                u(i,j) * (u(i+1,j)-u(i-1,j)) * 0.5 * dxi - ...
                v_here * (u(i,j+1)-u(i,j-1)) * 0.5 * dyi );
        end
    end
    for j=jmin+1:jmax
        for i=imin:imax
            u_here = 0.25*(u(i,j-1)+u(i,j)+u(i+1,j-1)+u(i+1,j));
            vs(i,j) = v(i,j) + dt * (...
                nu * ( (v(i-1,j)-2*v(i,j)+v(i+1,j)) * dxi^2 + ...
                       (v(i,j-1)-2*v(i,j)+v(i,j+1)) * dyi^2 ) - ...
                u_here * (v(i+1,j)-v(i-1,j)) * 0.5 * dxi - ...
                v(i,j) * (v(i,j+1)-v(i,j-1)) * 0.5 * dyi );
        end
    end

    %% Boundary Conditions on u* and v*
    u(:,jmin-1) = 2*u_bot - u(:,jmin);
    u(:,jmax+1) = 2*u_top - u(:,jmax);
    v(imin-1,:) = 2*v_left - v(imin,:);
    v(imax+1,:) = 2*v_right - v(imax,:);

    %% Poisson RHS
    R = zeros(nx*ny,1); n = 0;
    for j=jmin:jmax
        for i=imin:imax
            n = n + 1;
            R(n) = -rho/dt * ( ...
                (us(i+1,j)-us(i,j))*dxi + ...
                (vs(i,j+1)-vs(i,j))*dyi );
        end
    end

    %% Solve Pressure Poisson Equation
    pv = L\R;

    %% Map Pressure Vector to Matrix
    n = 0;
    for j=jmin:jmax
        for i=imin:imax
            n = n + 1;
            p(i,j) = pv(n);
        end
    end

    %% Corrector Step
    for j=jmin:jmax
        for i=imin+1:imax
            u(i,j) = us(i,j) - dt/rho * (p(i,j)-p(i-1,j)) * dxi;
        end
    end
    for j=jmin+1:jmax
        for i=imin:imax
            v(i,j) = vs(i,j) - dt/rho * (p(i,j)-p(i,j-1)) * dyi;
        end
    end
end