clear all; close all; clc;
syms u v c
syms nx ny

assume(c,'real')
assumeAlso(ny,'real')
assumeAlso(nx,'real')
assumeAlso(u,'real')
assumeAlso(v,'real')
assumeAlso(c>0)
assumeAlso(nx*nx+ny*ny==1)

K=[0 nx ny;...
    (-u*u+c*c)*nx-u*v*ny, 2*u*nx+v*ny,        u*ny;...
    -u*v*nx+(-v*v+c*c)*ny,       v*nx, u*nx+2*v*ny]

syms lambda
assumeAlso(lambda,'real')
K=K-lambda*eye(3);

K1=K(2:3,:);
K1(:,1)=[];

K2=K(2:3,:);
K2(:,2)=[];

K3=K(2:3,:);
K3(:,3)=[];

poly = -lambda*det(K1)-nx*det(K2)+ny*det(K3)

simplify(poly)
% version_date=sprintf(version -date);
% if(datenum(a,'mmm dd, yyyy') > 73515)
%     simplify(poly,'IgnoreAnalyticConstraints',true,'Criterion', 'preferReal','Steps',100)
% end
check_version=struct2cell(ver);
if ~strcmp(check_version{3},'(R2012b)')
    simplify(poly,'IgnoreAnalyticConstraints',true,'Criterion', 'preferReal','Steps',100)
end

fprintf('\n\n');
%%%%%%%%%%%
syms h
assumeAlso(h,'real');

M=[ 1 0 0; u h 0; v 0 h]
iM=inv(M)

syms A B
A=[0  1  0;...
    c*c-u*u, 2*u, 0;...
    -u*v, v, u]
B=[0  0  1;...
    -u*v, v, u;...
    c*c-v*v, 0, 2*v]

% check
% K-(nx*A+ny*B)

% transformed matrices
fprintf('\nsimplify(Atilde):: \n');
Atilde = simplify(iM*A*M)
fprintf('\nsimplify(Btilde):: \n');
Btilde = simplify(iM*B*M)

Ktilde = Atilde*nx + Btilde*ny;
fprintf('\nsimplify(Ktilde):: \n');
simplify(Ktilde)

% get eigenvalues for Ktilde:
% assume(nx,'real')
% assume(ny,'real')
assumeAlso(nx^2 + ny^2 == 1);
assumeAlso((nx^2 + ny^2)^(1/2) == 1);
fprintf('\nsimplify(eig(Ktilde)):: \n');
simplify(eig(Ktilde))
[right,VP,left]=eig(Ktilde);
fprintf('\nsimplify(right):: \n');
right = simplify(right)
% fprintf('\nsimplify(VP):: \n');
% simplify(VP)
% quick check: one should have Ktilde*right = right*VP
fprintf('\nquick check: one should have Ktilde*right - right*VP = 0\n');
if ~strcmp(check_version{3},'(R2012b)')
    simplify(Ktilde*right - right*VP, 'IgnoreAnalyticConstraints',true,'Criterion', 'preferReal','Steps',100)
else
    simplify(Ktilde*right - right*VP)
end

% [Ri,VP]=eig(K)
% [Le,VP]=eig(K')
% 
% 
% assume(c,'real')
% assume(nx,'real')
% assume(ny,'real')
% assume(u,'real')
% assume(v,'real')
% assume((nx^2 + ny^2) == 1);
% Ri=simplify(Ri)
% 
% %%%%%%%%%%%
% error('qqq')
% %%%%%%%%%%%

% % syms L
% % assumeAlso(L,'real')
% % % assumeAlso(nx,'real')
% % % assumeAlso(ny,'real')
% % % old and wrong: L = [0 1 -1; 1 h/c*nx h/c*ny; -1 h/c*nx h/c*ny]
% % L = [0 ny -nx; c h*nx h*ny; -c h*nx h*ny]'
% % fprintf('\nsimplify(L):: \n');
% % simplify(L)

[left,VP]=eig(Ktilde');
fprintf('\nsimplify(left):: \n');
left = simplify(left)
fprintf('\nsimplify(VP from Ktilde transposed):: \n');
simplify(VP)

% save values before playing with them
right_sav=right;
left_sav=left;

right(:,1)=right(:,1)*nx;
right(:,2)=right(:,2)*c*ny;
right(:,3)=right(:,3)*c*ny*(-1);
% 
left(:,1)=left(:,1)*nx;
left(:,2)=left(:,2)*h*ny;
left(:,3)=left(:,3)*c*ny*(-1);
% 
right(:,2)=right(:,2)/c/h/2;
right(:,3)=right(:,3)/c/c/2;
fprintf('\nNew right:: \n');
right
fprintf('\nNew left'':: \n');
left'



LR =left'*right;
% assumeAlso(c,'real')
% assume(nx,'real')
% assume(ny,'real')
% assumeAlso(u,'real')
% assumeAlso(v,'real')
assumeAlso((nx^2 + ny^2) == 1);
LR=simplify(LR)

%%%% do an in book: The Development and Testing of Characteristic-based Semi-Lagrangian
right(:,1)=right(:,1)*(-1);
right(:,2)=right(:,2)*(c^2);
right(:,3)=right(:,3)*(-c^3/h);
fprintf('\nNew right:: \n');
syms g
assumeAlso(g,'real');
assumeAlso(g==c^2/h);
right=simplify(right)

fprintf('\nNew left'':: \n');
left(:,1)=left(:,1)*(-1);
left(:,2)=left(:,2)/(c*c);
left(:,3)=left(:,3)/(-c*c*c/h);
left=simplify(left); left'

LR =left'*right;
assumeAlso((nx^2 + ny^2) == 1);
LR=simplify(LR)

LAR=simplify(left'*Atilde*right)
LBR=simplify(left'*Btilde*right)

Cx=LAR-diag(diag(LAR))
Cy=LBR-diag(diag(LBR))
syms V W Vx Vy
assumeAlso(V,'real');
assumeAlso(Vx,'real');
assumeAlso(Vy,'real');
assumeAlso(W,'real');
V=[h;u;v];
W = left' * V;
W=simplify(W)

syms hx ux vx hy uy vy;
assumeAlso(hx,'real');
assumeAlso(ux,'real');
assumeAlso(vx,'real');
assumeAlso(hy,'real');
assumeAlso(uy,'real');
assumeAlso(vy,'real');
Vx=[hx;ux;vx];
Vy=[hy;uy;vy];
simplify(Cx*left'*Vx)
simplify(Cy*left'*Vy)

error('q')

fprintf('\nquick check: one should have Ktilde^T*left - left*VP = 0\n');
% assume(nx,'real')
% assume(ny,'real')
% assume(nx^2 + ny^2 == 1);
assumeAlso((nx^2 + ny^2)^(1/2) == 1);
simplify(left'*Ktilde-VP*left')
simplify(Ktilde'*left-left*VP)

%%%%%%%%%%%
error('qqq')
%%%%%%%%%%%
clear u v nx ny c lambda

c=20;
u=5;
v=3;
nx=0.6;
ny=sqrt(1-nx*nx);


%%%
lambda = u*nx+v*ny;

K=[0 nx ny;(-u*u+c*c)*nx-u*v*ny, 2*u*nx+v*ny,u*ny;-u*v*nx+(-v*v+c*c)*ny,v*nx,u*nx+2*v*ny];
K=K-lambda*eye(3);
K1=K(2:3,:); K1(:,1)=[];
K2=K(2:3,:); K2(:,2)=[];
K3=K(2:3,:); K3(:,3)=[];
poly = -lambda*det(K1)-nx*det(K2)+ny*det(K3)

%%%
lambda = u*nx+v*ny + c;

K=[0 nx ny;(-u*u+c*c)*nx-u*v*ny, 2*u*nx+v*ny,u*ny;-u*v*nx+(-v*v+c*c)*ny,v*nx,u*nx+2*v*ny];
K=K-lambda*eye(3);
K1=K(2:3,:); K1(:,1)=[];
K2=K(2:3,:); K2(:,2)=[];
K3=K(2:3,:); K3(:,3)=[];
poly = -lambda*det(K1)-nx*det(K2)+ny*det(K3)

%%%
lambda = u*nx+v*ny - c;

K=[0 nx ny;(-u*u+c*c)*nx-u*v*ny, 2*u*nx+v*ny,u*ny;-u*v*nx+(-v*v+c*c)*ny,v*nx,u*nx+2*v*ny];
K=K-lambda*eye(3);
K1=K(2:3,:); K1(:,1)=[];
K2=K(2:3,:); K2(:,2)=[];
K3=K(2:3,:); K3(:,3)=[];
poly = -lambda*det(K1)-nx*det(K2)+ny*det(K3)
