clear all; close all; clc;
syms u v c
syms nx ny

assume(c,'real')
assume(ny,'real')
assume(nx,'real')
assume(u,'real')
assume(v,'real')
assume(c>0)
assume(nx*nx+ny*ny==1)

K=[0 nx ny;...
    (-u*u+c*c)*nx-u*v*ny, 2*u*nx+v*ny,        u*ny;...
    -u*v*nx+(-v*v+c*c)*ny,       v*nx, u*nx+2*v*ny]

syms lambda
assume(lambda,'real')
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
B=struct2cell(ver);
if ~strcmp(B{3},'(R2012b)')
    simplify(poly,'IgnoreAnalyticConstraints',true,'Criterion', 'preferReal','Steps',100)
end

fprintf('\n\n');
%%%%%%%%%%%
syms h

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
assume(nx,'real')
assume(ny,'real')
assume(nx^2 + ny^2 == 1);
assume((nx^2 + ny^2)^(1/2) == 1);
fprintf('\nsimplify(eig(Ktilde)):: \n');
simplify(eig(Ktilde))
[right,VP,left]=eig(Ktilde)
fprintf('\nsimplify(right):: \n');
simplify(right)
fprintf('\nsimplify(VP):: \n');
simplify(VP)
% simplify(left)
% quick check: one should have Ktilde*right = right*VP
fprintf('\nquick check: one should have Ktilde*right = right*VP \n');
simplify(Ktilde*right - right*VP, 'IgnoreAnalyticConstraints',true,'Criterion', 'preferReal','Steps',100)

%%%%%%%%%%%
error('qqq')
%%%%%%%%%%%

syms h L
assume(h,'real')
assume(L,'real')
assume(nx,'real')
assume(ny,'real')
L = [0 1 -1; 1 h/c*nx h/c*ny; -1 h/c*nx h/c*ny]
simplify(L')

simplify(L*Ktilde)
simplify(L'*Ktilde)

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
