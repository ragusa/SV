% clear all; close all; 
clc;
syms u v c
syms nx ny

assume(c,'real')
assumeAlso(nx,'real')
assumeAlso(ny,'real')
assumeAlso(u,'real')
assumeAlso(v,'real')
assumeAlso(c>0)
assumeAlso(nx*nx+ny*ny==1)

% A = dF/dU
A = [ 0 1 0 ;...
    -u^2+c^2, 2*u, 0;...
    -u*v, v, u];
% B= dG/dU
B = [ 0 0 1 ;...
    -u*v, v, u; ...
    -v^2+c^2, 0, 2*v ];
% build K
K = nx*A +ny*B
% find eigenvalues and right eigenvectors
[R,VP]=eig(K);
VP=simplify(VP)
R=simplify(R);
R(:,1)=R(:,1)*nx/c/2; %%% jcr
R(:,2)=R(:,2)*(v+c*ny)/c/2;
R(:,3)=R(:,3)*(v-c*ny)/c/2;
R=simplify(R)
% check-1
aux=K*R-R*VP;
% assume(nx,'real')
% assume(ny,'real')
aux=simplify(aux);
assumeAlso((nx*nx+ny*ny)^(1/2)==1)
aux=simplify(aux)
%% R=R*2*c

% find eigenvalues and left eigenvectors
% assume(nx,'real')
% assume(ny,'real')
[Lt,VP]=eig(K');
assume(nx*nx+ny*ny==1)
VP=simplify(VP);
% assume(ny,'real')
% assume(nx,'real')
% assume(nx*nx+ny*ny==1)
Lt=simplify(Lt);
Lt(:,1)=Lt(:,1)*nx*2*c;
Lt(:,2)=Lt(:,2)*ny;
Lt(:,3)=Lt(:,3)*(-ny);
L=Lt';
% assume(nx,'real')
% assume(ny,'real')
L=simplify(L)
% L(:,1)=L(:,1)/c
% check-1
aux=L*K-VP*L;
% assume(nx,'real')
% assume(ny,'real')
aux=simplify(aux);
assume((nx*nx+ny*ny)^(1/2)==1)
aux=simplify(aux)


a=L*K*R;
assume(ny,'real')
assume(nx,'real')
a=simplify(a);
assumeAlso((nx^2 + ny^2) == 1);
a=simplify(a)



LR=L*R
% assume(ny,'real')
% assume(nx,'real')
LR=simplify(LR);
% assume(nx,'real')
% assume(ny,'real')
assumeAlso(nx*nx+ny*ny==1)
LR=simplify(LR);
% assume(nx,'real')
% assume(ny,'real')
LR=simplify(LR)

% LR=R*L
% assume(ny,'real')
% assume(nx,'real')
% LR=simplify(LR);
% assume(nx,'real')
% assume(ny,'real')
% assume(nx*nx+ny*ny==1)
% LR=simplify(LR);
% assume(nx,'real')
% assume(ny,'real')
% LR=simplify(LR)

simplify(L*A*R)
simplify(L*B*R)


return

syms h
assume(h,'real');
U=[h;h*u;h*v];
W=L*U;
W=simplify(W)

return


%%%%%%%%%%%
error('qqq')
%%%%%%%%%%%
clear u v nx ny c lambda

c=20;
u=5;
v=3;
nx=0.6;
ny=sqrt(1-nx*nx);


