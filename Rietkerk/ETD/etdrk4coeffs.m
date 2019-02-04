function [E, E2, Q,f1,f2,f3] = etdrk4coeffs(L,dT)


% Precompute various ETDRK4 scalar quantities:
A=L*dT;
A2=A.*A;
A2c=A2.*L;
E = exp(A); E2 = exp(A/2);



%%%%%%%%%%%%%%%%%%%%%%
%coefficients for stepping
M = 256; % no. of points for complex means
roots1 = exp(1i*pi*((1:M)-.5)/M); % roots of unity
%I = eye(Nx,); Z = zeros(N-1);

Q=zeros(size(A));
f1 = zeros(size(A)); 
f2 = zeros(size(A)); 
f3 = zeros(size(A)); 

%f1 = Z; f2 = Z; f3 = Z; Q = Z;
mn = find(abs(A)<=1e-1);
if ~(isempty(mn))
    zIA=zeros(size(A));
    for j = 1:M
        z = roots1(j);
        zIA(mn) = 1./(z-A(mn));
        Q(mn) = Q(mn) + dT*zIA(mn)*(exp(z/2)-1);
        f1(mn) = f1(mn) + dT*zIA(mn)*(-4-z+exp(z)*(4-3*z+z^2))/z^2;
        f2(mn) = f2(mn) + dT*zIA(mn)*(2+z+exp(z)*(z-2))/z^2;
        f3(mn) = f3(mn) + dT*zIA(mn)*(-4-3*z-z^2+exp(z)*(4-z))/z^2;
    end
    %LR = h*KSEL(:,ones(M,1)) + r(ones(N,1),:);
    %Q = h*real(mean( (exp(LR/2)-1)./LR ,2));
    %f1 = h*real(mean( (-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3 ,2));
    %f2 = h*real(mean( (2+LR+exp(LR).*(-2+LR))./LR.^3 ,2));
    %f3 = h*real(mean( (-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3 ,2));
    f1 = real(f1/M); f2 = real(f2/M); f3 = real(f3/M); Q = real(Q/M);
end


m = find(abs(A)>1e-1);
Q(m) = dT*(E2(m)-1)./A(m);
f1(m) = (-(4+A(m))+E(m).*(4-3*A(m)+A2(m)))./(A2c(m));
f2(m) = (2+A(m)+E(m).*(-2+A(m)))./(A2c(m));
f3(m) = (-(4+3*A(m)+A2(m))+E(m).*(4-A(m)))./(A2c(m));
