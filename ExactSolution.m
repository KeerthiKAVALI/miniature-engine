function [U] = ExactSolution(x,t,N)
k = 1;
X = zeros(numel(x),N);
T = zeros(numel(t),N);
L = x(end);
n = 1:2:(2*(N-1));
X =bsxfun(@(x,n) sin(n*pi*x/L),x(:),n);
T =bsxfun(@(t,n) exp(-k*n^2*pi^2*t/L^2),t(:),n);
A = -4*L^3/(k*pi^3)*n.^(-3);

U = bsxfun(@times,X,A)*T' +  bsxfun(@times,(x(:).*(L-x(:)))/2/k,ones(1,numel(t)));
end
