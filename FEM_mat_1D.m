function A = FEM_mat_1D(x,op1,op2,lin_fields)    
if nargin == 4
elseif nargin == 3
    lin_fields = [];
else
    error('Too many or too few input arguments')
end
x = x(:);
Nnodes = numel(x);
Nelem = Nnodes-1;
T10 = [1:(Nnodes-1);2:Nnodes]';
Lin_shp{1}= @(xi,Xe) [(1-xi)/2;(1+xi)/2];
Lin_shp{2} = @(xi,Xe) [-ones(1,numel(xi));ones(1,numel(xi))]/(Xe(2)-Xe(1)); 
Lin_eval = @(xi,Xe) (Xe(1)*(1-xi) + Xe(2)*(1+xi))/2;
Jac = @(Xe) (Xe(2)-Xe(1))/2;
xi_integ = [0.774596669241483 0.0 -0.774596669241483];
w_integ = [0.555555555555556 0.888888888888889  0.555555555555556];
A = sparse(Nnodes,Nnodes);
for i=1:Nelem
    map = T10(i,:);
    Xe = x(map);
    loc_jac = Jac(Xe);
    M1 = Lin_shp{op1+1}(xi_integ,Xe);
    M2 = Lin_shp{op2+1}(xi_integ,Xe);
    D = w_integ;
    for j=1:size(lin_fields,2)
        D = D.*Lin_eval(xi_integ,lin_fields(map,j));
    end
    A(map,map) = A(map,map)+(M1*diag(D)*M2')*loc_jac;
end
end





