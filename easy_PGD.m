function FF = easy_PGD(AA,BB,GG,Nmodes,Niter,Bords);

% compute some scalars
Ndims = size(AA,1);
Nterms = size(AA,2);
Nb = size(BB{1},2);
Nmode_0 = size(GG{1},2);

Nvars = zeros(1,Ndims);
for i = 1:Ndims
    Nvars(i) = size(AA{i,1},1);
end


FF = cell(Ndims,1);
% add GG to FF
for i = 1:Ndims
        FF{i} = GG{i};
end

% PGD loop

for mode = (Nmode_0+1):Nmodes
    % RS initialization
    RS = cell(1,Ndims);
    for i = 1:Ndims
        RS{i} = rand(Nvars(i),1);
        RS{i}(Bords{i}) = 0;
    end
    
    % fixed point iterations
    for iter = 1:Niter
        % update each RS
        for i = 1:Ndims
            mask = [ 1:i-1 i+1:Ndims];
            % construction of the system
            K = sparse(Nvars(i),Nvars(i));
            b = zeros(Nvars(i),1);
            for j = 1:Nterms
                factor_K = 1;
                for tmp_dim = mask
                    factor_K = factor_K * ((RS{tmp_dim}')*AA{tmp_dim,j}*RS{tmp_dim});
                end
                K = K + factor_K*AA{i,j};
            end
            
            % compute the right-hand side
            factors_b = ones(1,Nb);
            for tmp_dim = mask
                factors_b = factors_b .* ((RS{tmp_dim}')*BB{tmp_dim});
            end
            b = BB{i}*factors_b';
            
            if (mode>1)
                for j = 1:Nterms
                    factors_b = ones(1,mode-1);
                    for tmp_dim = mask
                        factors_b = factors_b .* ((RS{tmp_dim}')*AA{tmp_dim,j}*FF{tmp_dim});
                    end
                    b = b - (AA{i,j}*FF{i})*factors_b';
                end
            end
            
            % solve
            b(Bords{i}) = 0;
            dK = max(K(:));
            K(Bords{i},:) = 0;
            K(Bords{i},Bords{i}) = dK*speye(length(Bords{i}));
            RS{i}= K\b;
        end
    end
    % update of FF
    for i = 1:Ndims
        FF{i} = [FF{i} RS{i}];
    end
end
end
