function f = mesh_lap_int(L,knownidx,f2)
    % solve the linear system belonging to the laplace equation
    % input:
    % L: laplace matrix
    % knownidx: index of points with known values (boundary values)
    % f2: vector with known values
    % output:
    % f: knitting time as a result of the laplace equation
    
N = size(L,1);
[KnownIndex,sortidx] = sort(knownidx);
f2 = f2(sortidx);
UnknownIndex = setdiff((1:N)', KnownIndex);
K = length(UnknownIndex);
% reshuffle rows & columns of lap matrix
lapi = [UnknownIndex; KnownIndex];
L = L(lapi, :); % rows
L = L(:, lapi); % columns

% Segregate known/unknown portions of lap
L11 = L(1:K    ,1:K    );
L12 = L(1:K    ,(K+1):N);
%L21 = L((K+1):N,1:K    );
%L22 = L((K+1):N,(K+1):N);

%f1 = [L11; L21]\(-[L12;L22]*f2);
f1 = L11\(-L12*f2);
f = [f1;f2];

[~, order] = sort(lapi);
f = f(order);


