% shapeops.m
function [N,S] = inferNSfromTRADE(TRADE)
% TRADE is (N*S)×N or (N*S)×N×1 as in DATA.mat
    N = size(TRADE,2);
    S = size(TRADE,1) / N;
    assert(abs(S - round(S)) < 1e-9, 'TRADE first dim not divisible by N.');
    S = round(S);
end

function tf = is_bilat_cube(X), tf = (ndims(X)==3); end               % N×N×S
function tf = is_mfn_matrix(X,N), tf = isequal(size(X),[NaN N]); end % S×N (S free)

function C = to_bilat_cube(X, TARIFFs)
% Coerce any of: N×N×S, S×N, (N-1)S×N into N×N×S (zeros on diag)
    tgt = size(TARIFFs); N=tgt(1); S=tgt(3);
    if isequal(size(X), tgt), C = X; return; end
    if isequal(size(X), [S N])
        % MFN import tariffs: column j applies to all i≠j
        C = zeros(N,N,S);
        for j=1:N, col=reshape(X(:,j),[1 1 S]); % 1×1×S slice for replication
            C(:,j,:) = repmat(col,[N,1,1]); C(j,j,:)=0; end
        return
    end
    if isequal(size(X), [(N-1)*S, N])
        C = zeros(N,N,S);
        for j=1:N
            tmp = reshape(X(:,j), [N-1, 1, S]); % (N-1)×1×S vector without self-row
            C(:,j,:) = [tmp(1:j-1,:,:); zeros(1,1,S); tmp(j:end,:,:)];
        end
        return
    end
    error('to_bilat_cube: incompatible size %s', mat2str(size(X)));
end

function M = importer_stack_from_cube(C)
% N×N×S → (N-1)×N×S by removing the diagonal rows country-by-country
    [N,~,S] = size(C);
    M = zeros(N-1,N,S);
    for j=1:N
        M(:,j,:) = C([1:j-1, j+1:N], j, :);
    end
end

function M = mfn_matrix_from_cube(C)
% N×N×S → S×N MFN import tariffs by importer j (avg over i≠j)
    [N,~,S] = size(C);
    M = zeros(S,N);
    for j=1:N
        M(:,j) = reshape(mean(C([1:j-1, j+1:N], j, :),1), [S,1]); % average over exporters -> S×1
    end
end
