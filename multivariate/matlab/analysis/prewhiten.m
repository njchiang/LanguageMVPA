function [X, W] = prewhiten(X)
% apply prewhitening
    % Compute and apply the ZCA mapping
    X = X - repmat(mean(X, 1), [size(X, 1) 1]);
    W = inv(sqrtm(cov(X)));
    X = X * W;   
end