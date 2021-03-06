function weights = PED_calculate_weights_final(data, classifications)            
% PED_calculate_weights_final
% Last edited March 12, 2014
% Samuel Clamons and Daniel Vasiliu
% 
% Use LBFGS to calculate optimal weights for each gene, using the exact 
% objective function. 
%
% Parameters:
%   data            - Matrix of data values, in 'fat' format (genes in
%                       columns, samples in rows).
%   classifications - Vector of sample classifications. Classifications
%                       should be 0 or 1, i.e. [0,0,0,1,1,1].
%
% Return:
%   vector of optimized weights for every gene.

    %% Append an intercept variable to the end of the data
    [n_samples, n_genes] = size(data);
    data = [data, ones(n_samples, 1)];
    n_genes = n_genes + 1;

    %% Declare various constants
    X      = data;
    Y      = pi * (classifications-1/2);
    lambda = 1e-8;

    %% HANSO - BFGS settings
    pars.X=X;
    pars.Y=Y;
    pars.nvar = n_genes;
    pars.fgname = 'fgtest2';
    pars.lambda = lambda;

    options.x0 = ones(n_genes,1);
    options.normtol = 1e-4;
    options.evaldist = 1e-4;
    options.maxit = 1200;
    options.nvec = 100;
    
    %% Objective function to optimize
%     objective_function = ...
%         @(b)(...
%              norm(Y-atan(X*b), 2) + ...
%              lambda * (norm(b,1)*norm(b,2))^(1/2) ...
%             );
%     objective_gradient = ...
%         @(b)(...
%              (norm(Y-atan(X*b), 2)^(-1))* X' * ...
%              ((atan(X*b)-Y) ./ (1+(X*b).^2)) + ...
%              (lambda/2) * sign(b) * sqrt(norm(b,2)/norm(b,1)) + ...
%              (lambda/2) * (b/norm(b,2)) * sqrt(norm(b,1)/norm(b,2))...
%             );
%     function_and_gradient = @(b)fminunc_wrapper(b, ...
%                                                 objective_function, ...
%                                                 objective_gradient);

    %% Run LBFGS to optimize weights:
    weights = hanso(pars, options);  
end