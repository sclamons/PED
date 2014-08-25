function weights = PED_select_genes_from_threshold(data, ...
                                                   classifications, ...
                                                   threshold, ...
                                                   threshold_type)
% PED_select_genes_from_threshold
% Last edited March 12, 2013
% Samuel Clamons
% 
% Select genes using PED regression, given either a selection size or a 
% weight threshold.
% 
% Parameters:
%   data                - Matrix of data values, in 'fat' format (genes 
%                           in columns, samples in rows).
%   classifications     - Vector of sample classifications. 
%                           Classifications should be 0 or 1, i.e. 
%                           [0,0,0,1,1,1].
%   threshold           - Constant determining how stringently genes
%                           will be filtered out. The higher the
%                           size_threshold, the larger the
%                           selection size and the higher the expected
%                           FDR. This threshold can either be a number
%                           of genes to select or a minimum allowed
%                           (absolute value of) weight threshold.
%   threshold_type      - String argument, either 'size' or 'weight'.
%                           Determines what kind of threshold was
%                           given.
    
    % Basic matrix size information
    [n_samples, n_genes] = size(data);
    
    % Run PED twice -- once on classifications as given, once on
    % classifications reversed.
    combined_weights = zeros(n_genes, 1);
    for parity = 0:1
        % Reverse, or don't, the classifications.
        classifications = mod(classifications + parity, 2);
        
        %%%%%%%%%%%%%%
        % FIRST PASS %
        %%%%%%%%%%%%%%

        %% Optimize weights with LBFGS
        weights = PED_calculate_weights(data, classifications);
        
        if strcmpi(threshold_type, 'size')
            % Select only the top genes.
            [~, sort_indices] = sort(abs(weights), 'descend');
            remaining_genes = sort_indices(1:threshold);
        elseif strcmpi(threshold_type, 'weight')
            % Select only genes with weights above a threshold.
            remaining_genes = find(abs(weights)/norm(weights) > ...
                                threshold);
        else
            error(['Argument threshold_type has invalid value ',...
                   threshold_type]);
        end
        

        %%%%%%%%%%%%%%%
        % SECOND PASS %
        %%%%%%%%%%%%%%%

        %% Optimize and filter again to clean up
        weights = PED_calculate_weights_final(data(:, remaining_genes), ...
                                              classifications);
        weights = weights(1:end-1);
        for i = 1:length(remaining_genes)
            % i = index in returned weights
            % remaining_genes(i) = index of the same gene in original data.
            if abs(weights(i)) >= 10^-6
                if abs(weights(i)) > ...
                   abs(combined_weights(remaining_genes(i)))
                    combined_weights(remaining_genes(i)) = weights(i);
                end                           
            end
        end
    end    
    
    % Return
    weights = combined_weights;
end