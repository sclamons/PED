function permuted_classifications = ...
                             PED_permute_classifications(Y, n_permutations)
    % PED_permute_classifications
    % Last edited March 12, 2014
    % Samuel Clamons
    %
    % Computes permutations of a classification vector of 0s and 1s. 
    % 
    % This function uses different algorithms for calculating permutations
    % for different numbers of permutations. For small classification
    % vectors, uses an algorithm based on generating all possible
    % combinations of samples from each condition to avoid
    % double-selection. If the classification vector is large, instead
    % randomly picks permutations.
    % 
    % Parameters:
    %   Y               - Original vector of classifications for an
    %                       experiment. Must consist entirely of 0s and 1s.
    %   n_permutations  - The number of permutations to create.
    %
    % Returns: A matrix wherein each column is a unique permutation of the
    %   classification vector.
    
    n_samples = length(Y);
    for i = 1:n_samples
        if Y(i) ~= 0 && Y(i) ~= 1
            error(['in permute_classifications: Y(%d)=%d, should be',...
                  ' 0 or 1'], i, Y(i));
        end
    end
    group_1 = find(Y);
    group_2 = setdiff(1:n_samples, group_1);
    n_swapped = min(...
                    floor(length(group_1)/2), ...
                    floor(length(group_2)/2)...
                    );
    
    % For relatively small numbers of samples, list out every possible
    % combination of samples from each condition, then choose two to swap
    % to create each permutation.
    if n_samples < 12
        % List all of the permutations (technically combinations) of
        % samples from each group. Each permutation is a swappable set of
        % samples. Choosing two such sets and swapping their vales creates
        % a permuted experiment.
        all_permutations_1 = nchoosek(group_1, floor(length(group_1)/2));
        all_permutations_2 = nchoosek(group_2, floor(length(group_2)/2));       
        
        % Encode each pair of sets of samples to swap as a base-D index,
        % where D is the number of samples swapped, the least significant 
        % digit of the index is the index of the set in the first group to 
        % swap, and the most significant digit of the index is the index of
        % the set in the second group to swap.
        %
        % Pick 'n_permutations' distinct permutations. For each
        % permutation, swap the corresponding indices of Y.
        base_1 = length(all_permutations_1);
        base_2 = length(all_permutations_2);
        n_permutations = min(base_1 * base_2, n_permutations);
        permutations = randperm(base_1 * base_2, n_permutations);
        
        permuted_classifications = zeros(n_samples, n_permutations);
        for i = 1:n_permutations
            % Indices run from 0 to base^2-1, but arrays run from 1
            % to base^2. Thus the egregious adding of 1.
            coded_indices = permutations(i);
            index_1 = mod(coded_indices, base_1) + 1;
            coded_indices = floor((coded_indices - index_1 + 1) / base_2);
            index_2 = mod(coded_indices, base_2) + 1;
            
            swap_indices_1 = all_permutations_1(index_1,:);
            swap_indices_2 = all_permutations_2(index_2,:);
            
            new_classifications = Y;
            new_classifications(swap_indices_1) = Y(swap_indices_2);
            new_classifications(swap_indices_2) = Y(swap_indices_1);
            permuted_classifications(:, i) = new_classifications';
        end
    else
    % For larger sample sizes, can't reasonably enumerate every
    % possibility. Instead, pick n_permutations at random -- at this point,
    % it is very unlikely that any identical permutations will be picked
    % unless you are making a LOT of permuted classifications (on the order
    % of (n_swapped!)).
        permuted_classifications = zeros(n_permutations);
        for i = 1:n_permutations
            swap_indices_1 = group_1(randperm(length(group_1), n_swapped));
            swap_indices_2 = group_2(randperm(length(group_2), n_swapped));
            new_classifications = Y;
            new_classifications(swap_indices_1) = Y(swap_indices_2);
            new_classifications(swap_indices_2) = Y(swap_indices_1);
            permuted_classifications(:, i) = new_classifications';
        end
    end
end