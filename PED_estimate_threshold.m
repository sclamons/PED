function selection_size = PED_estimate_threshold(filename, ...
                                                 num_simulations, fdr)
    % PED_extimate_threshold
    % Last edited March 12, 2014
    % Samuel Clamons
    %
    % Determine the optimal threshold constant to use for a dataset.
    % Determines a threshold constant to maximize selection size while
    % keeping the average FDR (estimated from simulations based on the 
    % data) below a user-defined value.
    %
    % Parameters:
    %   filename        - Name of the data file holding data to estimate
    %                       against.
    %   num_simulations - Number of simulations to generate and test on.
    %                       More simulations may make FDR estimation more
    %                       accurate, but increases run time.
    %   fdr             - Highest acceptable false discovery rate (fraction
    %                       of selected genes that are not actually
    %                       differentially expressed).
    %
    % Return: The maximum selection size that still adequately controls the
    %           false discovery rate.

    %%%%%%%%%%%%%%%%%%%%%%
    % SET UP SIMULATIONS %
    %%%%%%%%%%%%%%%%%%%%%%
    
    sim_location     = PED_generate_simulations(num_simulations, 1.2, ...
                                                filename, 0);
    data_location    = fullfile(sim_location, 'data');
    simulation_files = dir(fullfile(data_location, '*.csv'));
    num_simulations  = length(simulation_files);
    
    %%%%%%%%%%%%%%%%%%%
    % PED FIRST PASS  %
    %%%%%%%%%%%%%%%%%%%
    
    % Run first weight calculation pass, storing them so we don't have to
    % recalculate them for each new threshold. Each simulation gets two
    % entries in simulation_weights -- one for the classifications
    % as-given (stored in odd-numbered entries), and one for the 
    % classifications reversed (stored in even-numbered entries).
    simulation_weights = cell(2 * num_simulations, 1);
    simulation_labels  = cell(num_simulations, 1);
    num_genes = -1;
    for i = 1:length(simulation_files)
        [labels, M, y] = open_microarray_file(...
                        fullfile(data_location, simulation_files(i).name));  
        num_genes = size(M, 1)
        simulation_weights{2*i - 1} = PED_calculate_weights(M', y);
        y = mod(y + 1, 2);
        simulation_weights{2*i} = PED_calculate_weights(M', y);
        simulation_labels{i}    = labels;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%
    % OPTIMIZE THRESHOLD %
    %%%%%%%%%%%%%%%%%%%%%%
    
    selection_size  = max(floor(0.001 * num_genes), 1);
    SIZE_DELTA = selection_size;
    % Start threshold search with a step size of 0.1%; once the best 
    % threshold has been found with that granularity, step back and
    % optimize again one gene at a time.
    rough_tuning = true;
    while selection_size <= num_genes
        empirical_fdrs = -1 * ones(num_simulations, 1);
       
        for sim_num = 1:num_simulations
            [labels, M, y] = open_microarray_file(fullfile(...
                          data_location, simulation_files(sim_num).name));
            M = M';
            
            % Need to do the following twice -- one for classifications
            % as-written, one for classifications reversed -- and combine
            % the results
            both_selections = [];
            for parity = -1:0
                % Cut the results of the first pass with current threshold.            
                curr_weights = simulation_weights{sim_num*2 + parity};
                [~, indices] = sort(abs(curr_weights), 'descend');
                             
                selections   = indices(1:selection_size);

                % Apply second pass, then cut the results with the current
                % threshold.
                y = mod(y + 1 + parity, 2);
                curr_weights     = ...
                           PED_calculate_weights_final(M(:,selections), y);
                % Remove the intercept component of curr_weights.       
                curr_weights     = curr_weights(1:end-1);       
                % Retrieve the indices of the original data corresponding 
                % to the final selections.
                final_indices    = find(abs(curr_weights) >= 10^-6); 
                final_selections = selections(final_indices); 
                both_selections  = union(both_selections, ...
                                         final_selections);  
            end
            
            % Calculate accuracy
            disp(['During permutation simulations, selections with threshold ', num2str(selection_size), ':'])
            for i = 1:length(both_selections)
                disp(['\t', labels(both_selections(i))])
            end
            null_genes = ...
                cell2mat(strfind(labels(both_selections), 'no_diff'));
            empirical_fdrs(sim_num) = ...
                              length(null_genes) /length(both_selections);
            
            disp(['With selection size = ', ...
                  num2str(selection_size),...
                  ', average false discovery rate is ', ...
                  num2str(empirical_fdrs(sim_num))]);                
        end        
        
        % Termination condition -- if the FDR is higher than the target
        % FDR, then step back and continue with fine-grained delta. The 
        % second time the FDR target is hit, step back and return.
        if max(empirical_fdrs) >= fdr
            selection_size = selection_size - SIZE_DELTA;
            if rough_tuning % First time FDR threshold exceeded
                SIZE_DELTA = 1;
                rough_tuning = false;
            else            % Second time FDR threshold exceeded  
                return;
            end
        end
        
        selection_size = selection_size + SIZE_DELTA;   
    end
    error('Reached maximum threshold without reaching desired FDR.');
end

    