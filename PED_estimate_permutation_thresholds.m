function model_sizes = PED_estimate_permutation_thresholds(...
                                 input_filename, n_permutations, ...
                                 num_simulations, fdr)
    % PED_detect_expression
    % Last edited March 31, 2014
    % Samuel Clamons
    % 
    % Determine whether or not there is differential expression in a
    % dataset, using a permutation test on a data set with a
    % previously-determied estimated number of differentially-expressed
    % genes. Uses the original data, and the known selection size, to
    % determine an estimated weight threshold for differential expression.
    % The weight threshold is then used on permuted data to generate
    % null-hypothesis selections. If the real selection size is more
    % extreme than any of the null-hypothesis selection sizes, then the
    % data are considered differentially expressed.
    %
    % Parameters:
    %   input_filename  - Name of data file. Data file should be a CSV file
    %                       in 'tall' format, with genes in rows and 
    %                       samples in columns, comment lines should begin 
    %                       with '#', the first non-comment row should be a
    %                       classification row (0s and 1s for different
    %                       experimental conditions), and the first column
    %                       should hold gene names.
    %   classifications - Vector of sample classifications, either 0 or 1.
    %   n_permutations  - The number of permutations to generate and test.
    %   num_simulations - Number of simulations to generate. Simulations
    %                       are used to estimate the false discovery rate
    %                       of PED under different threshold constants.
    %                       Using more simulations may make FDR estimation
    %                       more accurate, but will increase runtime
    %                       linearly.
    %   fdr             - Average desired false discovery rate.
    %
    % Returns:
    %   model_sizes      - List of model sizes determined by selection on 
    %                       permutations.
    
    %% Create directory for permuted data
    dot_locs = strfind(input_filename, '.');
    last_dot = dot_locs(end);
    exp_location = input_filename(1:last_dot-1)
    permutation_location = [exp_location, '_permutations'];
    status = mkdir(permutation_location);
    if ~status
        disp(['mkdir failed with filename ', permutation_location])
        return;
    end
                            
    %% Get filename base
    filesep_locs = strfind(input_filename, filesep);
    if isempty(filesep_locs)
        last_filesep = 0;
    else
        last_filesep = filesep_locs(end);
    end
    filename_base = input_filename(last_filesep + 1:last_dot - 1);
    
    %% Create temporary files of permuted data
    [labels, data, classifications] = open_microarray_file(input_filename);  
    permutations = PED_permute_classifications(classifications, ...
                                               n_permutations);
    size(permutations)
    for permutation_num = 1:size(permutations, 2)
        disp(['permutation ', num2str(permutation_num)]);
        perm_classification = permutations(:, permutation_num);
        header = ['"# Permutation ', num2str(permutation_num), '"\n'];
        perm_name = ['permutation_', num2str(permutation_num), '.csv'];
        output_filename = fullfile(permutation_location, perm_name);
        write_microarray_from_data(data, output_filename, labels,...
                                   perm_classification, header);        
    end
    
    %% Find selection size for each permutation
    permutation_files    = dir(fullfile(permutation_location, '*.csv'));
    model_sizes          = -1 * ones(length(permutation_files), 1);
    
    for i = 1:length(permutation_files)
        model_sizes(i) = PED_estimate_threshold( ...
                            [permutation_location, filesep, ...
                             permutation_files(i).name], ...
                            num_simulations, fdr);
        disp(['Selection size of permutation #', num2str(i), ' = ', ...
              num2str(model_sizes(i))]);
    end       
end
