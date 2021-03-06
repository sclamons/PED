function simulation_dir = PED_generate_simulations(num_microarrays, ...
                                                   min_fold_diff, ...
                                                   microarray_filename, ...
                                                   condition_number)
% PED_generate_simulations                                            
% Last edited March 12, 2014
% Samuel Clamons
% 
% Creates and writes microarray using another microarray data file (usually
% a real one) as the basis. Condition 0 of the new microarray is identicial
% to the condition taken from the other microarray file. Condition 1 is
% generated in a gaussian fashion with the mean and standard deviation
% estimated from condition 0.
%
% A differential expression distribution is estimated from the
% fold-difference distribution of the given data, ignoring any
% fold-difference of absolute value below a given threshold min_fd. That
% distribution is used to inject fold-differences into the simulated data.
% 
%
% Parameters:
%   num_microarrays     - The number of microarray files to be written for 
%                           this experiment.
%   min_fold_diff       - Smallest fold-difference that will be counted as
%                           differentially expressed and added to the data.
%   microarray_filename - Filename of the microarray to be read.
%   condition_number    - Classification (usually 0 or 1) to be copied from
%                           the file specified by microarray_filename.
%   
%
% Returns: 
%   The location of the microarrays written by this function
%
% Side effects: 
%   Writes <num_microarrays> synthetic microarray files to a new experiment
%   folder under "./simulations/", as well as a text file of parameters
%   used.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% VERIFY DIRECTORY STRUCTURE %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Simulation folder name is based on the name of the data file from
    % which the simulations are generated. The workspace is the location of
    % that data file.
    data_file_filesep_locs = strfind(microarray_filename, filesep);
    if isempty(data_file_filesep_locs)
        workspace                = strcat('.', filesep);
        microarray_filename_only = microarray_filename;
    else
        workspace = microarray_filename(1:data_file_filesep_locs(end));
        microarray_filename_only = ...
                    microarray_filename(data_file_filesep_locs(end)+1:end);
    end
    
    data_file_dot_locs = strfind(microarray_filename_only, '.');
    data_file_no_dot = ...
                     microarray_filename_only(1:data_file_dot_locs(end)-1);
    
    simulation_dir = strcat(workspace, data_file_no_dot,'_simulations', ...
                            filesep);
    
    if exist(simulation_dir, 'dir') == 0
        mkdir(simulation_dir);
    else        
        disp(['Simulation directory ', simulation_dir, ' already exists']);
        disp('Skipping simulation generation.');
        return;
    end
        
    target_dir_data = fullfile(simulation_dir, 'data');
    target_dir_results = fullfile(simulation_dir, 'results');
    if exist(target_dir_data, 'dir') == 0
        mkdir(target_dir_data);
    end
    if exist(target_dir_results, 'dir') == 0
        mkdir(target_dir_results);
    end

    exp_location = simulation_dir;
    
    %%%%%%%%%%%
    %% SETUP %%
    %%%%%%%%%%%
    
    % Write parameter report
    param_report_string = [...
         '# Parameters for ', ...
            num2str(strrep(exp_location, '\', '/')), ':',...
         '\n#\tnum_microarrays = ', ...
            num2str(num_microarrays),...
         '\n#\tmin_threshold = ', ...
            num2str(min_fold_diff),...
         '\n#\tmicroarray_filename = ', ...
            microarray_filename,...
         '\n#\tcondition_number = ', ...
            num2str(condition_number),...
         '\n'];
    
    [param_file, msg] = ...
               fopen(strcat(exp_location, filesep, 'parameters.txt'), 'w');
    if param_file < 0
        disp(pwd())
        disp(msg)
    end
    fprintf(param_file, param_report_string);
    
    % Generate microarray header
    header = ['"# Synthetic microarray data generated by the matlab ',...
              ' script" \n# "generate_microarray_data_mimic.m".',...
              '\n"# Gene distributions and fold-difference',...
              ' distributions are based on data from file "\n#',...
              microarray_filename, '\n"# from condition ', ...
              num2str(condition_number), ', with a minimum ',...
              'fold-difference of ', num2str(min_fold_diff), '"\n'];
            
    
    % Read out microarray
    [~, target_data, classifications] = ...
                                open_microarray_file(microarray_filename);
    [num_genes, ~] = size(target_data);
    
    % Use data from the given condition of the sample data for condition 0
    target_condition_indices = find(classifications == condition_number);
    samples_per_condition    = length(target_condition_indices);
    num_samples              = 2 * samples_per_condition;    
    cond_0_indices           = 1:samples_per_condition;
    cond_1_indices           = samples_per_condition+1 : num_samples;
    data                     = zeros(num_genes, num_samples);
    data(:, cond_0_indices)  = target_data(:, target_condition_indices);
    
    % Set up new classification vector
    new_classifications = zeros(1, num_samples);
    new_classifications(cond_1_indices) = 1;
    
    %%%%%%%%%%%
    %% PRINT %%
    %%%%%%%%%%%
        
    % For each microarray to print....
    for microarray_num = 1:num_microarrays
        sprintf('Writing!')
        % Generate data for condition 1 based on condition 0
        for i = 1:num_genes            
            gene_mean = mean(data(i, cond_0_indices));
            gene_std  =  std(data(i, cond_0_indices));
            data(i, cond_1_indices) = ...
                    normrnd(gene_mean, gene_std, 1, samples_per_condition);
        end

        % Set up label array
        labels = cell(1, num_genes);

        % Add differential expression to data.
        no_diff_label_counter = 1;
        for i = 1:num_genes
            gene_fd = fold_difference(target_data(i, cond_0_indices),...
                                      target_data(i, cond_1_indices));
            if abs(gene_fd) >= min_fold_diff
                if gene_fd > 0
                    data(i, cond_1_indices) = ...
                                       data(i, cond_1_indices) .* gene_fd;
                    labels{i} = strcat('upreg_by_', num2str(gene_fd));
                else
                    data(i, cond_1_indices) = ...
                                       data(i, cond_1_indices) ./ -gene_fd;
                    labels{i} = strcat('downreg_by_', num2str(-gene_fd));
                end            
            else
                labels{i} = ...
                        strcat('no_diff_', num2str(no_diff_label_counter));
                no_diff_label_counter = no_diff_label_counter + 1;
            end
        end

        sim_name = ['simulation_', num2str(microarray_num), '.csv'];
        output_filename = fullfile(exp_location, 'data', sim_name);
        
        % Actually print the data!
        write_microarray_from_data(data, output_filename, labels, ...
                                              new_classifications, header);
        sprintf('microarray %d printed\n', microarray_num)
    end
end
