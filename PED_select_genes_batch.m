function PED_select_genes_batch(input_directory, output_directory, fdr)
    % PED_select_genes_batch
    % Last edited April 1, 2014
    % Samuel Clamons
    %
    % Run PED on every data file in a given directory. Prints a report of
    % the weights discovered for each file, plus an overview of the
    % accuracies and selection sizes of all of the files. 
    %
    %   input_directory     - String of the name of the folder containing 
    %                           data.
    %   output_directory    - String of the name of the folder in which to 
    %                           write individual reports summary report.
    %   fdr                 - Maximum allowed false discovery rate.
    %
    % Usage example: 
    %   To run PED with a maximum FDR of 5% on every CSV file in 
    %   'data/experiment_1_data', print the results to 
    %   'results/experiment_1_results:
    %       in_directory = 'data/experiment_1_data/';
    %       out_directory = 'results/experiment_1_results/';
    %       PED_select_genes_batch(in_directory, out_directory, 0.05);
    %   
    
    
    %% Set up variables, locate files    
    data_files = dir(fullfile(input_directory, '*.csv'));
    filenames     = cell(1,length(data_files));
    sizes         = -1 * ones(length(data_files), 1);
    diff_detected = -1 * ones(length(data_files), 1);
    expression_scores = -1 * ones(length(data_files), 1);
    
    %% Create an output directory if necessary
    if exist(output_directory) ~= 7
        warning(['Output directory ', output_directory, ...
                 ' does not exist or is not a directory. ',...
                 'Creating directory.']);
        [s, mess, messid] = mkdir(output_directory);
    end
    
    %% Run PED on each file
    for i = 1:length(data_files)
        filenames{i}  = data_files(i).name;
        in_filename   = fullfile(input_directory, filenames{i});
        dot_locations = strfind(filenames{i}, '.');
        last_dot      = dot_locations(end);

        out_filename  = strcat(output_directory, filesep, ...
                        filenames{i}(1:last_dot-1), '_results.csv');
        [sizes(i), diff_detected(i)] = ...
                     PED_select_genes(in_filename, out_filename, 10, fdr);     
    end
    
    %% Open report file for writing
    [file_handle, msg] = fopen(strcat(output_directory, filesep, ...
                                'PED_batch_results.csv'), 'w');
    if file_handle < 0
        disp(pwd())
        disp(msg)
        return;
    end
    
    %% Write summary
    fprintf(file_handle, ['# With FDR ', num2str(fdr), ...
                          '\n# Dataset,Selection Size,',...
                          'Differential Expression Detected?,']);
    for i = 1:length(data_files)
        fprintf(file_handle, '\n%s,%d,%d,%f', filenames{i}, ...
                sizes(i), diff_detected(i));
    end
    fprintf(file_handle,'\nAverage,%f',mean(sizes));
end