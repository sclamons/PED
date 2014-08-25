function [selection_size, is_differential_expression] = ...
                                 PED_select_genes(input_filename,...
                                                  output_filename,...
                                                  num_simulations, fdr)
    % PED_select_genes
    % Last edited May 1, 2014
    % Samuel Clamons
    %
    % Determine which genes are differentially expressed in a data set by
    % PED. Automatically determines a threshold constant, then applies that
    % threshold constant to select genes. Also runs a differential
    % expression test with that constant, and outputs the result of the
    % test.
    %
    % Requires a CSV file of the experiment data in "tall" format, with
    % features (genes) in rows and samples in columns. To make an
    % appropriate CSV from a DataMatrix object, use PED_DataMatrix_to_CSV.
    %
    % Parameters:
    %   input_filename  - Name of data file. Data file should be a CSV file
    %                       in 'tall' format, with genes in rows and 
    %                       samples in columns, comment lines should begin 
    %                       with '#', the first non-comment row should be a
    %                       classification row (0s and 1s for different
    %                       experimental conditions), and the first column
    %                       should hold gene names.
    %   output_filename - Name of output file, which will hold selected
    %                       genes after execution.
    %   num_simulations - Number of simulations to generate. Simulations
    %                       are used to estimate the false discovery rate
    %                       of PED under different threshold constants.
    %                       Using more simulations may make FDR estimation
    %                       more accurate, but will linearly increase 
    %                       runtime.
    %   fdr             - Maximum allowed false discovery rate.
    %
    % Returns: 
    %   selection_size              - Number of genes selected as 
    %                                   differentially expressed.
    %   is_differential_expression  - 1 if differential expression is 
    %                                   detected, 0 otherwise.
    %
    % Outputs: Writes a list of differentially-expressed genes to
    % output_file.
    
    %% Find optimal selection size
    threshold = PED_estimate_threshold(input_filename, num_simulations, ...
                                       fdr);
                         
    %% Optimize by PED regression and select genes
    [labels, data, classifications] = open_microarray_file(input_filename);
    data = data';
    weights = PED_select_genes_from_threshold(data, classifications, ...
                                              threshold, 'size');
    selections = find(weights);
    
    %% Test for differential expression
    n_permutations = 9;
    %n_permutations = 1;
    model_sizes = PED_estimate_permutation_thresholds(input_filename,...
                                                      n_permutations, ...
                                                      num_simulations, ...
                                                      fdr);
    is_differential_expression = length(selections) > max(model_sizes);
                                        
    %% Open output file for writing
    [outfile, msg] = fopen(output_filename, 'w');
    if outfile < 0
        disp(pwd())
        disp(['Opening file ', output_filename, ' for writing'])
        disp(msg)
        return;
    end
    
    %% Print output
    fprintf(outfile, ['"# Genes selected from file:"\n"#', ...
                      input_filename, '"\n"# using FDR = ', ...
                      num2str(fdr), ' and ', num2str(num_simulations), ...
                      ' simulations"\n"# Differential expression ']);
    if ~is_differential_expression
        fprintf(outfile, ' NOT ');
    end
    fprintf(outfile, ['detected."']);
    fprintf(outfile, ['\n"# Real data model size: ', num2str(length(weights)),...
                      '"\n"# Permutation model sizes:"']);
    for i = 1:length(model_sizes)
        fprintf(outfile, ['\n"#\t', num2str(model_sizes(i)), '"']);
    end
    fprintf(outfile, '\n"Gene","Weight","Regulation Direction"');
    
    for i = 1:length(selections)
        gene_idx = selections(i);
        fprintf(outfile, ['\n', labels{gene_idx}, ',', ...
                          num2str(weights(gene_idx)), ',',...
                         ]);
        if mean(data(find(classifications), gene_idx)) < ...
           mean(data(find(classifications - 1), gene_idx))
            fprintf(outfile, 'downregulated');
        else
            fprintf(outfile, 'upregulated');
        end
    end
        
    %% Set output variable
    selection_size = nnz(weights);
end