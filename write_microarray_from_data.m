function write_microarray_from_data(data, output_filename, labels,...
                                    classifications, header)
% write_microarray_from_data
% Last edited March 12, 2014
% Samuel Clamons
% 
% Writes data for a microarray experiment to CSV file in "tall" format with
% the following specifications:
%   * Basic format is '\n' separated, comma-delineated flat text.
%   * Features (genes, miRNAs, etc) are in rows
%   * Samples (microarrays) are in columns
%   * Header lines begin with a '#' character, or the two-character
%       combination '"#' if the line is quoted to prevent mid-line breaks
%       when displayed in spreadsheet software
%   * The first non-comment row should be a classification row (0s and 1s 
%       for differentexperimental conditions)
%   * The first column should hold feature (gene) names
%
% Parameters:
%   data                    - matrix containing the microarray data, in
%                               tall format (genes in rows, samples in
%                               columns).
%   output_filename         - the name of the CSV output file. 
%   labels                  - cell vector of string names of genes.
%   classifications         - vector of classifications (0 or 1) of columns
%                               in the experiment.
%   header                  - string to add to the beginning of the 
%                               microarray as a commented header. Used for 
%                               adding descriptive information to the 
%                               array.
%
% Output: 
%   No return value 
%
% Side effects: 
%    Prints a microarray to <filename> in 'tall' format (genes in rows,
%    replicates in columns).
%
% Use example: 
%   labels = cell(1,2);
%   labels{1} = 'gene1';
%   labels{2} = 'gene2';
%   data = [.1, .8, 2.1, .2; .3, 1.2, 1.7, .15];
%   write_microarray_from_data(data, './example_output.csv', labels,...
%                               [0, 0, 1, 1]', ...
%                               '"# This text appears as a header"')
%                      
    %% Sanity check input
    [num_genes, num_samples] = size(data);
    if num_genes ~= length(labels)
        error(strcat('Warning: Data matrix and gene labels imply ',...
                     'different numbers of genes.'))
    end
            
    sprintf('Starting to write microarray data...\n')
   
    %% Open output file for writing
    [microarrayFile, msg] = fopen(output_filename, 'w');
    if microarrayFile < 0
        disp(pwd())
        disp(msg)
    end

    %% Print header
    sprintf('Printing header\n')
    fprintf(microarrayFile, header);

    %% Print classification row
    sprintf('Printing classification row\n')
    fprintf(microarrayFile, 'Gene Name');
    for i = 1:num_samples
        if classifications(i) == 0
            fprintf(microarrayFile, ',0');
        else
            fprintf(microarrayFile, ',1');
        end
    end
   
    %% Print data
    sprintf('Printing genes...\n')
    if ismember(NaN, data)
        disp(strcat('NaN detected in file ', microarray_filename))
    end
    for i = 1:num_genes
        fprintf(microarrayFile, '\n%s', labels{i});
        for j = 1:num_samples
            fprintf(microarrayFile, ',%f', data(i,j));
        end
    end
    fclose(microarrayFile);
    sprintf('Finished writing to %s', output_filename)                    
end
                       