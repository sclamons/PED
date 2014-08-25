function [labels, M, y] = open_microarray_file(filename)
% open_microarray_file
% Last edited March 31, 2014
% Samuel Clamons
%
% Read a CSV file of expression data in "tall" format written by the
% function "write_microarray_from_data". 
%
% Parameters:
%   filename - Name of the CSV file with experimental expression data.
%
% Return:
%   labels  - Cell vector of the names of features (genes) as discovered in
%               the first column of the data file.
%   M       - Expression data matrix, in the same format (tall or fat) as
%               the data file.
%   y       - Column vector of experimental conditions, either 0s or 1s. 
%
    file_data = importdata(filename);
        
    while(strfind(file_data.textdata{1}, '#') == 1 | ...
          strfind(file_data.textdata{1}, '#') == 2)
        file_data.textdata(1) = [];
    end
    
    y = file_data.data(1, :)';
    M = file_data.data(2:end, :);
    if size(M,1) > 0
        labels = file_data.textdata(2:1 + size(M,1));        
    else
        labels = {};
    end
end