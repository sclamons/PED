function fd = fold_difference(a, b)
% Calculate the fold-difference between two lists a and b of numeric 
% values. 

    mean_a = mean(a);
    mean_b = mean(b);
    if mean_b > mean_a
        fd = mean_b / mean_a;
    else
        fd = -1.0 * mean_a / mean_b;
    end
end