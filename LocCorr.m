function [Fin_tab, Fin_lab, Fin_corr] = LocCorr(FileName, LocThresh)
%LOCCORR is a Function for the analysis of squencing data correlations
%   --------------------INPUTS---------------------------------------------
%   FileName - insert name of file from the directory, eventually the path,
%   it the file is located elsewhere. File type: FASTA.
%   LocThresh - threshold for number of sequences we wnat to have a
%   mutation before we deem it significant. In decimal numbers. Recommended
%   0.005 (any lower will most likely fill memory).
%   ---------------------OUTPUTS-------------------------------------------
%   Fin_tab - cell array of tables showing the counts of mutated sequences
%   for the given loci and mutation types
%   Fin_lab - cell array of tables showing the specific mutation in integer
%   form based on MATLAB's nt2int function (1 - A, 2 - C, 3 - G, 4 -T)    
%   Fin_corr - matrix showing the loci, where mutations have occured
    %% Load data, change nucleotides to integers
    [name, sequence] = fastaread(FileName);
    int = [];
    for i = 1:size(sequence, 2)
        s_temp = char(sequence(1, i)); 
        int(i, :) = nt2int(s_temp);
    end

    %% Define reference as mode at all loci 
    ref = mode(int, 1);

    %% Define mutated loci as 1 and original as 0
    muts = [];
    for i = 1:size(int, 1)
        for j = 1:size(int, 2)
            if int(i, j) == ref(1, j)
                muts(i, j) = 0;
            else
                muts(i, j) = 1;
            end
        end
    end

    %% Ignore all loci, where the number of loactions is smaller than 0.5 %
    thresh = find(sum(muts) < LocThresh*size(sequence, 2));
    labels = 1:size(int, 2);
    new = muts; 
    new(:, thresh) = []; % remove from the matrix of mutations
    new_int = int;
    new_int(:, thresh) = []; % remove from the matrix of nucleotide integers
    mut_nkt = new_int.*new; % remove those integers, which aren't mutated
    labels(thresh) = []; % remove these from the labels of loci as well

    %% Get correlation of loci (mutation or not)
    [correlation, pvals] = corrcoef(new);
    alpha1 = 0.05/((size(new, 2))^2);
    [x, y] = find(pvals < alpha1); % Use only those, where the p value is smaller than 0.05 with Bonferroni correction - (0.05/(598^2)));

    %% Use Chi-square test to find correlation of nucleotides at mutated loci
    chi = [];
    p = [];
    for i = 1:length(x)
        x_temp = mut_nkt(:, x(i));
        y_temp = mut_nkt(:, y(i));
        % REMOVE NONMUTATED LOCI
        x_temp((x_temp == 0)) = NaN;
        x_temp((y_temp == 0)) = NaN;
        y_temp(isnan(x_temp)) = NaN;
        [~, chi(i), p(i), ~] = crosstab(x_temp, y_temp);
    end

    %% 
    k = find(p < (alpha1/size(new, 2))); % Find p values smaller than 0.05 with Bonferroni correction
    xi = x(k);
    yi = y(k);
    correlated_loci = [labels(xi); labels(yi)];


    X_tables = {};
    X_lables = {};
    for j = 1:length(k)
        k_temp = k(j);
        x_temp = mut_nkt(:, x(k_temp));
        y_temp = mut_nkt(:, y(k_temp));
        % REMOVE NONMUTATED LOCI
        x_temp(find(x_temp == 0)) = NaN;
        x_temp(find(y_temp == 0)) = NaN;
        % Remove unconcrete mutations 
        x_temp(x_temp > 4) = NaN;
        x_temp(y_temp > 4) = NaN;
        y_temp(isnan(x_temp)) = NaN;
        [X_tables{j}, ~, ~, X_lables{j}] = crosstab(x_temp, y_temp);
    end

    %% Take into account only combinations with one mutatuion type for each locus
    Fin_tab = {};
    Fin_lab = {};
    Fin_corr = [];
    for i = 1:size(X_tables, 2)
        if size(X_tables{i}) == [1, 2]
            Fin_tab{end + 1} = X_tables{i};
            Fin_lab{end + 1} = X_lables{i};
            Fin_corr = [Fin_corr, correlated_loci(:, i)];
        end
    end

end

