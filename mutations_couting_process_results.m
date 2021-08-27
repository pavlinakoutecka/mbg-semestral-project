
% Script for processing results of the functions in mutations_counting.py

clear all; clc; close all;


%% Numbers of mutations and diff. sequences
num_sequences = '10000';
M = readtable(sprintf('%ss_mutations_counts.csv', num_sequences));
M = table2array(M);
D = readtable(sprintf('%ss_diff_seq_counts.csv', num_sequences));
D = table2array(D);

% filtration - thresholding using numbers of different sequences
% threshold set as 12 % of all sequences
% threshold 11 % does not lead to significant filtration
th = 0.15 * str2double(num_sequences); 
M_filt = M(:, D(2,:) > th);

%%
figure;
scatter(M_filt(1,:), M_filt(2,:), 10, 'filled');
title('Numbers of significant mutations on positions');
xlabel('Position');
ylabel('Number of mutations');

figure;
bar(M(1,:), M(2,:));
title('Numbers of mutations on all positions');
xlabel('Position');
ylabel('Number of mutations');

figure;
bar(D(1,:), D(2,:));
title('Numbers of different sequences on all positions');
xlabel('Position');
ylabel('Number of diff. sequences');

%% Occurrences of nucleotides
num_sequences = '10000';
N = readtable(sprintf('%ss_nucleotides_occurrences.csv', num_sequences));
N = table2array(N);






