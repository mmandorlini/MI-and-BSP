classdef FunctionsBSP
% Collection of functions necessary for the code :
properties
end
methods (Static)
    
function [ e, A, B ] = multiscaleSampleEntropy( x, m, r, tau )
% MULTISCALE SAMPLE ENTROPY
%
% Based on "Multiscale entropy analysis of biological signals"
% By Madalena Costa, Ary L. Goldberger, and C.-K. Peng
% Published on 18 February 2005 in Phys. Rev. E 71, 021906.
%
% This code was implemented by John Malik on 26 April 2017.
% Contact: john.malik@duke.edu
switch nargin
    case 1
        m = 2;
        r = 0.15;
        tau = 1;
    case 2
        r = 0.15;
        tau = 1;
    case 3
        tau = 1;
end
% coarse signal
y = mean(buffer(x(:), tau), 1);
% (m+1)-element sequences
X = buffer(y, m + 1, m, 'nodelay')';
% matching (m+1)-element sequences
A = sum(pdist(X, 'chebychev') < r * nanstd(x, 1));
% matching m-element sequences
X = X(:, 1:m);
B = sum(pdist(X, 'chebychev') < r * nanstd(x, 1));
% take log
if A == 0 || B == 0
    e = NaN;
    return
end
e = log(B / A);
end

function [apen] = approx_entropy(n,r,a)
%% Code for computing approximate entropy for a time series: Approximate
% Entropy is a measure of complexity. It quantifies the unpredictability of
% fluctuations in a time series
% To run this function- type: approx_entropy('window length','similarity measure','data set')
% i.e  approx_entropy(5,0.5,a)
% window length= length of the window, which should be considered in each iteration
% similarity measure = measure of distance between the elements
% data set = data vector
% small values of apen (approx entropy) means data is predictable, whereas
% higher values mean that data is unpredictable
% concept boorowed from http://www.physionet.org/physiotools/ApEn/
% Author: Avinash Parnandi, parnandi@usc.edu, http://robotics.usc.edu/~parnandi/
%%
data =a;
for m=n:n+1; % run it twice, with window size differing by 1
set = 0;
count = 0;
counter = 0;
window_correlation = zeros(1,(length(data)-m+1));
for i=1:(length(data))-m+1,
    current_window = data(i:i+m-1); % current window stores the sequence to be compared with other sequences
    
    for j=1:length(data)-m+1,
    sliding_window = data(j:j+m-1); % get a window for comparision with the current_window
    
    % compare two windows, element by element
    % can also use some kind of norm measure; that will perform better
    for k=1:m,
        if((abs(current_window(k)-sliding_window(k))>r) && set == 0)
            set = 1; % i.e. the difference between the two sequence is greater than the given value
        end
    end
    if(set==0) 
         count = count+1; % this measures how many sliding_windows are similar to the current_window
    end
    set = 0; % reseting 'set'
    
    end
   counter(i)=count/(length(data)-m+1); % we need the number of similar windows for every cuurent_window
   count=0;
i;
end  %  for i=1:(length(data))-m+1, ends here
counter;  % this tells how many similar windows are present for each window of length m
%total_similar_windows = sum(counter);
%window_correlation = counter/(length(data)-m+1);
correlation(m-n+1) = ((sum(counter))/(length(data)-m+1));
 end % for m=n:n+1; % run it twice   
   correlation(1);
   correlation(2);
apen = log(correlation(1)/correlation(2));
end   

% trainClassifier with leave one out approach
function [trainedClassifier, validationAccuracy] = trainClassifier(trainingData)
% [trainedClassifier, validationAccuracy] = trainClassifier(trainingData)
% returns a trained classifier and its accuracy. This code recreates the
% classification model trained in Classification Learner app. Use the
% generated code to automate training the same model with new data, or to
% learn how to programmatically train models.
%
%  Input:
%      trainingData: a matrix with the same number of columns and data type
%       as imported into the app.
%
%  Output:
%      trainedClassifier: a struct containing the trained classifier. The
%       struct contains various fields with information about the trained
%       classifier.
%
%      trainedClassifier.predictFcn: a function to make predictions on new
%       data.
%
%      validationAccuracy: a double containing the accuracy in percent. In
%       the app, the History list displays this overall accuracy score for
%       each model.
%
% Use the code to train the model with new data. To retrain your
% classifier, call the function from the command line with your original
% data or new data as the input argument trainingData.
%
% For example, to retrain a classifier trained with the original data set
% T, enter:
%   [trainedClassifier, validationAccuracy] = trainClassifier(T)
%
% To make predictions with the returned 'trainedClassifier' on new data T2,
% use
%   yfit = trainedClassifier.predictFcn(T2)
%
% T2 must be a matrix containing only the predictor columns used for
% training. For details, enter:
%   trainedClassifier.HowToPredict

% Auto-generated by MATLAB on 14-May-2020 23:06:59


% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18'};
predictors = inputTable(:, predictorNames);
response = inputTable.column_19;
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];

% Train a classifier
% This code specifies all the classifier options and trains the classifier.
subspaceDimension = max(1, min(9, width(predictors) - 1));
classificationEnsemble = fitcensemble(...
    predictors, ...
    response, ...
    'Method', 'Subspace', ...
    'NumLearningCycles', 30, ...
    'Learners', 'knn', ...
    'NPredToSample', subspaceDimension, ...
    'ClassNames', [0; 1]);

% Create the result struct with predict function
predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
ensemblePredictFcn = @(x) predict(classificationEnsemble, x);
trainedClassifier.predictFcn = @(x) ensemblePredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedClassifier.ClassificationEnsemble = classificationEnsemble;
trainedClassifier.About = 'This struct is a trained model exported from Classification Learner R2019b.';
trainedClassifier.HowToPredict = sprintf('To make predictions on a new predictor column matrix, X, use: \n  yfit = c.predictFcn(X) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nX must contain exactly 18 columns because this model was trained using 18 predictors. \nX must contain only predictor columns in exactly the same order and format as your training \ndata. Do not include the response column or any columns you did not import into the app. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18'};
predictors = inputTable(:, predictorNames);
response = inputTable.column_19;
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];

% Perform cross-validation
partitionedModel = crossval(trainedClassifier.ClassificationEnsemble, 'KFold', 10);

% Compute validation predictions
[validationPredictions, validationScores] = kfoldPredict(partitionedModel);

% Compute validation accuracy
validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');

end

end
end