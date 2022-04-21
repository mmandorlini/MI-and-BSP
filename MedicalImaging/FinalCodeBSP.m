%% BSP Topic 1 -- The potentialities of HRV as a tool for assessing fetal wellbeing.
%  The value of cardiotocographic signal for the diagnosis of pathological fetal condition
%  of Intrauterine Growth Restriction
%  _______________________________________________________________________________________
%  For all the functions detail see the file FunctionsBSP
clc;
clear;
close all;
ctg=[];

%% Import Data

% Dataset recordings are saved as arrays in cells from 1 to 10
ctg{1}=struct2array(load('Healthy_1.mat'));
ctg{2}=struct2array(load('Healthy_2.mat'));
ctg{3}=struct2array(load('Healthy_3.mat'));
ctg{4}=struct2array(load('Healthy_4.mat'));
ctg{5}=struct2array(load('Healthy_5.mat'));
ctg{6}=struct2array(load('IUGR_1.mat'));
ctg{7}=struct2array(load('IUGR_2.mat'));
ctg{8}=struct2array(load('IUGR_3.mat'));
ctg{9}=struct2array(load('IUGR_4.mat'));
ctg{10}=struct2array(load('IUGR_5.mat'));


%% Constants Definition

fs=2; %sampling frequency

VLF=[0 0.03]; %very low frequency range
LF=[0.03 0.15]; % low frequency range
MF=[0.15 0.5]; % mid frequency range
HF=[0.5 1]; % high frequency range

frequency_perc=[]; % array to store freqency percentage
frequency_area=[]; % array to store frequency areas
ctg_orig=ctg; % saving original ctg with no preprocessing

%% Preprocessing 

%     figure();
% Preprocessing of Healthy and IUGR data from cell 1 to cell 10 
for n=1:10 
    t=0:1/fs:length(ctg{n})/fs-1/fs; % time vector definition
    
    % Setting outliers as NaN
    outliers = find((ctg{n}<50)|(ctg{n}>200)); % outliers: abnormal values over 200 and under 50
    ctg{n}(outliers) = NaN;
    indici{n} = isnan(ctg{n}); % logical array of NaN values position 

    % Resampling - NaN removal
    ctg{n}=resample(ctg{n},t);

    % Moving average 
    ctg{n} = movmean(ctg{n},5);
    
    % CTG segmentation -- partition of the signal vector in nonoverlapping segments
    % [long segments] samples length=360 ; time length= 180s
    temp_sample_end_long=floor(length(ctg{n})/360)*360;
    segm_ctg_long{n} = buffer(ctg{n}(1:temp_sample_end_long),360,0);
    segm_indici_long{n} = buffer(indici{n}(1:temp_sample_end_long),360,0);
    
    % Removal of samples where NaN values are > 5% of the segment samples
    del_index_long = zeros*size(segm_indici_long{n},2);
    for i = 1:size(segm_indici_long{n},2)
        if sum(sum(segm_indici_long{n}(:,i)==1))>0.05*360
            del_index_long(i) = 1;
        end
    end
    segm_ctg_long{n}(:,del_index_long==1) = [];
    
    % CTG segmentation
    % [short segments] samples length=120 ; time length= 60s
    temp_sample_end_short=floor(length(ctg{n})/120)*120;
    segm_ctg_short{n}= buffer(ctg{n}(1:temp_sample_end_short),120,0);
    segm_indici_short{n} = buffer(indici{n}(1:temp_sample_end_short),120,0);
    
    % Removal of samples where NaN values are > 5% of the segment samples
    del_index_short = zeros*size(segm_indici_long{n},2);
    for i = 1:size(segm_indici_short{n},2)
        if sum(sum(segm_indici_short{n}(:,i)==1))>0.05*120
            del_index_short(i) = 1;
        end
    end
    segm_ctg_short{n}(:,del_index_short==1) = [];
end

%% Frequency Domain Analysis

for n=1:10
    [row_long, col_long]= size(segm_ctg_long{n});
    % row_long = number of elements for each long segment
    % col_long = total number of segments
    
    % arrays in which the spectral information will be saved for each segment of the different CTG recordings
    perc_VLF=[]; perc_LF=[]; perc_MF=[]; perc_HF=[]; 
    perc_LF_rid=[]; perc_MF_rid=[]; perc_HF_rid=[];
    LF_over_MFHF=[];
    area_VLF=[]; area_LF=[]; area_MF=[]; area_HF=[]; area_TOT=[]; area_TOT_rid=[]; 
    
%     figure();
    
    % PSD estimation for every segment
    for i=1:col_long
        segm_i=segm_ctg_long{n}(:,i)-mean(segm_ctg_long{n}(:,i)); % Mean removal
       
        % Parametric PSD estimation [YuleAR]
        nfft=1024;
        AR_order=16;
        [PSD_YAR, f_YAR] = pyulear(segm_ctg_long{n}(:,i),AR_order,nfft,fs);
        
        % Plotting together the estimations of the different segments:
        % subplot(ceil(col_long/3), 3, i), plot(f_YAR, PSD_YAR);
        % title(['Healthy-',num2str(n),' Segment-',num2str(i)]);

        % Frequency components
        freq= f_YAR;
        PSD_est= PSD_YAR;
        
        % Indexing frequency ranges VLF, LF, MF, HF
        f_psd_LF=and(ge(freq,LF(1)),le(freq,LF(2)));
        f_psd_HF=and(ge(freq,HF(1)),le(freq,HF(2)));
        f_psd_VLF=and(ge(freq,VLF(1)),le(freq,VLF(2)));
        f_psd_MF=and(ge(freq,MF(1)),le(freq,MF(2)));
        
        % Computing and storing frequency areas
        area_VLF_temp=trapz(PSD_est(f_psd_VLF));
        area_VLF=[area_VLF, area_VLF_temp];
        
        area_LF_temp=trapz(PSD_est(f_psd_LF));
        area_LF=[area_LF, area_LF_temp];
        
        area_MF_temp=trapz(PSD_est(f_psd_MF));
        area_MF=[area_MF, area_MF_temp];
        
        area_HF_temp= trapz(PSD_est(f_psd_HF));
        area_HF=[area_HF, area_HF_temp];
        
        area_TOT_temp=area_VLF_temp+area_LF_temp+area_MF_temp+area_HF_temp;
        area_TOT=[area_TOT, area_TOT_temp];
        area_TOT_rid_temp=area_LF_temp+area_MF_temp+area_HF_temp;
        area_TOT_rid=[area_TOT_rid, area_TOT_rid_temp];


        % Computing and storing frequency percentage
        perc_VLF=[perc_VLF, area_VLF_temp/area_TOT_temp*100];
        perc_LF=[perc_LF, area_LF_temp/area_TOT_temp*100];
        perc_MF=[perc_MF, area_MF_temp/area_TOT_temp*100];
        perc_HF=[perc_HF, area_HF_temp/area_TOT_temp*100];
        perc_LF_rid=[perc_LF_rid, area_LF_temp/area_TOT_rid_temp*100];
        perc_MF_rid=[perc_MF_rid, area_MF_temp/area_TOT_rid_temp*100];
        perc_HF_rid=[perc_HF_rid, area_HF_temp/area_TOT_rid_temp*100];
        LF_over_MFHF=[LF_over_MFHF, area_LF_temp/(area_MF_temp+area_HF_temp)];
    end

% Storing the collected paramters into structures
frequency_area{n} = struct('area_VLF',area_VLF,'area_LF',area_LF, 'area_MF',area_MF , 'area_HF',area_HF, 'area_TOT',area_TOT);
frequency_perc{n}= struct('VLF',perc_VLF,'LF',perc_LF, 'MF',perc_MF , 'HF',perc_HF, 'LFoverMFHF',LF_over_MFHF);
frequency_perc_rid{n}= struct('LF',perc_LF_rid, 'MF',perc_MF_rid , 'HF',perc_HF_rid, 'LFoverMFHF',LF_over_MFHF);
end

%% Time Domain Analysis

for n=1:10   
    [row_short, col_short]= size(segm_ctg_short{n});
    % row_short = number of elements for each long segment
    % col_short = total number of segments
    
    % Arrays in which the short term variability, delta index and interval index
    % will be saved for each segment of the different CTG recordings
    STV_temp=[]; Delta_temp=[]; II_temp=[];

    for j = 1:col_short    
        T24=[]; Tinf=[]; Tdiff=[]; 
      
        % T24 = downsampled signal in [ms]
        for i = 1:row_short
            T24(i) = 30000/segm_ctg_short{n}(i,j);
        end
        
        % Tinf = averaged T24 values over 5 samples windows
        k=0;
        for m = 1:24
            Tinf(m)= sum(T24((1+k):1:(5+k)))/5;
            k=k+5;
        end
        
        %Tdiff = absolute value of the difference between two subsequent
        %values of Tinf
        for l=1:23
            Tdiff(l) = abs(Tinf(l+1)-Tinf(l));
        end
        
        % Parameters computation
        STV_temp(j) = mean(Tdiff);
        Delta_temp(j) = max(Tinf)- min(Tinf);
        II_temp(j) = std(abs(Tdiff))/STV_temp(j);
    end
    
    % Storing the collected paramters
    STV{n} = mean(STV_temp);
    DELTA{n} = mean(Delta_temp);
    II{n} = mean(II_temp);
end

for n=1:10    
    [row_long, col_long]= size(segm_ctg_long{n}); 
    % row_long= number of elements for each long segment
    % col_long = total number of segments
    
    % Arrays in which mean, range, standard deviation, long term irregularities, variance,
    % median, will be saved for each segment of the different CTG recordings   
    mean_temp=[]; range_temp=[]; std_temp=[]; LTI_temp=[]; var_temp=[]; med_temp=[]; Samp_En_temp=[]; ApEn_temp=[]; 
    
    for j = 1:col_long
        T24=[]; Tinf=[]; m24=[];
        
        % T24= undersampled signal in [ms]
        for i = 1:row_long
            T24(i) = 30000/segm_ctg_long{n}(i,j);
        end
        

        % Tinf = averaged T24 values over 5 samples windows
        k=0;
        for m = 1:fix(row_long/5)-1
            Tinf(m) = sum(T24((1+k):1:(5+k)))/5;
            Tinf(m+1) = sum(T24((1+k+5):1:(5+k+5)))/5;
            k=k+5;
            m24(m) = sqrt((Tinf(m))^2+(Tinf(m+1))^2);
        end
        
        % Parameters computation
        mean_temp(j) = mean(segm_ctg_long{n}(:,j));
        range_temp(j) = max(segm_ctg_long{n}(:,j))-min(segm_ctg_long{n}(:,j));
        std_temp(j) = std(segm_ctg_long{n}(:,j));
        var_temp(j) = var(segm_ctg_long{n}(:,j));
        med_temp(j) = median(segm_ctg_long{n}(:,j));
        LTI_temp(j) = iqr(m24);
    end
    
    % Storing the collected paramters
    MEAN{n} = mean(mean_temp);
    RANGE{n} = mean(range_temp);
    STD{n} = mean(std_temp);
    VAR{n} = mean(var_temp);
    MEDIAN{n} = mean(med_temp);
    LTI{n} = mean(LTI_temp);
end

%% Complexity Analysis
scales = 15;

for n=1:10
    
    % Multiscale Sample Entropy at different scales
    for tau = 1:scales
        [e,A,B] = FunctionsBSP.multiscaleSampleEntropy(ctg{n},1,0.1,tau);
        sample_entropy_scale(tau,n) = e;
    end
    % Sample entropy with tau= 8 -- Time scale that allows better distinction between
    % Healty and IUGR recordings
    SampEn{n} = FunctionsBSP.multiscaleSampleEntropy(ctg{n},1,0.1,8); 
    ApEn(n) = FunctionsBSP.approx_entropy(1,0.1,ctg{n});
end

% Multiscale Sample Entropy analysis
figure,
for i=1:scales
    % [H] = Healthy recordings
    % [I] = IGUR recordings
    mean_entropy_H(i) = mean(sample_entropy_scale(i,1:5));
    mean_entropy_I(i) = mean(sample_entropy_scale(i,6:10));
    err_H(i) = std(sample_entropy_scale(i,1:5));
    err_I(i) = std(sample_entropy_scale(i,6:10));
end
errorbar(mean_entropy_H,err_H, 'g'), hold on, errorbar(mean_entropy_I,err_I, 'r'),hold off, xlim([0 16]), title('Multiscale Sample Entropy Analysis'), 
legend; % plotting the multiscale sample entropy trend with respect to the scale factor value

%% Parameters collection

for n = 1:10
    % Storing all the computed parameters in a structure
    Parameters_struct{n} = struct('MEAN',MEAN{n},'RANGE',RANGE{n}, 'STD',STD{n}, 'VAR',VAR{n}, 'MEDIAN',MEDIAN{n}, 'DELTA',DELTA{n}, 'II',II{n}, 'STV',STV{n}, 'LTI',LTI{n}, 'area_VLF',mean(frequency_area{1,n}.area_VLF), 'area_LF',mean(frequency_area{1,n}.area_LF), 'area_HF',mean(frequency_area{1,n}.area_HF), 'area_MF',mean(frequency_area{1,n}.area_MF), 'LF_perc_rid',mean(frequency_perc_rid{1,n}.LF), 'HF_perc_rid',mean(frequency_perc_rid{1,n}.HF), 'MF_perc_rid',mean(frequency_perc_rid{1,n}.MF), 'LFoverMFHF',mean(frequency_perc_rid{1,n}.LFoverMFHF), 'SAMP_EN',SampEn{n});
    
    % Storing all parameters in a matrix for further analysis
    % Parameters_matrix:
    % [observations] -> 10 (5 Healthy + 5 IUGR in this order)
    % [variables]
    Parameters_matrix(n,:) = [MEAN{n},
        RANGE{n},
        STD{n},
        VAR{n},
        MEDIAN{n},
        DELTA{n},
        II{n},
        STV{n},
        LTI{n},
        mean(frequency_area{1,n}.area_VLF),
        mean(frequency_area{1,n}.area_LF),
        mean(frequency_area{1,n}.area_HF),
        mean(frequency_area{1,n}.area_MF),
        mean(frequency_perc_rid{1,n}.LF),
        mean(frequency_perc_rid{1,n}.HF),
        mean(frequency_perc_rid{1,n}.MF),
        mean(frequency_perc_rid{1,n}.LFoverMFHF),
        SampEn{n}];
end

% Number of defined parameters
num_param = size(Parameters_matrix,2);

% Adding the target class value to each observation: 
% 1 for IUGR
% 0 for Healthy
Parameters_matrix(1:5,num_param+1) = 0;
Parameters_matrix(6:10,num_param+1) = 1;

%% Statistical Analysis

% Linear correlation coefficient between all pairs of columns (variables + target)
CORR = corrcoef(Parameters_matrix);
indici = abs(CORR)>=0.9;

% Barplot of linear correlaation coefficients with respect to the target
figure, bar(CORR(num_param+1,1:num_param)), title('Target linear correlation coefficients');

% Boxplots
figure,subplot(1,2,1), boxplot(Parameters_matrix(1:5,3)), ylim([1 13]), title('Healthy Standard Deviaton');
subplot(1,2,2), boxplot(Parameters_matrix(6:10,3)), ylim([1 13]), title('IUGR Standard Deviaton');

figure,subplot(1,2,1), boxplot(Parameters_matrix(1:5,8)), ylim([1 10]), title('Healthy STV');
subplot(1,2,2), boxplot(Parameters_matrix(6:10,8)), ylim([1 10]), title('IUGR STV');

figure,subplot(1,2,1), boxplot(Parameters_matrix(1:5,9)), ylim([1 40]), title('Healthy LTI');
subplot(1,2,2), boxplot(Parameters_matrix(6:10,9)), ylim([1 40]), title('IUGR LTI');

figure,subplot(1,2,1), boxplot(Parameters_matrix(1:5,17)), ylim([3.8 4.3]), title('Healthy LF/(MF+HF)');
subplot(1,2,2), boxplot(Parameters_matrix(6:10,17)), ylim([3.8 4.3]), title('IUGR LF/(MF+HF)');

figure,subplot(1,2,1), boxplot(Parameters_matrix(1:5,18)), ylim([1.2 2.2]), title('Healthy Sample Entropy');
subplot(1,2,2), boxplot(Parameters_matrix(6:10,18)), ylim([1.2 2.2]), title('IUGR Sample Entropy');

figure,subplot(1,2,1), boxplot(ApEn(1:5)), ylim([1 2.3]), title('Healthy Approximate Entropy');
subplot(1,2,2), boxplot(ApEn(6:10)), ylim([1 2.3]), title('IUGR Approximate Entropy');

%% Variable scaling and PCA

% Mean scaling of the variables values
scaled_param = (Parameters_matrix(:,1:(size(Parameters_matrix,2)-1)) - mean(Parameters_matrix(:,1:(size(Parameters_matrix,2)-1))));

% Standard/MinMax scaling of the variables values
for i=1:size(Parameters_matrix,2)-1
    std_scaled_param(:,i) = scaled_param(:,i)/std(Parameters_matrix(:,i));
    MinMax_scaled_param(:,i) = (Parameters_matrix(:,i)-min(Parameters_matrix(:,i)))/(max(Parameters_matrix(:,i))-min(Parameters_matrix(:,i)));
end

% Principal Component Analysis
[PCA_coeff,score,latent,tsquared,PCA_explained,mu] = pca(MinMax_scaled_param);

% Barplot of variance fraction explained by the different PCs
figure, bar(PCA_explained), title('PC explained variance fraction');

% Number of principal component used for the classification
num_pca=4;

% Barplot of coefficients of the different selected PCAs
figure, title('PC explained variance fraction'),
for i = 1:num_pca
    subplot(1,num_pca,i), bar(PCA_coeff(:,i));
end

% Total variance explained by the first num_pca components
PCA_explained_tot = sum(PCA_explained(1:num_pca));

% Create a new matrix for the variables
for i=1:num_pca
    for j=1:10
        PCA_MATRIX(j,i) = std_scaled_param(j,:)*PCA_coeff(:,i);
    end
end
target = Parameters_matrix(:,size(Parameters_matrix,2));

% Add the target class value to each observation
std_scaled_param = [std_scaled_param, target];
MinMax_scaled_param = [MinMax_scaled_param, target];
PCA_MATRIX = [PCA_MATRIX, target];

%% Classification model selection

% Training of the selected best model
% Validation of the model with the leave one out approach
[trainedClassifier, validationAccuracy] = FunctionsBSP.trainClassifier(MinMax_scaled_param);

