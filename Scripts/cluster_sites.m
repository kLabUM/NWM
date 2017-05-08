clear all, close all
data = readtable('cluster_data.txt');
%% standardize data
%leave covs, granger, dynamics the same, but setting NaN's to be bad
standard = data;
ind = isnan(standard.covs_train);
standard.covs_train(ind) = -1;

standard.lags_train = (standard.lags_train - mean(standard.lags_train))/std(standard.lags_train);

ind = isnan(standard.covs_test);
standard.covs_test(ind) = -1;


standard.lags_test = (standard.lags_test - mean(standard.lags_test))/std(standard.lags_test);
standard.bottom_width = (standard.bottom_width - mean(standard.bottom_width))/std(standard.bottom_width);
standard.elevation = (standard.elevation - mean(standard.elevation))/std(standard.elevation);
standard.mannings = (standard.mannings - mean(standard.mannings))/std(standard.mannings);
standard.slope = (standard.slope - mean(standard.slope))/std(standard.slope);

standard.tf21 = standard.tf21/100;
ind = find(standard.tf21<-1);
standard.tf21(ind) = -1;
ind = isnan(standard.tf21);
standard.tf21(ind) = -1;

standard.tf31 = standard.tf31/100;
ind = find(standard.tf31<-1);
standard.tf31(ind) = -1;
ind = isnan(standard.tf31);
standard.tf31(ind) = -1;

standard.tf32 = standard.tf32/100;
ind = find(standard.tf32<-1);
standard.tf32(ind) = -1;
ind = isnan(standard.tf32);
standard.tf32(ind) = -1;

standard.ss1 = standard.ss1/100;
ind = find(standard.ss1<-1);
standard.ss1(ind) = -1;
ind = isnan(standard.ss1);
standard.ss1(ind) = -1;

standard.ss2 = standard.ss2/100;
ind = find(standard.ss2<-1);
standard.ss2(ind) = -1;
ind = isnan(standard.ss2);
standard.ss2(ind) = -1;

standard.ss3 = standard.ss3/100;
ind = find(standard.ss3<-1);
standard.ss3(ind) = -1;
ind = isnan(standard.ss3);
standard.ss3(ind) = -1;

standard.ss4 = standard.ss4/100;
ind = find(standard.ss4<-1);
standard.ss4(ind) = -1;
ind = isnan(standard.ss4);
standard.ss4(ind) = -1;

standard.oe221 = standard.oe221/100;
ind = find(standard.oe221<-1);
standard.oe221(ind) = -1;
ind = isnan(standard.oe221);
standard.oe221(ind) = -1;

standard.oe211 = standard.oe211/100;
ind = find(standard.oe211<-1);
standard.oe211(ind) = -1;
ind = isnan(standard.oe211);
standard.oe211(ind) = -1;

standard.usgs_dist = (standard.usgs_dist - mean(standard.usgs_dist))/std(standard.usgs_dist);
%% assign "good" label. 
standard.good_tf21 = zeros(length(data.ifis_id),1);
ind = find(standard.tf21>=.4);
standard.good_tf21(ind) = 1;

standard.good_tf31 = zeros(length(data.ifis_id),1);
ind = find(standard.tf31>=.4);
standard.good_tf31(ind) = 1;

standard.good_tf32 = zeros(length(data.ifis_id),1);
ind = find(standard.tf32>=.4);
standard.good_tf32(ind) = 1;

standard.good_ss1 = zeros(length(data.ifis_id),1);
ind = find(standard.ss1>=.4);
standard.good_ss1(ind) = 1;

standard.good_ss2 = zeros(length(data.ifis_id),1);
ind = find(standard.ss2>=.4);
standard.good_ss2(ind) = 1;

standard.good_ss3 = zeros(length(data.ifis_id),1);
ind = find(standard.ss3>=.4);
standard.good_ss3(ind) = 1;

standard.good_ss4 = zeros(length(data.ifis_id),1);
ind = find(standard.ss4>=.4);
standard.good_ss4(ind) = 1;

standard.good_oe221 = zeros(length(data.ifis_id),1);
ind = find(standard.oe221>=.4);
standard.good_oe221(ind) = 1;

standard.good_oe211 = zeros(length(data.ifis_id),1);
ind = find(standard.oe211>=.4);
standard.good_oe211(ind) = 1;

classification = table(standard.ifis_id,standard.covs_train,standard.lags_train,standard.covs_test,standard.lags_test,standard.granger,standard.bottom_width,standard.elevation,standard.mannings,standard.dynamics,standard.lat,standard.lon,standard.usgs_dist,standard.slope,standard.good_tf21,standard.good_tf31,standard.good_tf32,standard.good_ss1,standard.good_ss2,standard.good_ss3,standard.good_ss4,standard.good_oe221,standard.good_oe211,'VariableNames',standard(:,[1:10,20:22,25:34]).Properties.VariableNames);
writetable(classification,'/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/classification_inputs/all.txt','Delimiter','\t');

%% plot

newdata = table2array(data(:,[2:19,22]));

states = shaperead('usastatelo' , 'UseGeoCoords' , true);
flows = shaperead('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/NWM_parameters/NHDFlowline.shp');
k = 5;
map = colormap(hsv(k));
[test,centers] = kmeans(newdata,k);
close all
hold on
plot(states(15).Lon,states(15).Lat)
for i = 1:k
    ind = find(test==i);
    scatter(data.lon(ind),data.lat(ind),50,map(i,:),'filled')
end
scatter(data.usgs_lon,data.usgs_lat,'x')
leg = ['Iowa',string(1:k),'USGS'];
legend(leg)
hold off

%% train models and cluster sites
clear all,close all
data = readtable('classification_inputs/all.txt');
predictors = data(:,[7:9,13,14]);
data.good_all = zeros(length(data.good_oe211),1);
for i = 1:length(data.good_oe211)
    if sum(data{i,15:23})>0
        data.good_all(i) = 1;
    end
end
response = data.good_all;
predictorNames = {'slope','bottom_width', 'elevation', 'mannings','usgs_dist'};
template = templateTree('MaxNumSplits', 20);
classificationEnsemble = fitcensemble(...
    predictors, ...
    response, ...
    'Method', 'RUSBoost', ...
    'NumLearningCycles', 30, ...
    'Learners', template, ...
    'LearnRate', 0.1, ...
    'ClassNames', [0; 1]);
predictorExtractionFcn = @(t) t(:, predictorNames);
ensemblePredictFcn = @(x) predict(classificationEnsemble, x);
trainedClassifier.predictFcn = @(x) ensemblePredictFcn(predictorExtractionFcn(x));
trainedClassifier.ClassificationEnsemble = classificationEnsemble;
partitionedModel = crossval(trainedClassifier.ClassificationEnsemble, 'KFold', 25);
[validationPredictions, validationScores] = kfoldPredict(partitionedModel);
a = 1;
b = 1;
c = 1;
d = 1;
for i = 1:length(response)
    if response(i) == 0
        if validationPredictions(i) == 0
            trueneg(a) = i;
            a = a+1;
        else
            falsepos(b) = i;
            b = b+1;
        end
    else
        if validationPredictions(i) == 0
            falseneg(c) = i;
            c = c+1;
        else
            truepos(d) = i;
            d = d+1;
        end
    end
end

truepos_center = mean(predictors{truepos,:});
trueneg_center = mean(predictors{trueneg,:});
falsepos_center = mean(predictors{falsepos,:});
falseneg_center = mean(predictors{falseneg,:});

pca_data = table2array(data(:,[7:9,13,14]));
[coeff,score,sing_values] = pca(pca_data);
%coeff is the direction of the vectors for the principal components (orthogonal)
%score is the data matrix converted to orthogonal coordinates of coeff
%sing_values are the singular values that correspond to each principal component

k = 1;
s = 5;
for i = 1:s
    for j = 1:s
        figure(1)
        subplot(s,s,k)
        hold on
        scatter(score(trueneg,i),score(trueneg,j),'r')
        scatter(score(falsepos,i),score(falsepos,j),'rx')
        scatter(score(truepos,i),score(truepos,j),'b')
        scatter(score(falseneg,i),score(falseneg,j),'bx')
        hold off
        
        figure(2)
        subplot(s,s,k)
        hold on
        scatter(predictors{trueneg,i},predictors{trueneg,j},'r')
        scatter(predictors{falsepos,i},predictors{falsepos,j},'rx')
        scatter(predictors{truepos,i},predictors{truepos,j},'b')
        scatter(predictors{falseneg,i},predictors{falseneg,j},'bx')
        hold off
        
        k = k+1;
    end
end

%% post pca model
clear all,close all
data = readtable('classification_inputs/all.txt');
pre_predictors = data(:,[7:9,13,14]);
pca_data = table2array(data(:,[7:9,13,14]));
[coeff,score,sing_values] = pca(pca_data);
predictors = array2table(score);
data.good_all = zeros(length(data.good_oe211),1);
for i = 1:length(data.good_oe211)
    if sum(data{i,15:23})>0
        data.good_all(i) = 1;
    end
end
response = data.good_all;
predictorNames = {'component1','component2', 'component3', 'component4','component5'};
template = templateTree('MaxNumSplits', 20);
classificationEnsemble = fitcensemble(...
    predictors, ...
    response, ...
    'Method', 'RUSBoost', ...
    'NumLearningCycles', 30, ...
    'Learners', template, ...
    'LearnRate', 0.1, ...
    'ClassNames', [0; 1]);
predictorExtractionFcn = @(t) t(:, predictorNames);
ensemblePredictFcn = @(x) predict(classificationEnsemble, x);
trainedClassifier.predictFcn = @(x) ensemblePredictFcn(predictorExtractionFcn(x));
trainedClassifier.ClassificationEnsemble = classificationEnsemble;
partitionedModel = crossval(trainedClassifier.ClassificationEnsemble, 'KFold', 25);
[validationPredictions, validationScores] = kfoldPredict(partitionedModel);
a = 1;
b = 1;
c = 1;
d = 1;
for i = 1:length(response)
    if response(i) == 0
        if validationPredictions(i) == 0
            trueneg(a) = i;
            a = a+1;
        else
            falsepos(b) = i;
            b = b+1;
        end
    else
        if validationPredictions(i) == 0
            falseneg(c) = i;
            c = c+1;
        else
            truepos(d) = i;
            d = d+1;
        end
    end
end

truepos_center = mean(predictors{truepos,:});
trueneg_center = mean(predictors{trueneg,:});
falsepos_center = mean(predictors{falsepos,:});
falseneg_center = mean(predictors{falseneg,:});

k = 1;
s = 5;
for i = 1:s
    for j = 1:s
        figure(1)
        subplot(s,s,k)
        hold on
        scatter(score(trueneg,i),score(trueneg,j),'r')
        scatter(score(falsepos,i),score(falsepos,j),'rx')
        scatter(score(truepos,i),score(truepos,j),'b')
        scatter(score(falseneg,i),score(falseneg,j),'bx')
        hold off
        
        figure(2)
        subplot(s,s,k)
        hold on
        scatter(pre_predictors{trueneg,i},pre_predictors{trueneg,j},'r')
        scatter(pre_predictors{falsepos,i},pre_predictors{falsepos,j},'rx')
        scatter(pre_predictors{truepos,i},pre_predictors{truepos,j},'b')
        scatter(pre_predictors{falseneg,i},pre_predictors{falseneg,j},'bx')
        hold off
        
        k = k+1;
    end
end
%% kernel PCA
close all
% data = readtable('cluster_data.txt');
% pca_data = table2array(data(:,[7:9,22,25]));
test = pca_data';
K = zeros(size(test,2),size(test,2));
for row = 1:size(test,2)
    for col = 1:row
        temp = sum(((test(:,row) - test(:,col)).^2));
        K(row,col) = exp(-temp); % sigma = 1
%         K(row,col) = (test(:,row)'*test(:,col) + 1)^2;

    end
end

K = K + K'; 
% Dividing the diagonal element by 2 since it has been added to itself
for row = 1:size(test,2)
    K(row,row) = K(row,row)/2;
end
% We know that for PCA the data has to be centered. Even if the input data
% set 'X' lets say in centered, there is no gurantee the data when mapped
% in the feature space [phi(x)] is also centered. Since we actually never
% work in the feature space we cannot center the data. To include this
% correction a pseudo centering is done using the Kernel.
one_mat = ones(size(K));
V = K - one_mat*K - K*one_mat + one_mat*K*one_mat;
clear K

[coeff,eig_values] = pcacov(V);

k = 1;
s = 5;
for i = 1:s
    for j = 1:s
        subplot(s,s,k)
        hold on
%         scatter(coeff(trueneg,i),coeff(trueneg,j),'k')
%         scatter(coeff(truepos,i),coeff(truepos,j),'g')
%         scatter(coeff(falseneg,i),coeff(falseneg,j),'ro')
%         scatter(coeff(falsepos,i),coeff(falsepos,j),'rx')
        scatter(coeff(response==0,i),coeff(response==0,j),'r')
        scatter(coeff(response==1,i),coeff(response==1,j),'g')
        hold off
        
        k = k+1;
    end
end