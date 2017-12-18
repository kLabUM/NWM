clear all
close all
cd('~/Google Drive/Docs Kevin/National Water Model/NWM/Scripts');
load('mar29.mat');
% load('may09.mat');
load('dec5.mat');
load('granger.mat');
train = dec5;
test = mar29;
% test = Feb25toMay09;
granger_ids = granger;
clearvars mar29 dec5 granger Feb25toMay09
uniq_id = unique(train.ifis_id);
data = readtable('IFIS_NWM_data_new.txt');
diff = setdiff(data.ifis_id,uniq_id);
ind = find(data.ifis_id~=diff);
data = data(ind,:);

%% learn state space models
data.lags_train = zeros(length(uniq_id),1);
data.covs_train = zeros(length(uniq_id),1);
data.lags_test = zeros(length(uniq_id),1);
data.covs_test = zeros(length(uniq_id),1);
data.granger = zeros(length(uniq_id),1);
data.dynamics = zeros(length(uniq_id),1);
data.tf21 = zeros(length(uniq_id),1);
data.tf31 = zeros(length(uniq_id),1);
data.tf32 = zeros(length(uniq_id),1);
data.ss1 = zeros(length(uniq_id),1);
data.ss2 = zeros(length(uniq_id),1);
data.ss3 = zeros(length(uniq_id),1);
data.ss4 = zeros(length(uniq_id),1);
data.oe221 = zeros(length(uniq_id),1);
data.oe211 = zeros(length(uniq_id),1);
granger_ids = categorical(granger_ids.site(granger_ids.value==1));
usgs = readtable('USGS_sites.txt');
usgs = usgs(2:end,:);
usgs_lat = cellfun(@str2double,usgs.dec_lat_va);
usgs_lon = cellfun(@str2double,usgs.dec_long_va);
data.usgs_dist = ones(length(data.lat),1);
for i = 1:length(data.ifis_id)
    ind = find(train.ifis_id==data.ifis_id(i));
    data.dynamics(i) = train.dynamics(ind(1));
    R = 6371; %radius of earth
    lat_diff = (usgs_lat - data.lat(i))*pi()/180;
    lon_diff = (usgs_lon - data.lon(i))*pi()/180;
    a = sin(lat_diff/2).*sin(lat_diff/2) + cos(data.lat(i)*pi()/180).*cos(usgs_lat*pi()/180).*sin(lon_diff/2).*sin(lon_diff/2);
    c = 2*atan2(a.^(1/2),(1-a).^(1/2));
    d = R*c;
    data.usgs_dist(i) = min(d);
%     ind = find(train.ifis_id=='CEDARRV03');


    y = train.gage_height(ind);
    y(y<-2) = NaN;
    y = fillmissing(y,'linear');
    u = train.t0(ind);
    
    ind = find(test.ifis_id==uniq_id(i));
%     ind = find(test.ifis_id=='CEDARRV03');
    y_test = test.gage_height(ind);
    y_test(y_test<-2) = NaN;
    y_test = fillmissing(y_test,'linear');
    u_test = test.t0(ind);

    [c,lg] = xcov(y,u,100,'coeff');
    c = c(101:201); %only consider non-neg lags
    lg = lg(101:201);
    [data.covs_train(i),ind] = max(c); 
    data.lags_train(i) = lg(ind);
    
    [c_test,lg_test] = xcov(y_test,u_test,100,'coeff');
    c_test = c_test(101:201); %only consider non-neg lags
    lg_test = lg_test(101:201);
    [data.covs_test(i),ind] = max(c_test);
    data.lags_test(i) = lg_test(ind);
    
    if ismember(uniq_id(i),granger_ids)
        data.granger(i) = 1;
    end
    
    try
        toy_train = iddata(y,u,1);
        toy_train = misdata(toy_train);
        toy_train = detrend(toy_train);

        toy_test = iddata(y_test,u_test,1);
        toy_test = misdata(toy_test);
        toy_test = detrend(toy_test);
    catch
        data.tf21(i) = NaN;
        data.tf31(i) = NaN;
        data.tf32(i) = NaN;
        data.ss1(i) = NaN;
        data.ss2(i) = NaN;
        data.ss3(i) = NaN;
        data.ss4(i) = NaN;
        data.oe221(i) = NaN;
        data.oe211(i) = NaN;
        continue
    end
    
    % TF
    Options = tfestOptions;                                    
    Options.Display = 'off';                                    
    Options.WeightingFilter = [];                              
    Options.InitialCondition = 'backcast';                     
    try
        np = 2;                                                    
        nz = 1;                                                    
        num = arrayfun(@(x)NaN(1,x), nz+1,'UniformOutput',false);  
        den = arrayfun(@(x)[1, NaN(1,x)],np,'UniformOutput',false);                               
        iodValue = 0;                                              
        iodFree = true;                                            
        iodMin = 0;                                                
        iodMax = 30;                                               
        sysinit = idtf(num, den, 0);                               
        iod = sysinit.Structure.ioDelay;                           
        iod.Value = iodValue;                                      
        iod.Free = iodFree;                                        
        iod.Maximum = iodMax;                                      
        iod.Minimum = iodMin;                                      
        sysinit.Structure.ioDelay = iod;                           

        tf21 = tfest(toy_train, sysinit, Options);
        [~,fit,~] = compare(toy_test,tf21);
        data.tf21(i) = fit;
    catch
        data.tf21(i) = NaN;
    end
    
    try
        np = 3;                                                    
        nz = 1;                                                    
        num = arrayfun(@(x)NaN(1,x), nz+1,'UniformOutput',false);  
        den = arrayfun(@(x)[1, NaN(1,x)],np,'UniformOutput',false);
        sysinit = idtf(num, den, 0);                               
        iod = sysinit.Structure.ioDelay;                           
        iod.Value = iodValue;                                      
        iod.Free = iodFree;                                        
        iod.Maximum = iodMax;                                      
        iod.Minimum = iodMin;                                      
        sysinit.Structure.ioDelay = iod;           
        tf31 = tfest(toy_train, sysinit, Options);
        [~,fit,~] = compare(toy_test,tf31);
        data.tf31(i) = fit;
    catch
        data.tf31(i) = NaN;
    end
    
    try
        np = 3;                                                    
        nz = 2;                                                    
        num = arrayfun(@(x)NaN(1,x), nz+1,'UniformOutput',false);  
        den = arrayfun(@(x)[1, NaN(1,x)],np,'UniformOutput',false);
        sysinit = idtf(num, den, 0);                               
        iod = sysinit.Structure.ioDelay;                           
        iod.Value = iodValue;                                      
        iod.Free = iodFree;                                        
        iod.Maximum = iodMax;                                      
        iod.Minimum = iodMin;                                      
        sysinit.Structure.ioDelay = iod;           
        tf32 = tfest(toy_train, sysinit, Options);
        [~,fit,~] = compare(toy_test,tf32);
        data.tf32(i) = fit;
    catch
        data.tf32(i) = NaN;
    end
    
    %SS
    try
        ss1 = n4sid(toy_train, 1, 'DisturbanceModel', 'none', 'Ts', 0, Options);
        [~,fit,~] = compare(toy_test,ss1);
        data.ss1(i) = fit;
    catch
        data.ss1(i) = NaN;
    end
    
    try
        ss2 = n4sid(toy_train, 2, 'DisturbanceModel', 'none', 'Ts', 0, Options);
        [~,fit,~] = compare(toy_test,ss2);
        data.ss2(i) = fit;
    catch
        data.ss2(i) = NaN;
    end
    
    try
        ss3 = n4sid(toy_train, 3, 'DisturbanceModel', 'none', 'Ts', 0, Options);
        [~,fit,~] = compare(toy_test,ss3);
        data.ss3(i) = fit;
    catch
        data.ss3(i) = NaN;
    end
    
    try
        ss4 = n4sid(toy_train, 4, 'DisturbanceModel', 'none', 'Ts', 0, Options);
        [~,fit,~] = compare(toy_test,ss4);
        data.ss4(i) = fit;
    catch
        data.ss4(i) = NaN;
    end
    
    %OE
    Opt = oeOptions;                     
    Opt.InitialCondition = 'backcast';   
    Opt.Focus = 'simulation';            
    try
        nb = 2;                              
        nf = 2;                              
        nk = 1;                              
        oe221 = oe(toy_train,[nb nf nk], Opt);
        [~,fit,~] = compare(toy_test,oe221);
        data.oe221(i) = fit;
    catch
        data.oe221(i) = NaN;
    end
    
    try
        nb = 2;                              
        nf = 1;                              
        nk = 1;                              
        oe211 = oe(toy_train,[nb nf nk], Opt);
        [~,fit,~] = compare(toy_test,oe211);
        data.oe211(i) = fit;
    catch
        data.oe211(i) = NaN;
    end
    
    fprintf('%i/%i\r\n',i,length(uniq_id))
    
    
    
%     subplot(2,2,1)
%     plot(lg,c)
%     
%     subplot(2,2,3)
%     yyaxis left
%     plot(u)
%     yyaxis right
%     plot(y)
%     legend('Flow','Height')
%     
%     subplot(2,2,2)
%     plot(lg_test,c_test)
%     
%     subplot(2,2,4)
%     yyaxis left
%     plot(u_test)
%     yyaxis right
%     plot(y_test)
%     legend('Flow','Height')
%     
%     suptitle(sprintf('%s, Granger=%i',uniq_id(i),data.granger(i)))
    
%     pause
end

% save data table
%%
clear all, close all
cd('~/Google Drive/Docs Kevin/National Water Model/NWM');
lat = ncread('NWM_parameters/RouteLink_2016_04_07.nudgingOperational2016-04-08_chanparm3_mann_BtmWdth.nc','lat');
lon = ncread('NWM_parameters/RouteLink_2016_04_07.nudgingOperational2016-04-08_chanparm3_mann_BtmWdth.nc','lon');
ID = ncread('NWM_parameters/RouteLink_2016_04_07.nudgingOperational2016-04-08_chanparm3_mann_BtmWdth.nc','link');
to = ncread('NWM_parameters/RouteLink_2016_04_07.nudgingOperational2016-04-08_chanparm3_mann_BtmWdth.nc','to');
bottom_width = ncread('NWM_parameters/RouteLink_2016_04_07.nudgingOperational2016-04-08_chanparm3_mann_BtmWdth.nc','BtmWdth');
elevation = ncread('NWM_parameters/RouteLink_2016_04_07.nudgingOperational2016-04-08_chanparm3_mann_BtmWdth.nc','alt');
mannings = ncread('NWM_parameters/RouteLink_2016_04_07.nudgingOperational2016-04-08_chanparm3_mann_BtmWdth.nc','n');
slope = ncread('NWM_parameters/RouteLink_2016_04_07.nudgingOperational2016-04-08_chanparm3_mann_BtmWdth.nc','So');
musk_coeff = ncread('NWM_parameters/RouteLink_2016_04_07.nudgingOperational2016-04-08_chanparm3_mann_BtmWdth.nc','MusX');
order = ncread('NWM_parameters/RouteLink_2016_04_07.nudgingOperational2016-04-08_chanparm3_mann_BtmWdth.nc','order');
%rout_time = ncread('NWM_parameters/RouteLink_2016_04_07.nudgingOperational2016-04-08_chanparm3_mann_BtmWdth.nc','MusK');


states = shaperead('usastatelo' , 'UseGeoCoords' , true);
ind = inpolygon(lon,lat,states(15).Lon,states(15).Lat);

lat = lat(ind);
lon = lon(ind);
ID = ID(ind);
to = to(ind);
bottom_width = bottom_width(ind);
elevation = elevation(ind);
mannings = mannings(ind);
slope = slope(ind);
musk_coeff = musk_coeff(ind);
order = order(ind);
rout_time = rout_time(ind);

usgs = readtable('Scripts/USGS_sites.txt');
usgs = usgs(2:end,:);
usgs_lat = cellfun(@str2double,usgs.dec_lat_va);
usgs_lon = cellfun(@str2double,usgs.dec_long_va);
usgs_dist = ones(length(lat),1);
for i = 1:length(lat)
    R = 6371; %radius of earth
    lat_diff = (usgs_lat - lat(i))*pi()/180;
    lon_diff = (usgs_lon - lon(i))*pi()/180;
    a = sin(lat_diff/2).*sin(lat_diff/2) + cos(lat(i)*pi()/180).*cos(usgs_lat*pi()/180).*sin(lon_diff/2).*sin(lon_diff/2);
    c = 2*atan2(a.^(1/2),(1-a).^(1/2));
    d = R*c;
    usgs_dist(i,1) = min(d);
end

alldata = table(ID,lat,lon,to,bottom_width,elevation,mannings,slope,musk_coeff,order,rout_time,usgs_dist);
%save alldata
%%
clear all, close all
load alldata_may09

% load Scripts/cluster_data_2.mat
load Scripts/ss_data_TF_only.mat %change filename of data from first section as appropriate
data.ifis_id = string(data.ifis_id);
ifis_nwm = readtable('Scripts/IFIS_NWM_data_culled.txt');
k = 1;
omit_ind = [];
for i = 1:length(data.ifis_id)
    ind = find(ifis_nwm.ifis_id==data.ifis_id(i));
    if length(ind)>0
        nwm_id(k) = ifis_nwm.nwm_id(ind);
        k = k+1;
    else
        omit_ind = [omit_ind,i];
    end
end
data(omit_ind,:) = [];
data.nwm_id = nwm_id';

standard = alldata;
standard.order = double(standard.order);

standard.bottom_width = (standard.bottom_width - mean(standard.bottom_width))/std(standard.bottom_width);
standard.elevation = (standard.elevation - mean(standard.elevation))/std(standard.elevation);
standard.mannings = (standard.mannings - mean(standard.mannings))/std(standard.mannings);
standard.slope = (standard.slope - mean(standard.slope))/std(standard.slope);
standard.usgs_dist = (standard.usgs_dist - mean(standard.usgs_dist))/std(standard.usgs_dist);
standard.musk_coeff = (standard.musk_coeff - mean(standard.musk_coeff))/std(standard.musk_coeff);
standard.order = (standard.order - mean(standard.order))/std(standard.order);
% standard.rout_time = (standard.rout_time - mean(standard.rout_time))/std(standard.rout_time);

data.good_all = zeros(length(data.ifis_id),1);
threshold = 50;
% for i = 1:length(data.ifis_id)
%     if data.tf21(i)>threshold||data.tf31(i)>threshold||data.tf32(i)>threshold||data.ss1(i)>threshold||data.ss2(i)>threshold||data.ss3(i)>threshold||data.ss4(i)>threshold||data.oe221(i)>threshold||data.oe211(i)>threshold
%         data.good_all(i) = 1;
%     end
% end
for i = 1:length(data.ifis_id)
    if max(table2array(data(i,19:32)))>threshold
        data.good_all(i) = 1;
    end
end

%%
[~,ind1,ind2] = intersect(standard.ID,data.nwm_id);

pre_predictors = standard(:,[5:8,10,11]);
pca_data = table2array(pre_predictors);
[coeff,score,sing_values] = pca(pca_data);
predictors = array2table(score(ind1,:));
response = data.good_all(ind2);

good_ind = [];
bad_ind = [];
for i = 1:length(data.nwm_id)
    temp = find(alldata.ID==data.nwm_id(i));
    if data.good_all(i)==1
        good_ind = [good_ind;temp];
    else
        bad_ind = [bad_ind;temp];
    end
end

figure
bottom_width = NaN(length(bad_ind),2);
bottom_width(:,1) = pre_predictors.bottom_width(bad_ind);
bottom_width(1:length(good_ind),2) = pre_predictors.bottom_width(good_ind);

elevation = NaN(length(bad_ind),2);
elevation(:,1) = pre_predictors.elevation(bad_ind);
elevation(1:length(good_ind),2) = pre_predictors.elevation(good_ind);

mannings = NaN(length(bad_ind),2);
mannings(:,1) = pre_predictors.mannings(bad_ind);
mannings(1:length(good_ind),2) = pre_predictors.mannings(good_ind);

slope = NaN(length(bad_ind),2);
slope(:,1) = pre_predictors.slope(bad_ind);
slope(1:length(good_ind),2) = pre_predictors.slope(good_ind);

order = NaN(length(bad_ind),2);
order(:,1) = pre_predictors.order(bad_ind);
order(1:length(good_ind),2) = pre_predictors.order(good_ind);

usgs_dist = NaN(length(bad_ind),2);
usgs_dist(:,1) = pre_predictors.usgs_dist(bad_ind);
usgs_dist(1:length(good_ind),2) = pre_predictors.usgs_dist(good_ind);

h1 = cat(3, bottom_width', elevation', mannings',slope',order',usgs_dist');

aboxplot(h1,'labels',{'Channel Bottom Width','Elevation',"Manning's Roughness",'Stream Order','Slope','Distance to USGS Gauge'},'colorgrad','red_up');
legend('High Performing','Low Performing')
ylabel('Standard Deviations from Mean')
set(gca,'FontSize',20)

figure
good_ind = find(response==1);
bad_ind = find(response==0);

comp1 = NaN(length(bad_ind),2);
comp1(:,1) = predictors.Var1(bad_ind);
comp1(1:length(good_ind),2) = predictors.Var1(good_ind);

comp2 = NaN(length(bad_ind),2);
comp2(:,1) = predictors.Var2(bad_ind);
comp2(1:length(good_ind),2) = predictors.Var2(good_ind);

comp3 = NaN(length(bad_ind),2);
comp3(:,1) = predictors.Var3(bad_ind);
comp3(1:length(good_ind),2) = predictors.Var3(good_ind);

comp4 = NaN(length(bad_ind),2);
comp4(:,1) = predictors.Var4(bad_ind);
comp4(1:length(good_ind),2) = predictors.Var4(good_ind);

comp5 = NaN(length(bad_ind),2);
comp5(:,1) = predictors.Var5(bad_ind);
comp5(1:length(good_ind),2) = predictors.Var5(good_ind);

comp6 = NaN(length(bad_ind),2);
comp6(:,1) = predictors.Var6(bad_ind);
comp6(1:length(good_ind),2) = predictors.Var6(good_ind);

h2 = cat(3, comp1',comp2',comp3',comp4',comp5',comp6');

aboxplot(h2,'labels',{'Var1','Var2','Var3','Var4','Var5','Var6'},'colorgrad','red_up');
legend('High Performing','Low Performing')
ylabel('Standard Deviations from Mean')
set(gca,'FontSize',22)

% logistic regression model

beta = glmfit(table2array(predictors),response,'binomial');

% figure
% subplot(2,2,2)
% boxplot(table2array(pre_predictors(bad_ind,:)),'labels',{'Bottom Width','Elevation','Mannings','Order','Slope','USGS Distance'})
% title('Bad Sites')
% subplot(2,2,1)
% boxplot(table2array(pre_predictors(good_ind,:)),'labels',{'Bottom Width','Elevation','Mannings','Order','Slope','USGS Distance'})
% title('Good Sites')
% ylabel('Standardized Values')
% 
% good_ind = find(response==1);
% bad_ind = find(response==0);
% % figure
% subplot(2,2,4)
% boxplot(table2array(predictors(bad_ind,:)))
% % title('Bad Sites')
% subplot(2,2,3)
% boxplot(table2array(predictors(good_ind,:)))
% % title('Good Sites')
% ylabel('Standardized Values')


%%
predictorNames = {'Var1','Var2','Var3','Var4','Var5','Var6'};
% template = templateTree('MaxNumSplits', 25);
classificationEnsemble = fitcensemble(...
    predictors, ...
    response, ...
    'Method', 'LogitBoost', ...
    'NumLearningCycles', 500, ...
    'Learners', 'tree', ...
    'LearnRate', 0.1, ...
    'ClassNames', [0; 1]);
predictorExtractionFcn = @(t) t(:, predictorNames);
ensemblePredictFcn = @(x) predict(classificationEnsemble, x);
trainedClassifier.predictFcn = @(x) ensemblePredictFcn(predictorExtractionFcn(x));
trainedClassifier.ClassificationEnsemble = classificationEnsemble;
partitionedModel = crossval(trainedClassifier.ClassificationEnsemble, 'KFold', 25);
% partitionedModel = crossval(trainedClassifier.ClassificationEnsemble, 'leaveout','on');
[validationPredictions, validationScores] = kfoldPredict(partitionedModel);
[ROC_X,ROC_Y,T,AUC] = perfcurve(response,validationScores(:,2),1);
plot(ROC_X,ROC_Y)


% [check1,check2] = predict(classificationEnsemble,predictors);
% 
% confusion = confusionmat(response,double(check2(:,2)>0.7));
% [ROC_X,ROC_Y,T,AUC] = perfcurve(response,check2(:,2),1);
% 
% x = find(check1==1);
% y = find(response==1);
% truepos = intersect(x,y);
% y = find(response==0);
% falsepos = intersect(x,y);
% x = find(check1==0);
% trueneg = intersect(x,y);
% y = find(response==1);
% falseneg = intersect(x,y);
% 
% figure
% subplot(2,2,4)
% boxplot(table2array(predictors(truepos,:)));title(sprintf('True Positive: %i',length(truepos)))
% % ylim([-4.5,2.5])
% subplot(2,2,1)
% boxplot(table2array(predictors(trueneg,:)));title(sprintf('True Negative: %i',length(trueneg)))
% % ylim([-4.5,2.5])
% subplot(2,2,2)
% boxplot(table2array(predictors(falsepos,:)));title(sprintf('False Positive: %i',length(falsepos)))
% % ylim([-4.5,2.5])
% subplot(2,2,3)
% boxplot(table2array(predictors(falseneg,:)));title(sprintf('False Negative: %i',length(falseneg)))
% % ylim([-4.5,2.5])
% figure
% hold on
% scatter(check2(truepos,1),check2(truepos,2),'g')
% scatter(check2(trueneg,1),check2(trueneg,2),'r')
% scatter(check2(falsepos,1),check2(falsepos,2),'rx')
% scatter(check2(falseneg,1),check2(falseneg,2),'gx')
% legend('true pos','true neg', 'false pos','false neg')
% hold off

[prediction1,prediction2] = predict(classificationEnsemble,score);
% view(classificationEnsemble.Trained{10},'mode','graph')
%%
figure
hold on
states = shaperead('usastatelo' , 'UseGeoCoords' , true);
plot(states(15).Lon,states(15).Lat,'k')
count=0;
color = (prediction2(:,2)-min(prediction2(:,2)))/(max(prediction2(:,2))-min(prediction2(:,2)));
% % color = zeros(length(prediction2(:,1)),1);
% % ind = find(ROC_Y==1,1);
% % ind = find(prediction2(:,2)>T(ind));
% % color(ind) = 1;

% color = zeros(1,length(pca_data(:,1)));
% ind = find(pca_data(:,1)>=2);
% color(ind) = 1;

% color = glmval(beta,pca_data,'logit');
for i = 1:length(alldata.ID)
%     if inpolygon(alldata.lon(i),alldata.lat(i),states(15).Lon,states(15).Lat)
        to_ind = find(alldata.to==alldata.ID(i));
        if isempty(to_ind)==0
            for j=1:length(to_ind)
                plot([alldata.lon(i);alldata.lon(to_ind(j))],[alldata.lat(i);alldata.lat(to_ind(j))],'Color',[1-color(i) 0 color(i)],'LineWidth',alldata.order(i)/2)
    %             if prediction(i)==1
    %                 plot([alldata.lon(i);alldata.lon(to_ind(j))],[alldata.lat(i);alldata.lat(to_ind(j))],'r')
    %             else
    %                 plot([alldata.lon(i);alldata.lon(to_ind(j))],[alldata.lat(i);alldata.lat(to_ind(j))],'b')
    %             end
            end
        end
%     end
end
map = [linspace(1,0,100)',zeros(100,1),linspace(0,1,100)'];
colormap(map)
c = colorbar;
ylabel(c,'Probability of Performing Well')
set(gca,'FontSize',22)
hold off
%%
file = '/Users/kjfries/Google Drive/Docs Kevin/National Water Model/Figures/color_coded_streams_max_NSE.png';
export_fig(file,'-nocrop','-transparent','-r100',gcf)