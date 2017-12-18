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
% load('good.mat');
%%
% close all
% id_ind = [21,27:30,35,53,55:58,60,68,69,76,78,90,91,93,98,99,100,104,110,113,114,129,133,138,151,162,163,164,183,184,191,192,195,197]';
%%
% i = 22;
% ind = find(mar12.ifis_id==uniq_id(id_ind(i)));
% 
% y = mar12.gage_height(ind);
% y(y<-4) = NaN;
% u = mar12.t0(ind);
% 
% train = iddata(y(1:end-500),u(1:end-500),1);
% train = misdata(train);
% train = detrend(train);
% 
% test = iddata(y(end-499:end),u(end-499:end),1);
% test = misdata(test);
% test = detrend(test);
% 
% disp(uniq_id(id_ind(i)));
% plot(train);
% figure;
% plot(test);
%%
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
output = [];
%i=find(strcmp(data.ifis_id,'SUGARCR01'));
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
%     ind = find(train.ifis_id=='SUGARCR01');


    y = train.gage_height(ind);
    y(y<-2) = NaN;
    y = fillmissing(y,'linear');
    u = train.t0(ind);
    
    ind = find(test.ifis_id==data.ifis_id(i));
%     ind = find(test.ifis_id=='SUGARCR01');
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
    
    if ismember(data.ifis_id(i),granger_ids)
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
        [temp,fit,~] = compare(toy_test,tf21);
        output(:,1) = temp.OutputData;
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
        [temp,fit,~] = compare(toy_test,tf31);
        output(:,2) = temp.OutputData;
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
        [temp,fit,~] = compare(toy_test,tf32);
        output(:,3) = temp.OutputData;
        data.tf32(i) = fit;
    catch
        data.tf32(i) = NaN;
    end
    
    %SS
    try
        ss1 = n4sid(toy_train, 1, 'DisturbanceModel', 'none', 'Ts', 0, Options);
        [temp,fit,~] = compare(toy_test,ss1);
        output(:,4) = temp.OutputData;
        data.ss1(i) = fit;
    catch
        data.ss1(i) = NaN;
    end
    
    try
        ss2 = n4sid(toy_train, 2, 'DisturbanceModel', 'none', 'Ts', 0, Options);
        [temp,fit,~] = compare(toy_test,ss2);
        output(:,5) = temp.OutputData;
        data.ss2(i) = fit;
    catch
        data.ss2(i) = NaN;
    end
    
    try
        ss3 = n4sid(toy_train, 3, 'DisturbanceModel', 'none', 'Ts', 0, Options);
        [temp,fit,~] = compare(toy_test,ss3);
        output(:,6) = temp.OutputData;
        data.ss3(i) = fit;
    catch
        data.ss3(i) = NaN;
    end
    
    try
        ss4 = n4sid(toy_train, 4, 'DisturbanceModel', 'none', 'Ts', 0, Options);
        [temp,fit,~] = compare(toy_test,ss4);
        output(:,7) = temp.OutputData;
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
        [temp,fit,~] = compare(toy_test,oe221);
        output(:,8) = temp.OutputData;
        data.oe221(i) = fit;
    catch
        data.oe221(i) = NaN;
    end
    
    try
        nb = 2;                              
        nf = 1;                              
        nk = 1;                              
        oe211 = oe(toy_train,[nb nf nk], Opt);
        [temp,fit,~] = compare(toy_test,oe211);
        output(:,9) = temp.OutputData;
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
%     
%     pause
    for col=1:size(output,2)
        if max(output(:,col))>10*max(toy_test.OutputData)
            output(:,col) = nan(size(output(:,col)));
        end
    end
    meanline = nanmean(output,2);
    stdev = nanstd(output,0,2);
    
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(2,2,1)
    plot(u)
    title('Training Flows')
    
    subplot(2,2,3)
    plot(y)
    title('Training Heights')
    
    subplot(2,2,2)
    plot(u_test)
    title('Testing Flows')
    
    subplot(2,2,4)
    f = [meanline+2*stdev; flip(meanline-2*stdev)];
    temp = [1:length(meanline),flip(1:length(meanline))];
    fill(temp, f, [7 7 7]/8)
    hold on
    plot(toy_test.OutputData,'b')
    plot(meanline,'r')
    hold off
    title('Testing Heights')
    
    suptitle(sprintf('%s',string(data.ifis_id(i))));
    
    file = sprintf('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/Figures/%s',string(data.ifis_id(i)));
    saveas(gcf,file,'epsc')
    
    pause
    close all
end

%%

data.lags_train = zeros(length(uniq_id),1);
data.covs_train = zeros(length(uniq_id),1);
data.lags_test = zeros(length(uniq_id),1);
data.covs_test = zeros(length(uniq_id),1);
% data.granger = zeros(length(uniq_id),1);
data.dynamics = zeros(length(uniq_id),1);
data.tf01 = zeros(length(uniq_id),1);
data.tf11 = zeros(length(uniq_id),1);
data.tf02 = zeros(length(uniq_id),1);
data.tf12 = zeros(length(uniq_id),1);
data.tf22 = zeros(length(uniq_id),1);
data.tf03 = zeros(length(uniq_id),1);
data.tf13 = zeros(length(uniq_id),1);
data.tf23 = zeros(length(uniq_id),1);
data.tf33 = zeros(length(uniq_id),1);
data.tf04 = zeros(length(uniq_id),1);
data.tf14 = zeros(length(uniq_id),1);
data.tf24 = zeros(length(uniq_id),1);
data.tf34 = zeros(length(uniq_id),1);
data.tf44 = zeros(length(uniq_id),1);
% granger_ids = categorical(granger_ids.site(granger_ids.value==1));
usgs = readtable('USGS_sites.txt');
usgs = usgs(2:end,:);
usgs_lat = cellfun(@str2double,usgs.dec_lat_va);
usgs_lon = cellfun(@str2double,usgs.dec_long_va);
data.usgs_dist = ones(length(data.lat),1);
output = [];

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
%     ind = find(train.ifis_id=='SUGARCR01');


    y = train.gage_height(ind);
    y(y<-2) = NaN;
    y = fillmissing(y,'linear');
    u = train.t0(ind);
    
    ind = find(test.ifis_id==data.ifis_id(i));
%     ind = find(test.ifis_id=='SUGARCR01');
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
    
%     if ismember(data.ifis_id(i),granger_ids)
%         data.granger(i) = 1;
%     end
    
    try
        toy_train = iddata(y,u,1);
        toy_train = misdata(toy_train);
        toy_train = detrend(toy_train);

        toy_test = iddata(y_test,u_test,1);
        toy_test = misdata(toy_test);
        toy_test = detrend(toy_test);
    catch
        data.tf01(i) = NaN;
        data.tf11(i) = NaN;
        data.tf02(i) = NaN;
        data.tf12(i) = NaN;
        data.tf22(i) = NaN;
        data.tf03(i) = NaN;
        data.tf13(i) = NaN;
        data.tf23(i) = NaN;
        data.tf33(i) = NaN;
        data.tf04(i) = NaN;
        data.tf14(i) = NaN;
        data.tf24(i) = NaN;
        data.tf34(i) = NaN;
        data.tf44(i) = NaN;
        continue
    end
    
    % TF
    Options = tfestOptions;                                    
    Options.Display = 'off';                                    
    Options.WeightingFilter = [];                              
    Options.InitialCondition = 'backcast';
    k = 1;
    for np=0:4
        for nz=0:np
            try                                                   
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

                tf = tfest(toy_train, sysinit, Options);
                [temp,fit,~] = compare(toy_test,tf);
                output(:,k) = temp.OutputData;
                fits(k) = fit;
                k = k+1;
            catch
                fits(k) = NaN;
                k = k+1;
            end
        end
    end
    data(i,19:32) = fits;
    
    fprintf('%i/%i\r\n',i,length(data.ifis_id))
    
    
    
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
%     
%     pause
%     for col=1:size(output,2)
%         if max(output(:,col))>10*max(toy_test.OutputData)
%             output(:,col) = nan(size(output(:,col)));
%         end
%     end
%     meanline = nanmean(output,2);
%     stdev = nanstd(output,0,2);
%     
%     figure('units','normalized','outerposition',[0 0 1 1])
%     subplot(2,2,1)
%     plot(u)
%     title('Training Flows')
%     
%     subplot(2,2,3)
%     plot(y)
%     title('Training Heights')
%     
%     subplot(2,2,2)
%     plot(u_test)
%     title('Testing Flows')
%     
%     subplot(2,2,4)
%     f = [meanline+2*stdev; flip(meanline-2*stdev)];
%     temp = [1:length(meanline),flip(1:length(meanline))];
%     fill(temp, f, [7 7 7]/8)
%     hold on
%     plot(toy_test.OutputData,'b')
%     plot(meanline,'r')
%     hold off
%     title('Testing Heights')
%     
%     suptitle(sprintf('%s',string(data.ifis_id(i))));
%     
%     file = sprintf('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/Figures/TF_Only/%s',string(data.ifis_id(i)));
%     saveas(gcf,file,'epsc')
    
%     pause
%     close all
end