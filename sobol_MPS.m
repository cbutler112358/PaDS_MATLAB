%% Testing convergence of Azzini et al. Sobol estimators; computing and 
% plotting convergence for total and first-order indices. Bootstrapping 
% is also performed to compute confidence intervals. Finally, Sobol
% estimators for the simulations using beta distributed priors are also
% computed. 
%
% Data files for each gene drive are named similarly, e.g. "VBSA_TLU_data_array_4e4_v5.mat"
% for two-locus underdominance, "VBSA_TLUTH_data_array_4e4_v5.mat" for 
% tethered homing, etc.
%
% In the case of tethered homing, for the beta distributed parameter
% values, it was necessary to "manually parallelize" computations so that
% the .mat files are split up using the naming convention
% "VBSA_TLUTH_data_array_5e3_beta_rep1.mat", ...,
% "VBSA_TLUTH_data_array_5e3_beta_rep8.mat". This is handled specifically
% in the final section. 

clear
close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% calculate estimators %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data files for the other gene drives are similarly named:
% two-locus underdominance (TLU), tethered homing (TLUTH), TARE, and 
% TADE
load("VBSA_TLU_data_array_4e4_v5.mat")
data_matrix     = data_array.allele_freq_4yr; 
% load("VBSA_TARE_data_array_3e4_v6_mk2.mat")
% data_matrix     = data_array.allele_freq_6yr;
% no. of samples 
N               = data_array.N;
sample_matrix   = data_array.sample_matrix;

% construct data matrices
A(:,:,1)   = data_matrix(1:N,:);           % A 
A(:,:,2)   = data_matrix((2*N+1):(3*N),:); % A1, used in fitness cost
A(:,:,3)   = data_matrix((3*N+1):(4*N),:); % A2, used in SDD
A(:,:,4)   = data_matrix((4*N+1):(5*N),:); % A3, used in LDD

B(:,:,1)   = data_matrix((N+1):(2*N),:);   % B
B(:,:,2)   = data_matrix((5*N+1):(6*N),:); % B1
B(:,:,3)   = data_matrix((6*N+1):(7*N),:); % B2
B(:,:,4)   = data_matrix((7*N+1):(8*N),:); % B3

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% calculate indices  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% sample vec
sample_vec = 1:N; % 1:100:N;

% measure convergence
AZZ2021_FO = nan(length(sample_vec),3,3);
AZZ2021_TO = nan(length(sample_vec),3,3);
% measure variance
var_mat    = nan(length(sample_vec),3);
for j = 2:4
    disp(j)
    for i = 1:length(sample_vec)
        % S_1_AZZ2021_vec(i,:)    = 2*sum((B1(1:i,:) - B(1:i,:)).*(A(1:i,:) - A1(1:i,:))) ./ sum((A(1:i,:) - B(1:i,:)).^2 + (B1(1:i,:) - A1(1:i,:)).^2);
        % ST_1_AZZ2021_vec(i,:)   = sum((B(1:i,:) - B1(1:i,:)).^2 + (A(1:i,:) - A1(1:i,:)).^2) ./ sum((A(1:i,:) - B(1:i,:)).^2 + (B1(1:i,:) - A1(1:i,:)).^2);
        AZZ2021_FO(i,:,j-1)    = 2*sum((B(1:sample_vec(i),:,j) - B(1:sample_vec(i),:,1)).*(A(1:sample_vec(i),:,1) - A(1:sample_vec(i),:,j))) ./ sum((A(1:sample_vec(i),:,1) - B(1:sample_vec(i),:,1)).^2 + (B(1:sample_vec(i),:,j) - A(1:sample_vec(i),:,j)).^2);
        AZZ2021_TO(i,:,j-1)    = sum((B(1:sample_vec(i),:,1) - B(1:sample_vec(i),:,j)).^2 + (A(1:sample_vec(i),:,1) - A(1:sample_vec(i),:,j)).^2) ./ sum((A(1:sample_vec(i),:,1) - B(1:sample_vec(i),:,1)).^2 + (B(1:sample_vec(i),:,j) - A(1:sample_vec(i),:,j)).^2);
    end
end

for i = 1:length(sample_vec)
    var_mat(i,:)    = (1/(2*sample_vec(i)))*sum((A(1:sample_vec(i),:,1) - B(1:sample_vec(i),:,1)).^2);
end


%% Plot convergence of Azzini et al. indices computed in the section above.

%%%%%%%%%%%%
%%% plot %%%
%%%%%%%%%%%%

fig = figure; 
% plot comparison of first order indices 
subplot(231)
plot(AZZ2021_FO(:,1,1),'linewidth',1.5)
hold on
plot(AZZ2021_FO(:,1,2),'linewidth',1.5)
plot(AZZ2021_FO(:,1,3),'linewidth',1.5)
title('(a)')
xticks([1,10000,20000,30000])
xticklabels([0,1,2,3])
legend('fitness cost', 'SDD', 'LDD')
% xlabel('patch 1')
ylim([0,1])
set(gca,'fontsize',16)

subplot(232)
plot(AZZ2021_FO(:,2,1),'linewidth',1.5)
hold on
plot(AZZ2021_FO(:,2,2),'linewidth',1.5)
plot(AZZ2021_FO(:,2,3),'linewidth',1.5)
title('(b)')
xticks([1,10000,20000,30000])
xticklabels([0,1,2,3])
% title('First order indices, fitness cost')
% xlabel('patch 2')
ylim([0,1])
set(gca,'fontsize',16)

subplot(233)
plot(AZZ2021_FO(:,3,1),'linewidth',1.5)
hold on
plot(AZZ2021_FO(:,3,2),'linewidth',1.5)
plot(AZZ2021_FO(:,3,3),'linewidth',1.5)
title('(c)')
xticks([1,10000,20000,30000])
xticklabels([0,1,2,3])
% xlabel('patch 3')
ylim([0,1])
set(gca,'fontsize',16)

pos = [0, 20, 733, 528];
set(gcf, 'Position', pos); %// gives x left, y bottom, width, height

% figure

% plot comparison of total order indices 
subplot(234)
plot(AZZ2021_TO(:,1,1),'linewidth',1.5)
hold on
plot(AZZ2021_TO(:,1,2),'linewidth',1.5)
plot(AZZ2021_TO(:,1,3),'linewidth',1.5)
% plot(var_mat(:,1),'linewidth',1.5)
% legend('fitness cost', 'SDD', 'LDD')
title('(d)')
xticks([1,10000,20000,30000])
xticklabels([0,1,2,3])
% xlabel('patch 1')
ylim([0,1])
set(gca,'fontsize',16)

subplot(235)
plot(AZZ2021_TO(:,2,1),'linewidth',1.5)
hold on
plot(AZZ2021_TO(:,2,2),'linewidth',1.5)
plot(AZZ2021_TO(:,2,3),'linewidth',1.5)
% plot(var_mat(:,2),'linewidth',1.5)
title('(e)')
xticks([1,10000,20000,30000])
xticklabels([0,1,2,3])
% title('Total order indices, fitness cost')
% xlabel('patch 2')
ylim([0,1])
set(gca,'fontsize',16)

subplot(236)
plot(AZZ2021_TO(:,3,1),'linewidth',1.5)
hold on
plot(AZZ2021_TO(:,3,2),'linewidth',1.5)
plot(AZZ2021_TO(:,3,3),'linewidth',1.5)
% plot(var_mat(:,3),'linewidth',1.5)
% xlabel('patch 3')
title('(f)')
xticks([1,10000,20000,30000])
xticklabels([0,1,2,3])
ylim([0,1])

set(gca,'fontsize',16)

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'sensitivity indices','interpreter','latex');
xlabel(han,'simulation number $(10^4)$','interpreter','latex');
han.XLabel.Position(2) = -0.075; % han.XLabel.Position(2) + 0.3;
set(gca,'fontsize',16)

pos = [0, 20, 733, 528];
set(gcf, 'Position', pos); %// gives x left, y bottom, width, height

%% [v5] Bootstrapping estimates provided above

clear
close

% distribution size used for bootstrap estimates
numReps = 1000; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% calculate estimators %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load("VBSA_TLU_data_array_4e4_v5.mat")
data_matrix     = data_array.allele_freq_4yr; 
% load("VBSA_TARE_data_array_3e4_v6_mk2.mat")
% data_matrix     = data_array.allele_freq_6yr;
% no. of samples 
N               = data_array.N;
sample_matrix   = data_array.sample_matrix;
% "unpack" data matrices
OG_A(:,:,1)   = data_matrix(1:N,:);           % A 
OG_A(:,:,2)   = data_matrix((2*N+1):(3*N),:); % A1, used in fitness cost
OG_A(:,:,3)   = data_matrix((3*N+1):(4*N),:); % A2, used in SDD
OG_A(:,:,4)   = data_matrix((4*N+1):(5*N),:); % A3, used in LDD

OG_B(:,:,1)   = data_matrix((N+1):(2*N),:);   % B
OG_B(:,:,2)   = data_matrix((5*N+1):(6*N),:); % B1
OG_B(:,:,3)   = data_matrix((6*N+1):(7*N),:); % B2
OG_B(:,:,4)   = data_matrix((7*N+1):(8*N),:); % B3

% arrays to store distributions; numReps x numParams x numPatch
FO_boot_mat = nan(numReps, 3, 3);
TO_boot_mat = nan(numReps, 3, 3);

for ii = 1:numReps
    disp(ii)
    % construct new data matrices by sampling original data matrix with
    % replacement
    A_sample = randsample(1:N, N, true);
    A(:,:,1) = OG_A(A_sample,:,1);
    A(:,:,2) = OG_A(A_sample,:,2);
    A(:,:,3) = OG_A(A_sample,:,3);
    A(:,:,4) = OG_A(A_sample,:,4);
    B(:,:,1) = OG_B(A_sample,:,1);
    B(:,:,2) = OG_B(A_sample,:,2);
    B(:,:,3) = OG_B(A_sample,:,3);
    B(:,:,4) = OG_B(A_sample,:,4);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% calculate indices  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
        
    % measure variance
    % var_mat    = nan(length(sample_vec),3);
    for j = 2:4
        % AZZ2021_FO(i,:,j-1)    = 2*sum((B(1:sample_vec(i),:,j) - B(1:sample_vec(i),:,1)).*(A(1:sample_vec(i),:,1) - A(1:sample_vec(i),:,j))) ./ sum((A(1:sample_vec(i),:,1) - B(1:sample_vec(i),:,1)).^2 + (B(1:sample_vec(i),:,j) - A(1:sample_vec(i),:,j)).^2);
        % AZZ2021_TO(i,:,j-1)    = sum((B(1:sample_vec(i),:,1) - B(1:sample_vec(i),:,j)).^2 + (A(1:sample_vec(i),:,1) - A(1:sample_vec(i),:,j)).^2) ./ sum((A(1:sample_vec(i),:,1) - B(1:sample_vec(i),:,1)).^2 + (B(1:sample_vec(i),:,j) - A(1:sample_vec(i),:,j)).^2);
        FO_boot_mat(ii,:,j-1)    = 2*sum((B(1:N,:,j) - B(1:N,:,1)).*(A(1:N,:,1) - A(1:N,:,j))) ./ sum((A(1:N,:,1) - B(1:N,:,1)).^2 + (B(1:N,:,j) - A(1:N,:,j)).^2);
        TO_boot_mat(ii,:,j-1)    = sum((B(1:N,:,1) - B(1:N,:,j)).^2 + (A(1:N,:,1) - A(1:N,:,j)).^2) ./ sum((A(1:N,:,1) - B(1:N,:,1)).^2 + (B(1:N,:,j) - A(1:N,:,j)).^2);
    end
    
    % var_mat(i,:)    = (1/(2*N))*sum((A(1:N,:,1) - B(1:N,:,1)).^2);

end

% construct vectors for each 

% e.g., percentiles of fitness cost for patch 1:
[quantile(FO_boot_mat(:,1,1),0.025), quantile(FO_boot_mat(:,1,1),0.975)]
% and for patch 2:
[quantile(FO_boot_mat(:,2,1),0.025), quantile(FO_boot_mat(:,2,1),0.975)]
% and patch 3
[quantile(FO_boot_mat(:,3,1),0.025), quantile(FO_boot_mat(:,3,1),0.975)]

% and SDD for patch 1
[quantile(FO_boot_mat(:,1,2),0.025), quantile(FO_boot_mat(:,1,2),0.975)]
% patch 2
[quantile(FO_boot_mat(:,2,2),0.025), quantile(FO_boot_mat(:,2,2),0.975)]
% patch 3
[quantile(FO_boot_mat(:,3,2),0.025), quantile(FO_boot_mat(:,3,2),0.975)]

% and LDD for patch 1
[quantile(FO_boot_mat(:,1,3),0.025), quantile(FO_boot_mat(:,1,3),0.975)]*100
% patch 2
[quantile(FO_boot_mat(:,2,3),0.025), quantile(FO_boot_mat(:,2,3),0.975)]*100
% patch 3
[quantile(FO_boot_mat(:,3,3),0.025), quantile(FO_boot_mat(:,3,3),0.975)]*100

%% total order bootstraps

[quantile(TO_boot_mat(:,1,1),0.025), quantile(TO_boot_mat(:,1,1),0.975)]
% and for patch 2:
[quantile(TO_boot_mat(:,2,1),0.025), quantile(TO_boot_mat(:,2,1),0.975)]
% and patch 3
[quantile(TO_boot_mat(:,3,1),0.025), quantile(TO_boot_mat(:,3,1),0.975)]

% and SDD for patch 1
[quantile(TO_boot_mat(:,1,2),0.025), quantile(TO_boot_mat(:,1,2),0.975)]
% patch 2
[quantile(TO_boot_mat(:,2,2),0.025), quantile(TO_boot_mat(:,2,2),0.975)]
% patch 3
[quantile(TO_boot_mat(:,3,2),0.025), quantile(TO_boot_mat(:,3,2),0.975)]

% and LDD for patch 1
[quantile(TO_boot_mat(:,1,3),0.025), quantile(TO_boot_mat(:,1,3),0.975)]*100
% patch 2
[quantile(TO_boot_mat(:,2,3),0.025), quantile(TO_boot_mat(:,2,3),0.975)]*100
% patch 3
[quantile(TO_boot_mat(:,3,3),0.025), quantile(TO_boot_mat(:,3,3),0.975)]*100

%% Sobol indices of the data using non-uniform dispersal priors (beta data)

clear
close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% calculate estimators %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load("VBSA_TADE_data_array_4e4_beta.mat")
% data_matrix     = data_array.allele_freq_4yr; 
load("VBSA_TARE_data_array_3e4_beta.mat")
data_matrix     = data_array.allele_freq_6yr;
% no. of samples 
N               = data_array.N;
sample_matrix   = data_array.sample_matrix;

% construct data matrices
A(:,:,1)   = data_matrix(1:N,:);           % A 
A(:,:,2)   = data_matrix((2*N+1):(3*N),:); % A1, used in fitness cost
A(:,:,3)   = data_matrix((3*N+1):(4*N),:); % A2, used in SDD
A(:,:,4)   = data_matrix((4*N+1):(5*N),:); % A3, used in LDD

B(:,:,1)   = data_matrix((N+1):(2*N),:);   % B
B(:,:,2)   = data_matrix((5*N+1):(6*N),:); % B1
B(:,:,3)   = data_matrix((6*N+1):(7*N),:); % B2
B(:,:,4)   = data_matrix((7*N+1):(8*N),:); % B3

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% calculate indices  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% sample vec
sample_vec = 1:N; % 1:100:N;

% measure convergence
AZZ2021_FO = nan(length(sample_vec),3,3);
AZZ2021_TO = nan(length(sample_vec),3,3);
% measure variance
var_mat    = nan(length(sample_vec),3);
for j = 2:4
    disp(j)
    for i = 1:length(sample_vec)
        % S_1_AZZ2021_vec(i,:)    = 2*sum((B1(1:i,:) - B(1:i,:)).*(A(1:i,:) - A1(1:i,:))) ./ sum((A(1:i,:) - B(1:i,:)).^2 + (B1(1:i,:) - A1(1:i,:)).^2);
        % ST_1_AZZ2021_vec(i,:)   = sum((B(1:i,:) - B1(1:i,:)).^2 + (A(1:i,:) - A1(1:i,:)).^2) ./ sum((A(1:i,:) - B(1:i,:)).^2 + (B1(1:i,:) - A1(1:i,:)).^2);
        AZZ2021_FO(i,:,j-1)    = 2*sum((B(1:sample_vec(i),:,j) - B(1:sample_vec(i),:,1)).*(A(1:sample_vec(i),:,1) - A(1:sample_vec(i),:,j))) ./ sum((A(1:sample_vec(i),:,1) - B(1:sample_vec(i),:,1)).^2 + (B(1:sample_vec(i),:,j) - A(1:sample_vec(i),:,j)).^2);
        AZZ2021_TO(i,:,j-1)    = sum((B(1:sample_vec(i),:,1) - B(1:sample_vec(i),:,j)).^2 + (A(1:sample_vec(i),:,1) - A(1:sample_vec(i),:,j)).^2) ./ sum((A(1:sample_vec(i),:,1) - B(1:sample_vec(i),:,1)).^2 + (B(1:sample_vec(i),:,j) - A(1:sample_vec(i),:,j)).^2);
    end
end

for i = 1:length(sample_vec)
    var_mat(i,:)    = (1/(2*sample_vec(i)))*sum((A(1:sample_vec(i),:,1) - B(1:sample_vec(i),:,1)).^2);
end

%% Combine *.mat files for TLUTH VBSA simulations using beta distribution, 
% and then compute Sobol indices     

clear

i = 1;
filename = "VBSA_TLUTH_data_array_5e3_beta_rep" + i + ".mat";

% load the file and store the data
load(filename);
data_matrix     = data_array.allele_freq_4yr;

% no. of samples per data file
N = 5000;
% construct data matrices
A(:,:,1)   = data_matrix(1:N,:);           % A 
A(:,:,2)   = data_matrix((2*N+1):(3*N),:); % A1, used in fitness cost
A(:,:,3)   = data_matrix((3*N+1):(4*N),:); % A2, used in SDD
A(:,:,4)   = data_matrix((4*N+1):(5*N),:); % A3, used in LDD

B(:,:,1)   = data_matrix((N+1):(2*N),:);   % B
B(:,:,2)   = data_matrix((5*N+1):(6*N),:); % B1
B(:,:,3)   = data_matrix((6*N+1):(7*N),:); % B2
B(:,:,4)   = data_matrix((7*N+1):(8*N),:); % B3

for i = 2:8

    filename = "VBSA_TLUTH_data_array_5e3_beta_rep" + i + ".mat";
    data_matrix     = data_array.allele_freq_4yr;

    new_A = [];
    new_B = []; 

    % cat data matrices
    new_A(:,:,1)   = cat(1,A(:,:,1),data_matrix(1:N,:));
    new_A(:,:,2)   = cat(1,A(:,:,2),data_matrix((2*N+1):(3*N),:)); % A1, used in fitness cost
    new_A(:,:,3)   = cat(1,A(:,:,3),data_matrix((3*N+1):(4*N),:)); % A2, used in SDD
    new_A(:,:,4)   = cat(1,A(:,:,4),data_matrix((4*N+1):(5*N),:)); % A3, used in LDD

    A = new_A; 

    new_B(:,:,1)   = cat(1,B(:,:,1),data_matrix((N+1):(2*N),:));
    new_B(:,:,2)   = cat(1,B(:,:,2),data_matrix((5*N+1):(6*N),:)); % A1, used in fitness cost
    new_B(:,:,3)   = cat(1,B(:,:,3),data_matrix((6*N+1):(7*N),:)); % A2, used in SDD
    new_B(:,:,4)   = cat(1,B(:,:,4),data_matrix((7*N+1):(8*N),:)); % A3, used in LDD

    B = new_B;     

end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% calculate indices  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% sample vec
sample_vec = 1:40000; % 1:100:N;

% measure convergence
AZZ2021_FO = nan(length(sample_vec),3,3);
AZZ2021_TO = nan(length(sample_vec),3,3);
% measure variance
var_mat    = nan(length(sample_vec),3);
for j = 2:4
    disp(j)
    for i = 1:length(sample_vec)
        % S_1_AZZ2021_vec(i,:)    = 2*sum((B1(1:i,:) - B(1:i,:)).*(A(1:i,:) - A1(1:i,:))) ./ sum((A(1:i,:) - B(1:i,:)).^2 + (B1(1:i,:) - A1(1:i,:)).^2);
        % ST_1_AZZ2021_vec(i,:)   = sum((B(1:i,:) - B1(1:i,:)).^2 + (A(1:i,:) - A1(1:i,:)).^2) ./ sum((A(1:i,:) - B(1:i,:)).^2 + (B1(1:i,:) - A1(1:i,:)).^2);
        AZZ2021_FO(i,:,j-1)    = 2*sum((B(1:sample_vec(i),:,j) - B(1:sample_vec(i),:,1)).*(A(1:sample_vec(i),:,1) - A(1:sample_vec(i),:,j))) ./ sum((A(1:sample_vec(i),:,1) - B(1:sample_vec(i),:,1)).^2 + (B(1:sample_vec(i),:,j) - A(1:sample_vec(i),:,j)).^2);
        AZZ2021_TO(i,:,j-1)    = sum((B(1:sample_vec(i),:,1) - B(1:sample_vec(i),:,j)).^2 + (A(1:sample_vec(i),:,1) - A(1:sample_vec(i),:,j)).^2) ./ sum((A(1:sample_vec(i),:,1) - B(1:sample_vec(i),:,1)).^2 + (B(1:sample_vec(i),:,j) - A(1:sample_vec(i),:,j)).^2);
    end
end

for i = 1:length(sample_vec)
    var_mat(i,:)    = (1/(2*sample_vec(i)))*sum((A(1:sample_vec(i),:,1) - B(1:sample_vec(i),:,1)).^2);
end


