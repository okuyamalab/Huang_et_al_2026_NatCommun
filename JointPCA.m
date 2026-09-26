%% load data

% X_all (349x100)
% WT neurons(1-177)+KO neurons(178-349), self-shock(1-50)+average of other-shock(51-100) 
% 5Hz, -2-8 s (0-2 s self-shock/other-shock)
load('250807.mat')

%% joint subspace analysis
% This code is largely based on code from the following. 
% behret/paper_code_active_avoidance: analysis_code_v1 (Version v1). Zenodo. https://doi.org/10.5281/zenodo.11283463

%% prepare data for DR 
y_no_norm = X_all;
% cat subjects and normalize (mean subtract and normalize with Frobenius norm)
y = y_no_norm - mean(y_no_norm,2);
y = y./norm(y,'fro');
    
%% DR
% following the notation af Ani's paper:
% y: joint neural activities
% u: all coeffs from SVD
% u_i: coeffs per subject
% q: orthogonalized coeffs per subject 

% do pca
[u,s,v] = svd(y);
full_coeff = u;

% split up into subjects
q = {};
start_idx = 1;
cell_idx = {};
n_cells = [42 49 40 46 43 20 44 65];
n_pcs = 20;
for sub = 1:8
    % figure out indices of cells for this sub
    cell_idx{sub} = start_idx:start_idx+n_cells(sub)-1;
    start_idx = start_idx+n_cells(sub); 
    
    u_i = u(cell_idx{sub},1:n_pcs);
    [Q,R] = qr(u_i);
    q{sub} = Q(:,1:n_pcs);

    % qr decomposition flips some sings -> flip back
    ccs = [];
    for pc = 1:n_pcs
        cc = corrcoef(u_i(:,pc),q{sub}(:,pc));
        ccs(pc) = cc(1,2);
        if cc(1,2) < 0
            q{sub}(:,pc) = -q{sub}(:,pc);
        end
    end
end       

ve_real = calcVEall(y,full_coeff,n_pcs);
ve_sub_real = calcVEind(y,q,cell_idx,n_pcs);
dim_sim_real = calcDimSim(y_no_norm,q,cell_idx,n_pcs);

%% CONTROL: run PCA individually and save weigths
q_individual = {};
for sub = 1:8
    this_y = y(cell_idx{sub},:);
    [u,s,v] = svd(this_y);
    q_individual{sub} = u(:,1:n_pcs);
    
    % align sign of weigths such that q and q_ind have positive corr
    ccs = [];
    for pc = 1:n_pcs
        cc = corrcoef(q_individual{sub}(:,pc),q{sub}(:,pc));
        ccs(pc) = cc(1,2);
        if cc(1,2) < 0
            q_individual{sub}(:,pc) = -q_individual{sub}(:,pc);
        end
    end
end

ve_sub_ind = calcVEind(y,q_individual,cell_idx,n_pcs);
dim_sim_ind = calcDimSim(y_no_norm,q_individual,cell_idx,n_pcs);

%% CONTROL: shuffle time steps for every subject to destroy alignment info
rng(0)
for k = 1:1000 
    q_shuffle = {};
    y_shuffle = y;
    y_shuffle_no_norm = y_no_norm;

    % first shuffle time steps for all cells belonging to one subject
    for sub = 1:8
        rand_idx = randperm(size(y,2));        
        y_shuffle(cell_idx{sub},:) = y(cell_idx{sub},rand_idx);
        y_shuffle_no_norm(cell_idx{sub},:) = y_no_norm(cell_idx{sub},rand_idx);
    end
    
    % then run joint PCA and QR as above
    [u,s,v] = svd(y_shuffle);
    full_coeff_shuffle = u;

    for sub = 1:8
        u_i = u(cell_idx{sub},1:n_pcs);
        [Q,R] = qr(u_i);
        q_shuffle{sub} = Q(:,1:n_pcs);

        % qr decomposition flips some sings -> flip back
        ccs = [];
        for pc = 1:n_pcs
            cc = corrcoef(u_i(:,pc),q_shuffle{sub}(:,pc));
            ccs(pc) = cc(1,2);
            if cc(1,2) < 0
                q_shuffle{sub}(:,pc) = -q_shuffle{sub}(:,pc);
            end
        end

        org = y_shuffle(cell_idx{sub},:);
        for pc = 1:n_pcs
            rec = (org' * q_shuffle{sub}(:,pc) * q_shuffle{sub}(:,pc)')';
            ve_sub_shuffle(k,sub,pc) = 100*(1 - norm(org-rec,'fro')^2  / norm(org,'fro')^2);
        end
    end
    
    for pc = 1:n_pcs
        rec = (y_shuffle' * full_coeff_shuffle(:,pc) * full_coeff_shuffle(:,pc)')';
        ve_shuffle(pc,k) = 100*(1 - norm(y_shuffle-rec,'fro')^2  / norm(y_shuffle,'fro')^2);
    end

    this_y = y_shuffle_no_norm;
    this_q = q_shuffle;
    % calc cc between different PCs for all sub combinations
    cc_all_comb = [];
    for pc1 = 1:n_pcs
        for pc2 = 1:n_pcs
            n_comb = 1;
            for sub1 = 1:8
                sub1_pc = this_y(cell_idx{sub1},:)' * this_q{sub1}(:,pc1);
                for sub2 = sub1+1:8
                    sub2_pc = this_y(cell_idx{sub2},:)' * this_q{sub2}(:,pc2);
                    cc = corrcoef(sub1_pc,sub2_pc);
                    cc_all_comb(pc1,pc2,n_comb) = cc(1,2);
                    n_comb = n_comb+1;
                end
            end
        end
    end
    dim_sim_shuffle(:,:,k) = mean(cc_all_comb,3);
    disp(k)
end

%% SFig.7a
figure
plot(1:20,cumsum(ve_real),'k')
hold on
plot(1:20,mean(cumsum(ve_shuffle),2),'b')
legend('PCA+QR', 'Shuffle PCA+QR')
xlim([1 20])
ylim([0 100])
xlabel('Dimensions')
ylabel('VE (%)')

%% SFig.7b
figure
plot(1:20,cumsum(mean(ve_sub_real,1)),'k')  
hold on
plot(1:20,cumsum(mean(ve_sub_ind,1)),'r')
hold on
mean_ve_shuffle = cumsum(mean(mean(ve_sub_shuffle,1),2));
plot(1:20,mean_ve_shuffle(:),'b')
legend('PCA+QR', 'Shuffle PCA+QR')
xlim([1 20])
ylim([0 100])
xlabel('Dimensions')
ylabel('VE (%)')

%% SFig.7c
dim_sim = zeros(20,2);
for i= 1:20
    dim_sim(i,1) = dim_sim_real(i,i);
    dim_sim(i,2) = dim_sim_ind(i,i);
end
dim_sim_s = zeros(20,1000);
for k = 1:1000
    for i= 1:20
        dim_sim_s(i,k) = dim_sim_shuffle(i,i,k);
    end
end
figure
plot(1:20, dim_sim(:,1),'k')
hold on
plot(1:20, dim_sim(:,2),'r')
hold on
% plot(1:20, mean(dim_sim_s,2),'b')
shadedErrorBar(1:20, mean(dim_sim_s,2), std(dim_sim_s,[],2),'lineProps','b')
xlim([1 20])
ylabel('Mean Dim. Similarity')
xlabel('Dimensions')
legend('PCA+QR', 'ind PCA', 'Shuffle PCA+QR')

%% calculate PCA-QR projections for all subs 
% (for plotting projections and example of dims-sim measure)
% actual dim-sim is calculated above for all 3 cases (PCA-QR and 2 controls)

projs = [];
for sub = 1:8
    projs(sub,:,:) = y_no_norm(cell_idx{sub},:)' * q{sub};
end

%% Fig.3f
for i = 1:4
figure
    plot(1:50,mean(projs(1:4,51:100,i),1)','b')
    shadedErrorBar(1:50, mean(projs(1:4,51:100,i),1), std(projs(1:4,51:100,i),1)/sqrt(4),'lineProps','b')
%     plot(1:50,mean(projs(1:4,1:50,i),1)','r')
%     shadedErrorBar(1:50, mean(projs(1:4,1:50,i),1), std(projs(1:4,1:50,i),1)/sqrt(4),'lineProps','r')
hold on
    plot(1:50,mean(projs(5:8,51:100,i),1)','k')
    shadedErrorBar(1:50, mean(projs(5:8,51:100,i),1), std(projs(5:8,51:100,i),1)/sqrt(4),'lineProps','k')
%     plot(1:50,mean(projs(5:8,1:50,i),1)','k')
%     shadedErrorBar(1:50, mean(projs(5:8,1:50,i),1), std(projs(5:8,1:50,i),1)/sqrt(4),'lineProps','k')
xline(10.5)
xline(20.5)
xlim([1 50])
%     ylim([-7 4]) % self-shock
     ylim([-0.7 0.5]) % other-shock
end

%% AUC diff 

auc_diff = zeros(8,2,20);
% first identify up or down by comparing mean value of -2~0s vs 0~2s
for ii = 1:8
    for k = 1:20
        auc_diff(ii,1,k) = trapz(projs(ii,11:20,k))-trapz(projs(ii,1:10,k));
        auc_diff(ii,2,k) = trapz(projs(ii,61:70,k))-trapz(projs(ii,51:60,k));
    end
end

%　calculate p-value of corr between WT and KO
corr_diff_all = zeros(10000,1);
rng default
for k = 1:10000
    auc_tmp = zeros(8,2,20);
    r_arr = randperm(8);
    for t = 1:20   
        auc_tmp(1:4,:,t) = auc_diff(r_arr(1:4),:,t);
        auc_tmp(5:8,:,t) = auc_diff(r_arr(5:8),:,t);
    end    
    corr_WT = corr(reshape(auc_tmp(1:4,1,:),[],1),reshape(auc_tmp(1:4,2,:),[],1));
    corr_KO = corr(reshape(auc_tmp(5:8,1,:),[],1),reshape(auc_tmp(5:8,2,:),[],1));
    corr_diff_all(k) = corr_WT - corr_KO;
end
corr_diff_real = corr(reshape(auc_diff(1:4,1,:),[],1),reshape(auc_diff(1:4,2,:),[],1)) - corr(reshape(auc_diff(5:8,1,:),[],1),reshape(auc_diff(5:8,2,:),[],1));
p_val = sum(corr_diff_all>=corr_diff_real)/10000;

%% Fig.3g
figure
for t = 1:20
    for n = 1:4 %WT
        scatter(auc_diff(n,1,t),auc_diff(n,2,t),'k','filled')
        hold on
    end
%     for n = 5:8 %KO
%         scatter(auc_diff(n,1,t),auc_diff(n,2,t),'k','filled')
%         hold on
%     end
end

x2 = reshape(auc_diff(1:4,1,:),1,[]);
y2 = reshape(auc_diff(1:4,2,:),1,[]);
c2 = polyfit(x2,y2,1);
y_est_ko = polyval(c2,x2);
% Add trend line to plot
hold on
plot(x2,y_est_ko,'k-','LineWidth',2)
hold off
xlim([-40 40])
ylim([-5 5])
xline(0)
yline(0)

%% Functions

function ve = calcVEall(y,coeff,n_pcs)
    for pc = 1:n_pcs
        rec = (y' * coeff(:,pc) * coeff(:,pc)')';
        ve(pc) = 100*(1 - norm(y-rec,'fro')^2  / norm(y,'fro')^2);
    end
end

function ve = calcVEind(y,q,cell_idx,n_pcs)
    for sub = 1:8
        org = y(cell_idx{sub},:);
        for pc = 1:n_pcs
            rec = (org' * q{sub}(:,pc) * q{sub}(:,pc)')';
            ve(sub,pc) = 100*(1 - norm(org-rec,'fro')^2  / norm(org,'fro')^2);
        end
    end
end

function dim_sim = calcDimSim(y,q,cell_idx,n_pcs) 
    this_y = y;
    this_q = q;
    cc_all_comb = [];
    for pc1 = 1:n_pcs
        for pc2 = 1:n_pcs
            n_comb = 1;
            for sub1 = 1:8
                sub1_pc = this_y(cell_idx{sub1},:)' * this_q{sub1}(:,pc1);
                for sub2 = sub1+1:8
                    sub2_pc = this_y(cell_idx{sub2},:)' * this_q{sub2}(:,pc2);
                    cc = corrcoef(sub1_pc,sub2_pc);
                    cc_all_comb(pc1,pc2,n_comb) = cc(1,2);
                    n_comb = n_comb+1;
                end
            end
        end
    end
    dim_sim = mean(cc_all_comb,3);
end