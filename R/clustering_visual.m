TCG1=xlsread('norm_DEG_S.csv');
TCG2=xlsread('norm_DEG_R.csv');

%%% Interpolation
TCG1=pchip([0 6 12 24 48],TCG1,0:1:48);
%%% Interpolation
TCG2=pchip([0 6 12 24 48],TCG2,0:1:48);

%%%%%%%%%%%%%%%%%% hier clustering
%%%%%%%%%%%%%%%  Sensitive genes
corrDist = pdist(TCG1,'euclidean');
clusterTree = linkage(corrDist,'average');
clusters = cluster(clusterTree,'maxclust',9);
for j=1:9
    num1(j)=sum(clusters==j);
end

[a,order]=sort(num1);

figure
figind=9;
for c = order
    subplot(3,3,figind);
    plot(0:1:48,TCG1((clusters == c),:)','Color',[0.8 0.8 0.8],'MarkerSize',2);
    hold on,
    plot([0:1:48],mean(TCG1((clusters == c),:)',2),'k-','LineWidth',3);
    plot([0 6 12 24 48],mean(TCG1((clusters == c),[0 6 12 24 48]+1)',2),'b .','MarkerSize',30);
    axis tight
    figind=figind-1;
end
suptitle('Hierarchical Clustering of DBTRG Sensitive TCGs Profiles');


%%%%%%%%%%%%%%%%%% Resistant genes
corrDist = pdist(TCG2,'euclidean');
clusterTree = linkage(corrDist,'average');

clusters = cluster(clusterTree,'maxclust',9);

for j=1:9
    num2(j)=sum(clusters==j);
end

[a,order]=sort(num2);

figure
figind=9;
for c = order
    subplot(3,3,figind);
    plot(0:1:48,TCG2((clusters == c),:)','Color',[0.8 0.8 0.8],'MarkerSize',2);
    hold on,
    plot([0:1:48],mean(TCG2((clusters == c),:)',2),'k-','LineWidth',3);
    plot([0 6 12 24 48],mean(TCG2((clusters == c),[0 6 12 24 48]+1)',2),'b .','MarkerSize',30);
    axis tight
    figind=figind-1;
end
suptitle('Hierarchical Clustering of LN18 (Resistant) TCGs Profiles');

%%%%%%%%%%%%%%%%%% K-means
rng('default');

TCG1_new=TCG1(sum(TCG1,2)~=0,:);
TCG2_new=TCG2(sum(TCG2,2)~=0,:);



[cidx, ctrs] = kmeans(TCG1_new,9,'dist','sqeuclidean','rep',5,'disp','final');

for j=1:9
    num1(j)=sum(cidx==j);
end

[a,order]=sort(num1);

figure
figind=9;
for c = order
    subplot(3,3,figind);
    plot(0:1:48,TCG1_new((cidx == c),:)','Color',[0.8 0.8 0.8],'MarkerSize',2);
    hold on,
    plot([0:1:48],mean(TCG1_new((cidx == c),:)',2),'k-','LineWidth',3);
    plot([0 6 12 24 48],mean(TCG1_new((cidx == c),[0 6 12 24 48]+1)',2),'b .','MarkerSize',30);
    axis tight
    figind=figind-1;
end
suptitle('K-Means Clustering of Sensitive TCGProfiles');


[cidx, ctrs] = kmeans(TCG2_new,9,'dist','sqeuclidean','rep',5,'disp','final');
for j=1:9
    num2(j)=sum(cidx==j);
end

[a,order]=sort(num2);

figure
figind=9;
for c = order
    subplot(3,3,figind);
    plot(0:1:48,TCG2_new((cidx == c),:)','Color',[0.8 0.8 0.8],'MarkerSize',2);
    hold on,
    plot([0:1:48],mean(TCG2_new((cidx == c),:)',2),'k-','LineWidth',3);
    plot([0 6 12 24 48],mean(TCG2_new((cidx == c),[0 6 12 24 48]+1)',2),'b .','MarkerSize',30);
    axis tight
    figind=figind-1;
end
suptitle('K-Means Clustering of Resistant TCGProfiles');