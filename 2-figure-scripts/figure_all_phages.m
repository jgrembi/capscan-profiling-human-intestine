clear all;
close all;

load('Mapping_vOTU/viral_contig_prophage_coverage_summary.mat');
reads_threshold = 1e6;

sampletable = sortrows(sampletable, 'Abrreviation');

stoolind = intersect(find(cellfun(@(x) strcmpi(x, 'stool'), ...
    sampletable.PutativeLocation)), ...
    find(sampletable.reads_humanRemoved >= reads_threshold));
capsuleind = intersect(find(cellfun(@(x) strcmpi(x, 'capsule'), ...
    sampletable.PutativeLocation)), ...
    find(sampletable.reads_humanRemoved >= reads_threshold));
salivaind = intersect(find(cellfun(@(x) strcmpi(x, 'saliva'), ...
    sampletable.PutativeLocation)), ...
    find(sampletable.reads_humanRemoved >= reads_threshold));

excludeind1 = union(find(cellfun(@(x) strcmpi(x, 'set1'), ...
    sampletable.Set)), ...
    find(cellfun(@(x) strcmpi(x, '1'), sampletable.Set)));

excludeind2 = intersect(find(arrayfun(@(x) x == 3, ...
    sampletable.Subject)), ...
    find(cellfun(@(x) strcmpi(x, 'capsule'), ...
    sampletable.PutativeLocation)));

excludeind3 = unique([find(cellfun(@(x) strcmpi(x, '6'), ...
    sampletable.Set)); ...
    find(cellfun(@(x) strcmpi(x, '7'), sampletable.Set)); ...
    find(cellfun(@(x) strcmpi(x, '8'), sampletable.Set)); ...
    find(cellfun(@(x) strcmpi(x, '9'), sampletable.Set)); ...
    find(cellfun(@(x) strcmpi(x, 'set10'), sampletable.Set)); ...
    find(cellfun(@(x) strcmpi(x, 'set11'), sampletable.Set)); ...
    find(cellfun(@(x) strcmpi(x, 'set12'), sampletable.Set)); ...
    find(cellfun(@(x) strcmpi(x, 'set13'), sampletable.Set))]);

excludeind = unique([excludeind1; excludeind2; excludeind3]);

stoolind = setdiff(stoolind, excludeind);
capsuleind = setdiff(capsuleind, excludeind);
salivaind = setdiff(salivaind, excludeind);

allinds = [stoolind; capsuleind; salivaind];


%%
phagereads = cellfun(@(x) sum(x), sampletable.realreads);
totalreads = sampletable.reads_humanRemoved;

ratios = phagereads ./ totalreads ./ 2 .* 100;


%%

stooldata = ratios(stoolind);
capsuledata = ratios(capsuleind);
salivadata = ratios(salivaind);

colors = lines(7);

figure('color', 'w');
hold on;
for i = 1 : 3
    switch i
        case 1
            data = stooldata;
        case 2
            data = capsuledata;
        case 3
            data = salivadata;
    end
    violinplot(data, i .* ones(size(data)), ...
        'ViolinColor', colors(i, :));
end
hold off;
set(gca, 'XTick', 1 : 3);
set(gca, 'XTickLabel', {'stool', 'capsule', 'saliva'});
ylabel('Percentage of reads mapped as viral');
[h, p] = ttest2(stooldata, capsuledata)
[h, p] = ttest2(stooldata, salivadata)
[h, p] = ttest2(salivadata, capsuledata)


%%


phagedepth = cell2mat(cellfun(@(x) x, sampletable.realdepth, 'UniformOutput', false)');
phagereads = cell2mat(cellfun(@(x) x, sampletable.realreads, 'UniformOutput', false)');

reads_threshold = 100;
depth_threshold = 1;

phagedepth(phagereads <= reads_threshold) = 0;
phagedepth(phagedepth <= depth_threshold) = 0;

present_phage_ind = find(sum(phagedepth, 2) > 0);

% identify phages by subject ID and location
phages_stool = find(sum(phagedepth(:, stoolind), 2) > 0);
phages_capsule = find(sum(phagedepth(:, capsuleind), 2) > 0);
phages_saliva = find(sum(phagedepth(:, salivaind), 2) > 0);


figure('color', 'w');
v1 = length(phages_stool);
v2 = length(phages_capsule);
v3 = length(phages_saliva);
v12 = length(intersect(phages_stool, phages_capsule));
v13 = length(intersect(phages_stool, phages_saliva));
v23 = length(intersect(phages_capsule, phages_saliva));
v123 = length(intersect(intersect(phages_stool, phages_capsule), ...
    phages_saliva));
[~, S] = venn([v1 - v12, v2 - v12, v12]);
pos = S.Position(:, 1);
r = S.Radius;
offset = 0.8 * diff(pos);
text((pos(1) - r(1) + pos(2) - r(2)) / 2 - offset, 0, num2str(v1 - v12));
text((pos(1) + r(1) + pos(2) + r(2)) / 2 - offset, 0, num2str(v2 - v12));
text((pos(1) + r(1) + pos(2) - r(2)) / 2 - offset, 0, num2str(v12));
legend('Stool', 'Capsule');
axis off;
axis equal;
title('All phages');

%%
phagedepth = cell2mat(cellfun(@(x) x, sampletable.realdepth, 'UniformOutput', false)');
phagereads = cell2mat(cellfun(@(x) x, sampletable.realreads, 'UniformOutput', false)');

reads_threshold = 2;
depth_threshold = 0;

phagedepth(phagereads <= reads_threshold) = 0;
phagedepth(phagedepth <= depth_threshold) = 0;

phage_stool = phagedepth(:, stoolind);
phage_capsule = phagedepth(:, capsuleind);
phage_saliva = phagedepth(:, salivaind);
for i = 1 %: 3
    switch i
        case 1
            data1 = phage_stool;
            data2 = phage_capsule;
        case 2
            data1 = phage_stool;
            data2 = phage_saliva;
        case 3
            data1 = phage_capsule;
            data2 = phage_saliva;
    end
    data1_full = median(data1, 2);
    data2_full = median(data2, 2);
    delind = find(data1_full .* data2_full == 0);
    
    ratiodata_full = data1_full ./ data2_full;
    ratiodata = ratiodata_full;
    ratiodata(delind) = [];
%     figure('color', 'w'); histogram(log10(ratiodata));
%     xlabel('log10(ratio)');
%     ylabel('count');
    data1 = data1_full;
    data2 = data2_full;
    data1(delind) = [];
    data2(delind) = [];
    figure('color', 'w');
    hold on;
    plot(log10(data1), log10(data2), '.');
    plot([-4 3], [-4 3], '--k');
    hold off;
    switch i
        case 1
            xlabel('log10(sequencing depth), stool');
            ylabel('log10(sequencing depth), capsule');
        case 2
            xlabel('log10(sequencing depth), stool');
            ylabel('log10(sequencing depth), saliva');
        case 3
            xlabel('log10(sequencing depth), capsule');
            ylabel('log10(sequencing depth), saliva');
    end
    axis([-4 3 -4 3]);
    axis square;
    
    disp(length(data1));
    [r, p] = corrcoef(log10(data1), log10(data2));
    disp(r(2));
    disp(p(2));
end