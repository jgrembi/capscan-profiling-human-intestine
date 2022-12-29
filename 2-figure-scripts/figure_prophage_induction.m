clear all;
close all;

load('Mapping_vOTU/viral_contig_prophage_propagate_summary.mat');
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
reads_threshold = 1000;

induction_ind = cell2mat(cellfun(@(x) x', sampletable.prophagestate, 'UniformOutput', false));
induction_ind = induction_ind(:, prophageinds);
induction_ind = induction_ind';

induction_count = cellfun(@(x) length(find(x == 1)), sampletable.prophagestate);
%%

stooldata = induction_count(stoolind);
capsuledata = induction_count(capsuleind);
salivadata = induction_count(salivaind);

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
ylabel('Number of induced prophages');
[h, p] = ttest2(stooldata, capsuledata)
[h, p] = ttest2(stooldata, salivadata)
[h, p] = ttest2(salivadata, capsuledata)

%%
induction_ind(find(induction_ind < 0)) = 0;
induce_prophage_ind = find(max(induction_ind, [], 2) == 1);


%%
% identify induced prophages by subject ID and location
induce_prophage_ind_stool = find(sum(induction_ind(:, stoolind), 2) > 0);
induce_prophage_ind_capsule = find(sum(induction_ind(:, capsuleind), 2) > 0);

figure('color', 'w');
v1 = length(induce_prophage_ind_stool);
v2 = length(induce_prophage_ind_capsule);
v12 = length(intersect(induce_prophage_ind_stool, induce_prophage_ind_capsule));
[~, S] = venn([v1 - v12, v2 - v12, v12]);
pos = S.Position(:, 1);
r = S.Radius;
offset = 0.05 * diff(pos);
text((pos(1) - r(1) + pos(2) - r(2)) / 2 - offset, 0, num2str(v1 - v12));
text((pos(1) + r(1) + pos(2) + r(2)) / 2 - offset, 0, num2str(v2 - v12));
text((pos(1) + r(1) + pos(2) - r(2)) / 2 - offset, 0, num2str(v12));
legend('Stool', 'Capsule');
axis off;
title('Induced prophages');


%%
pH = arrayfun(@(x) x, sampletable.pH);
cov = arrayfun(@(x) x, sampletable.coverage);

figure('color', 'w');
plot(pH(capsuleind), induction_count(capsuleind) + 0.5 * (rand(length(capsuleind), 1) - 0.5), '.');
xlabel('pH of capsule sample');
ylabel('# of induced prophages');

pHdata = pH(capsuleind);
inductiondata = induction_count(capsuleind);
covdata = cov(capsuleind);
inductiondata(isnan(pHdata)) = [];
covdata(isnan(pHdata)) = [];
pHdata(isnan(pHdata)) = [];
% 

[r, p] = corrcoef(pHdata, inductiondata)
% 
text(3.2, 45, strcat('r = ', num2str(r(2), '%.2f')));
text(3.2, 42, strcat('p = ', num2str(p(2), '%.2f')));

axis([3 9 -1 50]);

inductiondata(pHdata < 5.5) = [];
pHdata(pHdata < 5.5) = [];
[r, p] = corrcoef(pHdata, inductiondata)