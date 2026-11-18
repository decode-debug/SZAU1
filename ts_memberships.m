function weights = ts_memberships(h2,centers)
% Znormalizowane wagi regul TS; jeden wiersz odpowiada jednemu pomiarowi.
    h2 = min(max(h2(:),centers(1)),centers(end));
    count = numel(centers);
    if count==1
        weights = ones(numel(h2),1);
        return
    end
    weights = zeros(numel(h2),count);
    for i = 1:count
        left = centers(max(i-1,1));
        right = centers(min(i+1,count));
        weights(:,i) = trimf(h2,[left,centers(i),right]);
    end
    weights = weights./sum(weights,2);
end
