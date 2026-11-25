function [weights,rawMu] = ts_bell_memberships(h2,centers,sigma)
% Dzwonowe (Gaussa) funkcje przynaleznosci o stalej wysokosci 1 i wspolnej szerokosci sigma [cm]
% z plaskimi ramionami dla skrajnych regul; weights sa znormalizowane do sumy 1.
    x = min(max(h2(:),centers(1)),centers(end));
    count = numel(centers);
    if count==1
        rawMu = ones(numel(x),1);
        weights = rawMu;
        return
    end
    rawMu = zeros(numel(x),count);
    for i = 1:count
        c = centers(i);
        mu = exp(-0.5*((x-c)/sigma).^2);
        if i==1
            mu(x<=c) = 1;
        elseif i==count
            mu(x>=c) = 1;
        end
        rawMu(:,i) = mu;
    end
    weights = rawMu./sum(rawMu,2);
end
