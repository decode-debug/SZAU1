function weights = ts_memberships(h2,centers)
% Znormalizowane trojkatne funkcje przynaleznosci TS z ramionami na brzegach.
% Jeden wiersz odpowiada jednej wartosci wektora h2; wynik jest rowny trimf.
    x = min(max(h2(:),centers(1)),centers(end));
    count = numel(centers);
    if count==1
        weights = ones(numel(x),1);
        return
    end
    weights = zeros(numel(x),count);
    for i = 1:count
        c = centers(i);
        left = centers(max(i-1,1));
        right = centers(min(i+1,count));
        if i==1
            mu = (right-x)/(right-c);
            mu(x<=c) = 1;
        elseif i==count
            mu = (x-left)/(c-left);
            mu(x>=c) = 1;
        else
            mu = min((x-left)/(c-left),(right-x)/(right-c));
        end
        weights(:,i) = max(mu,0);
    end
    weights = weights./sum(weights,2);
end
