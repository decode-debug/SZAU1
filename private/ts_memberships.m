function [weights,rawMu] = ts_memberships(h2,centers,sigma)
% Funkcje przynaleznosci TS z plaskimi ramionami na krancach zakresu:
%   bez sigma lub sigma <= 0 -> trojkatne,
%   sigma > 0 -> dzwonowe (Gaussa) o wysokosci 1 i szerokosci sigma [cm].
    x = min(max(h2(:),centers(1)),centers(end));
    count = numel(centers);
    if count==1
        rawMu = ones(numel(x),1);
        weights = rawMu;
        return
    end
    useBell = nargin>=3 && sigma>0;
    rawMu = zeros(numel(x),count);
    for i = 1:count
        c = centers(i);
        if useBell
            mu = exp(-0.5*((x-c)/sigma).^2);
            if i==1
                mu(x<=c) = 1;
            elseif i==count
                mu(x>=c) = 1;
            end
            rawMu(:,i) = mu;
        else
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
            rawMu(:,i) = max(mu,0);
        end
    end
    weights = rawMu./sum(rawMu,2);
end
