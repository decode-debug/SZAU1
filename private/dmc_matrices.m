function [M,Mp] = dmc_matrices(s,N,Nu,D)
% Macierze dynamiczne M i Mp regulatora DMC. Od konca odpowiedzi skokowej przyjmowana jest stala rzedna.
    s = s(:);
    if nargin<4, D = numel(s); end
    nS = numel(s);
    history = 1:D-1;
    M = zeros(N,Nu);
    Mp = zeros(N,D-1);
    for p = 1:N
        moves = 1:min(p,Nu);
        M(p,moves) = s(min(p-moves+1,nS)).';
        Mp(p,:) = (s(min(p+history,nS))-s(history)).';
    end
end
