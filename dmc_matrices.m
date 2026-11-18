function [M,Mp] = dmc_matrices(s,N,Nu)
% Macierze (4.100) i (4.101). Od probki D model przyjmuje stala rzedna s_D.
    s = s(:);
    D = numel(s);
    history = 1:D-1;
    M = zeros(N,Nu);
    Mp = zeros(N,D-1);
    for p = 1:N
        moves = 1:min(p,Nu);
        M(p,moves) = s(min(p-moves+1,D)).';
        Mp(p,:) = (s(min(p+history,D))-s(history)).';
    end
end
