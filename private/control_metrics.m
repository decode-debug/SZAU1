function m = control_metrics(P,y,u,ref,taskMask,distMask)
% Wskazniki jakosci regulacji dla pojedynczego scenariusza:
%   [MSE_zad, MSE_zakl, MSE_laczny, Przeregulowanie, Maks_blad_FD, Min_u, Max_u, Maks_du, Suma_du]
    err = y - ref;
    dir = sign(ref(end) - ref(1));
    overshoot = max(0, max(dir*err(taskMask)));
    du = diff([P.u0; u]);
    m = [mean(err(taskMask).^2), mean(err(distMask).^2), mean(err(taskMask|distMask).^2), ...
         overshoot, max(abs(err(distMask))), min(u), max(u), max(abs(du)), sum(abs(du))];
end
