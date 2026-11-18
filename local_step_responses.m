function S = local_step_responses(P,model,D)
% S(j,i) jest odpowiedzia modelu i na jednostkowy skok F1in w chwili j*Ts.
    t = (0:D)'*P.Ts;
    S = zeros(D,model.count);
    for i = 1:model.count
        sys = ss(model.A(:,:,i),model.B,[0 1],0,'InputDelay',P.tau);
        response = step(c2d(sys,P.Ts,'zoh'),t);
        S(:,i) = response(2:end);
    end
end
