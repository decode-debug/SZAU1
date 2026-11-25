function h = plant_step(P,h,u,d)
% Jeden krok dyskretny (Ts) obiektu nieliniowego: 5 podkrokow RK4 w poziomach [h1; h2].
    dt = P.Ts/5;
    for j = 1:5
        a = rhs(P,h,u,d);
        b = rhs(P,h+dt*a/2,u,d);
        c = rhs(P,h+dt*b/2,u,d);
        k4 = rhs(P,h+dt*c,u,d);
        h = h + dt*(a+2*b+2*c+k4)/6;
    end
    h = max(h,0.05);
end

function dh = rhs(P,h,u,d)
    f2 = P.a1*sqrt(max(h(1),0));
    f3 = P.a2*sqrt(max(h(2),0));
    cap2 = 3*P.C2*max(h(2),0.05)^2;
    dh = [(u+d-f2)/P.A1; (f2-f3)/cap2];
end
