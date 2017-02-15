function yvec = spring_damper_model(t, y, param, Fpert, M, tminmax)
    ks=param(1);
    cs=param(2);
    if t>tminmax(1) && t<tminmax(2)
        yvec=[y(2) ; (1/(Fpert-M))*(Fpert-ks*y(1)-cs*y(2)-M)];
    else
        if M>0
            yvec= [-ks*y(1)/cs ; (1/M)*(-ks*y(1)-cs*y(2)-M)];
        else yvec= [-ks*y(1)/cs ; 0];
        end
    end
 
    
