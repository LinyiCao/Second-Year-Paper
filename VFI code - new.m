clc,clear
rs = 0;
rb = 1;
c = 0.02
F = 0.01
% c = 1e-5;
% F = 1e-5;
sig = 1;

Nk = 101;
Np = 101;
Nd = 101;

k = linspace(0,1,Nk)';      % the third parameter is the number of the discretization
p = linspace(0,1,Np)';
delta = linspace(0,10,Nd)';

Tol_C = 1.0E-12;
Tol_V = 1.0E-12;

V0 = k / 4;
% V0 = k / 2;
C0 = zeros(Np,Nd);
% C0 = ones(Np,Nd);
% C0(:,1) = p;
C1 = zeros(Np,Nd);
for id = 1 : Nd
    C0(:,id) = p;
end


E_C = 1E10;
Vk = zeros(Nk,1);
p_star = zeros(Nk,1);
p_star_index = ones(Nk,1);
d_star = zeros(Nk,1);
d_star_index = ones(Nk,1);
Ns = 0;

while(E_C > Tol_C)
    
    E_Vk = 1E10;
    Nt = 0;
    while (E_Vk > Tol_V) 
        for ik = 2 : Nk
            ki = k(ik);
            V_p_d = zeros(Np,Nd);
            for ip = 1 : Np
                pi = p(ip);
                for id = 1 : Nd
                    di = delta(id);
                    kp = min(C0(ip,id) , ki);
                    ikp = find(k == kp);
                    if (di > 0)
                        V_p_d(ip,id) = (1-(kp/ki)^sig)*pi + (kp/ki)^sig * (exp(-rs*di)*V0(ikp)-c*di) - F;       
%                         V_p_d(ip,id) = (1-(kp/ki)^sig)*pi + (kp/ki)^sig * (exp(-rs*di)*V0(ikp)-(1-exp(-rs*di))*c/rs) - F; 
%                         c*di is the limit of (1-exp(-rs*di))*c/rs as
%                         rs goes to 0
                    else
                        V_p_d(ip,id) = (1 - (kp/ki)^sig) * pi;
                    end
                end
            end
            Vk(ik) = max(max(V_p_d));
            [i,j] = find(V_p_d == Vk(ik));
            p_star_index(ik) = i(1);
            d_star_index(ik) = j(1);
            p_star(ik) = p(p_star_index(ik));
            d_star(ik) = delta(d_star_index(ik));
        end
        E_Vk = sqrt(sum((Vk-V0).^2));
        V0 = Vk;
        Nt = Nt+1
    end

    for ip = 1 : Np
        pi = p(ip);
        for id = 2 : Nd
            di = delta(id);
            Val_k = abs(k - pi - exp(-rb*di)*(k-p_star));
            val_min = min(Val_k);
            index = find( Val_k == val_min);
            C1(ip,id) = k(index(1));
        end
    end
    C1(:,1) = p;
    E_C = sqrt(sum(sum((C1-C0).^2)));
    C0 = C1;
    Ns = Ns + 1
end
%%
fid = fopen('output_C.txt','w+');
str = '\n';
for i = 1 : Nd
    str = ['%15.5f',str];
end
fprintf(fid,str,C1');
fclose(fid);

fid = fopen('output_p_delta_V.txt','w+');
fprintf(fid,'%s\n','        p_star        d_star          Vk');
for i = 1 : Nk
    fprintf(fid,'%15.5f%15.5f%15.5f\n',p_star(i),d_star(i),Vk(i));
end
fclose(fid);

%%
% [dd,pp] = meshgrid(delta,p);
% surf(dd,pp,C1)

figure,plot(k,Vk)

figure,plot(k,d_star)
