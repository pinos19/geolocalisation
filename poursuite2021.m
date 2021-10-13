clear
close all
clc
load res_acq
fe=16.368e6; % fréquence d'échantillonnage
fi=2.046e6; % fréquence intermédiaire
Tc=1e-3; % période du code C/A
nbr_ech=fe*Tc; % nombre d'échantillons par période de code
nbrms=4000;%nombre de ms sur lesquelles se fera le traitement
prn=3;% définit le satellitte à utiliser par rapport au vecteur "sats"

ca = 2*gcacode(sats(prn))-1; % code C/A du satellite à étudier
sca = sampleca(ca,fe); % code rééchantillonné à la fréquence fe


scaP =[sca(nbr_ech-deccode(prn)+1:nbr_ech) sca(1:nbr_ech-deccode(prn))];
dsur2 =5;
scaE =[scaP(dsur2+1:nbr_ech) scaP(1:dsur2)]; %code généré en avance
scaL =[scaP(nbr_ech-dsur2+1:nbr_ech) scaP(1:nbr_ech-dsur2)];% code généré en retard
phi(j)
fid =fopen('rec.bin.091','r','l'); %ouverture du fichier en lecture
for j=1:nbrms
    % coder le schéma de poursuite de code
    data01 =fread(fid,nbr_ech,'ubit1','l');
    data = 2*data01 -1;
    
    t=(j-1)*1e-3:1/fe:j*1e-3-1/fe;
    porteusecos=cos(2*pi*(fi+dopp(prn))*t+phi(j));
    porteusesin=sin(2*pi*(fi+dopp(prn))*t+phi(j));
    IE(j)=sum(data.'.*porteusesin.*scaE);
    IP(j)=sum(data.'.*porteusesin.*scaP);
    IL(j)=sum(data.'.*porteusesin.*scaL);
    
    QE(j)=sum(data.'.*porteusecos.*scaE);
    QP(j)=sum(data.'.*porteusecos.*scaP);
    QL(j)=sum(data.'.*porteusecos.*scaL);
    deltaphi = atan(QP(j)/IP(j));
    porteusecos2=cos(2*pi*(fi+dopp(prn))*t+deltaphi);
    porteusesin2=sin(2*pi*(fi+dopp(prn))*t+deltaphi);
    IP2(j)=sum(data.'.*porteusesin2.*scaP);
    QP2(j)=sum(data.'.*porteusecos2.*scaP);
    discELpow=(IE(j)^2+QE(j)^2)-(IL(j)^2+QL(j)^2);
    if discELpow > 0
        % générer les prochains codes (scaE scaP scaL) plus en avance
        scaE=[scaE(2:nbr_ech) scaE(1)];
        scaP=[scaP(2:nbr_ech) scaP(1)];
        scaL=[scaL(2:nbr_ech) scaL(1)];
    elseif discELpow<0
        % générer les prochains codes (scaE scaP scaL) plus en retard
        scaE=[scaE(nbr_ech) scaE(1:nbr_ech-1)];
        scaP=[scaP(nbr_ech) scaP(1:nbr_ech-1)];
        scaL=[scaL(nbr_ech) scaL(1:nbr_ech-1)];
    end
    
end
figure
plot(IP.^2+QP.^2,'.');
hold all
plot(IE.^2+QE.^2,'.k');
plot(IL.^2+QL.^2,'.g');
legend('IP²+QP²','IE²+QE²','IL²+QL²');
figure
plot(IP);
hold all
plot(QP);
legend('IP','QP');
figure
plot(IP2);
hold all
plot(QP2);
legend('IP2','QP2');
fclose(fid)
