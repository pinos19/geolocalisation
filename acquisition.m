clear
close all
clc

fe=16.368e6; % fr�quence d'�chantillonnage
fi=2.046e6; % fr�quence interm�diaire
Tc=1e-3; % p�riode du code C/A
nbr_ech=fe*Tc; % nombre d'�chantillons par p�riode de code
sats=[];
deccode=[];
dopp=[];



for prn=1:32 % num�ro PRN du satellite � �tudier
    ca = 2*gcacode(prn)-1; % code C/A du satellite � �tudier
    % figure 
    % plot(ca)
    sca = sampleca(ca,fe); % code r��chantillonn� � la fr�quence fe

    % figure
    % plot(sca)
    % axis([0 length(sca) -1.2 1.2])
    % xlabel('Echantillon')
    % ylabel('Code C/A')
    % title('Code C/A de PRN 10, �chantillonn� � 16.368 MHz')

    fid = fopen('rec.bin.091','r','l'); % ouverture du fichier en lecture, format little-endian

    N=1;
    K=1;
    scaN=[];
    for n=1:N
        scaN = [scaN sca];
    end

    resultK=0;
    for k=1:K

        data01=fread(fid,N*nbr_ech,'ubit1','l'); % lecture du signal sur 1 p�riode de code
        data=2*data01-1; % conversion de binaire en -1 / 1

        % figure
        % plot(data)
        % axis([0 length(data) -1.2 1.2])
        % xlabel('Echantillon')
        % ylabel('Signal d''antenne')
        % title('Signal brut en sortie d''antenne, �chantillonn� � 16.368 MHz')

    %     fclose(fid); % fermeture du fichier

        ind1 = 1;
        vectdopp =-5000:20:5000;

        for fdoppler=vectdopp
            t=0 : 1/fe : N*1e-3-1/fe;
            porteusecos=cos(2*pi*(fi+fdoppler)*t);
            porteusesin=sin(2*pi*(fi+fdoppler)*t);
            dataport=data.*porteusecos.'+1i*data.*porteusesin.'; 
            prodfft=fft(dataport).*conj(fft(scaN)).';
            r=ifft(prodfft);
            result(:,ind1)=(abs(r)).^2;
            ind1=ind1+1;
        end
        resultK = resultK+result;
    end
    fclose(fid);
    maxr =max(max(result));
    moyr = mean(mean(result));
    if maxr/moyr >20
        figure
        % plot((abs(r)).^2)
        mesh(vectdopp,1:nbr_ech,resultK(1:nbr_ech,:));
        title(['PRN' num2str(prn)]);
        [absc,ord]=find(result==maxr);
        sats=[sats prn];
        deccode=[deccode absc(1)];% le d�calage entre la data et le signal r�plique qu'on a g�n�r�: 
        dopp=[dopp vectdopp(ord(1))]; %le doppler du signal re�u
    end
    

end


clearvars -except sats deccode dopp
save res_acq




