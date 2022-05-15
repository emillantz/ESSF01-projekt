%clear
clc
ws = warning('off','all');  % Turn off warning
%data
v_T = 25*1e-3; %v_T 25 mV @ 25C
%specification
Bf = 200;
C1 = 100*1e-9; %100 nF 
C2 = 2.2*1e-6; %2,2 µF 
R_s = 100; %100 Ω
R_L = 1000; %1 kΩ
R1 = 1000; %1 kΩ
R2 = 10000; %10 kΩ
A_ti = 11; %A_t∞
I_AS = (9.4*1e-3); %9,4 mA, max current in AS-stage
%calculations
rprim_pi = 2*(Bf * v_T) / (I_AS / 2); %calculated by hand as 2127,66 Ω
I_CE = 2*1e-3; %2 mA, max current in CE-stage
R_pi2 = (Bf * v_T) / I_CE;
%gm2 = I_CE / v_T;
AB0 = -1060.11; %calculated by hand
%ABtest = -(Bf^2)*(R1/(R1+R2))*(R_L/((R1*R2/(R1+R2))+rprim_pi+R_s));

%poles
p1 = -1/(((rprim_pi * (((R1*R2)/(R1+R2))+R_s)) / (rprim_pi+(((R1*R2)/(R1+R2))+R_s)))*C1);
p2 = -1 / (C2 * R_pi2);
%p1 = -14 536 rad/s, p2 = -181,82 rad/s
LP = abs((1 - AB0)*p1*p2); %2,799e9
w0 = sqrt(LP);
sumP_prim = -sqrt(2)*w0; %-74 822 rad/s
sumP = p1 + p2; %-14 718 rad/s

%Butterworth positioning
N = -(LP)/(sqrt(2)*w0+p1+p2); %-46 572 rad/s

%freq. compensation
s = zpk('s');
ABs = AB0 / ((1 - s/p1)*(1 - s/p2));
A_t = A_ti * -ABs / (1 - ABs);

%4 phantom-zero implementations
A_ti_ph = A_ti * (1 - s/(A_ti*N)) / (1 - s/N); %A_t∞_ph
%inductance in series w/ R1
L1_ph = -R1 / N;
delta_L1 = 1 + ((R_s*(R2+R_L) / (R_s+R2+R_L)) / R1);

ABs_L1 = AB0*(1 - s/N) / ((1 - s/(delta_L1*N)) * ((1 - s/p1)*(1 - s/p2)));
At_L1 = A_ti_ph * -ABs_L1 / (1 - ABs_L1);

%capacitance in parallel w/ R2 [best]
C1_ph = -1 / (N*R2);
delta_C1 = 1 + (R2 /(R_L + (R1*R_s/(R1+R_s))));

ABs_C1 = AB0*(1 - s/N) / ((1 - s/(delta_C1*N)) * ((1 - s/p1)*(1 - s/p2)));
At_C1 = A_ti_ph * -ABs_C1 / (1 - ABs_C1);

%capacitance in parallel w/ R_s
C2_ph = -1 / (N*R_s);
delta_C2 = 1 + (R_s / (R1*(R2+R_L)/(R1+R2+R_L)));

ABs_C2 = AB0*(1 - s/N) / ((1 - s/(delta_C2*N)) * ((1 - s/p1)*(1 - s/p2)));
At_C2 = A_ti_ph * -ABs_C2 / (1 - ABs_C2);

%inductance in series w/ R_L
L2_ph = -R_L / N;
delta_L2 = 1 + ((R2 + (R1*R_s/(R1+R_s)))/R_L);

ABs_L2 = AB0*(1 - s/N) / ((1 - s/(delta_L1*N)) * ((1 - s/p1)*(1 - s/p2)));
At_L2 = A_ti_ph * -ABs_L2 / (1 - ABs_L2);

%plot MATLAB-simulation figures
figure(1);bode((-1).*ABs,'b',(-1).*ABs_C1,'k--',(-1).*ABs_C2,'r--',(-1).*...
    ABs_L1,'g--',(-1).*ABs_L2,'y--');
title('Slingförstärkning, Aβ(s), före och efter fantomnollor'); legend('Aβ(s)','Aβ_Cph1(s)',...
    'Aβ_Cph2(s)','Aβ_Lph1(s)','Aβ_Lph2(s)','Location','Best')

figure(2);bode(A_t,'b',At_C1,'k--',At_C2,'r--',At_L1,'g--',At_L2,'y--');...
hold on;
title('Förstärkning, At, före och efter fantomnollor'); legend('At','At_Cph1','At_Cph2','At_Lph1'...
,'At_Lph2','Location','Best')

figure(3); step(A_t,'b');hold on;step(At_C1,'k--'); step(At_C2,'r--');step(...
At_L1,'g--');step(At_L2,'y--');
title('Stegsvar före och efter fantomnollor');legend('At','At_Cph1',...
'At_Cph2','At_Lph1','At_Lph2');

%% plot lab-values
%sgtitle('Mätresultat för At(s)')
%convert from V_out to A_db
comp_trans = lab_comp(:,2) ./ 0.1; %100 mV V_in
db_comp = mag2db(comp_trans);
nocomp_trans = lab_nocomp(:,2) ./ 0.1;
db_nocomp = mag2db(nocomp_trans);


figure(4)
%plot amplitude
subplot(2,1,1)
semilogx((2*pi) .* lab_comp(:,1), db_comp, 'b');
xlim([630 1*1e7])
title('Uppmätt överföringsfunktion i den kompenserade kretsen'); xlabel('Frekvens (rad/s)'); ylabel('A (dB)')

%plot phase
subplot(2,1,2)
phase_comp = 360 .* 1*1e-6 .* phase(:,1) .* phase(:,2);
semilogx((2*pi) .* phase(:,1), -phase_comp, 'b.')
hold on

poly_c = polyfit(phase(:,1), -phase_comp, 8);
y_c = polyval(poly_c, phase(:,1));
semilogx((2*pi) .* phase(:,1), y_c, '-b')
xlim([0 1*1e6])
title('Uppmätt fasfunktion i den kompenserade kretsen'); xlabel('Frekvens (rad/s)'); ylabel('Fas (deg)')


figure(5)
subplot(2,1,1)
semilogx((2*pi) .* lab_nocomp(:,1), db_nocomp, '-r')
xlim([630 1*1e7])
title('Uppmätt överföringsfunktion i den okompenserade kretsen'); xlabel('Frekvens (rad/s)'); ylabel('A (dB)')

subplot(2,1,2)
phase_nocomp = 360 .* 1*1e-6 .* phase(:,1) .* phase(:,3);
semilogx((2*pi) .* phase(:,1), -phase_nocomp,'r.')
hold on

poly_nc = polyfit(phase(:,1), -phase_nocomp, 8);
y_nc = polyval(poly_nc, phase(:,1));
semilogx((2*pi) .* phase(:,1), y_nc, '-r')
xlim([0 1*1e6])
title('Uppmätt fasfunktion i den okompenserade kretsen'); xlabel('Frekvens (rad/s)'); ylabel('Fas (deg)')

%plot step responses
figure(6)
plot(sr_nocomp(:,1), sr_nocomp(:,2), '-r')
title('Uppmätt stegsvar i den okompenserade kretsen'); xlabel('tid (s)'); ylabel('V_{ut} (V)')

figure(7)
plot(sr_comp(:,1), sr_comp(:,2), '-b')
title('Uppmätt stegsvar i den kompenserade kretsen'); xlabel('tid (s)'); ylabel('V_{ut} (V)')


%% plot LTspice-sim figures
figure(8)
sgtitle('Simuleringsresultat från LTspice (kompenserad krets)')
hold on
%compensated dB
subplot(2,1,1)
semilogx(ltspice_comp(:,1), ltspice_comp(:,2), '-b')
xlim([0 1e6])
title('Överföringsfunktion'); xlabel('Frekvens (rad/s)'); ylabel('A (dB)')
%compensated phase
subplot(2,1,2)
semilogx(ltspice_comp(:,1), ltspice_comp(:,3), '-b')
title('Fasfunktion'); xlabel('Frekvens (rad/s)'); ylabel('Fas (deg)')
%xlim([0 2.6*1e4])

figure(9)
sgtitle('Simuleringsresultat från LTspice (okompenserad krets)')
hold on
%uncompensated dB
subplot(2,1,1)
semilogx(ltspice_nocomp(:,1), ltspice_nocomp(:,2), '-r')
xlim([0 1e6])
title('Överföringsfunktion'); xlabel('Frekvens (rad/s)'); ylabel('A (dB)')
%uncompensated phase
subplot(2,1,2)
semilogx(ltspice_nocomp(:,1), ltspice_nocomp(:,3), '-r')
title('Fasfunktion'); xlabel('Frekvens (rad/s)'); ylabel('Fas (deg)')
%xlim([0 8*1e4])