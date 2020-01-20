clc, clear all, close all
%% Kvalita a model neurcitosti
%% I) Urceni nominalniho prenosu
K0 = 10;
T10 = 0.3;
T20 = 0.4;
P_0 = tf(K0, [(T10*T20), (T10+T20),1]);

% Urceni parametru pro neurcity model(prenos), ureal funkce vytvori
% neurcity realny parametr systemu, byla vyuzita vlastnost percentage a
% hodnota 1
K = ureal('K', K0, 'Percentage', 1);
T1 = ureal('T1', T10, 'Percentage', 1);
T2 = ureal('T2', T20, 'Percentage', 1);

P = tf(K, [(T1*T2), (T1+T2),1]);

%% II) Volba vahove funkce W1 

% !!!!!!!!!!!!!!!!!! NA TOHLE JE POTREBA SE PODIVAT A NEBO ZEPTAT
% KONIGS!!!!!!!!!!!

K_i = 0.0943211712628501; % Parametry ziskane pomoci toolu Control system tuner
K_p = 0.0708411782367862;

C_0 = tf([K_p, K_i],[1, 0]); % Prenos PI regulatoru s nalezenymi parametry

L_0 = C_0*P_0; % Prenos otevrene smycky
S_0 = 1/(1+L_0); % Citlivostni funkce
T_0 = L_0/(1+L_0); % Komplementarni citlivostni funkce

omega = logspace(-2,4,1000);
omega_1 = logspace(-4,4,1000);
W_1 = minreal(tf([1 4], [8, 0]));

S0_norma = norm(S_0, 'inf');
W1_S0_norma = norm(W_1*S_0,'inf');

% Sirka pasma
sirka_pasma_S0 = bandwidth(S_0); 
sirka_pasma_W1_S0 = bandwidth(W_1*S_0);

freqresp_S0 = abs(squeeze(freqresp(S_0, omega)));
freqresp_W1_S0 = abs(squeeze(freqresp(W_1*S_0, omega)));

figure;
semilogx(omega, freqresp_S0); 
xlabel('\omega[rad/s]');
ylabel('|S_0(j \omega)|');
grid on;

figure;
semilogx(omega_1, freqresp_W1_S0);
xlabel('\omega[rad/s]');
ylabel('|W_1(j \omega) * S_0(j \omega)|');
grid on;

figure;
semilogx(omega, 20*log10(freqresp_W1_S0), 'b');
hold on;
grid on;
semilogx(omega, 20*log10(freqresp_S0), 'r')
xlabel('\omega[rad/s]');
ylabel('[db]');
legend('20log(|W_1(j \omega) * S_0(j \omega)|)', '20log(|S_0(j \omega)|)');

%% III) Volba vahove funkce W2
K_for_W2 = 1.01*K0;
T1_for_W2 = 0.99*T10;
T2_for_W2 = 0.99*T20;

%Nejvetsi mozny prenos (okraj kruznice)
P_for_W2 = tf(K_for_W2, [(T1_for_W2*T2_for_W2), (T1_for_W2+T2_for_W2),1]);
% Vypocet nejvetsi mozne W2
W2 = minreal((P_for_W2/P_0)-1);
% Kontrola normy pro robustni stabilitu
W2_T0_norm = norm(W2*T_0,'inf');

freqresp_T0 = abs(squeeze(freqresp(T_0, omega)));
freqresp_W2_T0 = abs(squeeze(freqresp(W2*T_0, omega)));

figure;
semilogx(omega, freqresp_W2_T0, 'b');
hold on;
semilogx(omega, freqresp_T0, 'r');
xlabel('\omega[rad/s]');
legend('|W_2(j \omega) * T_0(j \omega)|','|T_0(j \omega)|');
grid on;


%% Nyquist pro 10 prenosu
N = 10;
P_10 = usample(P, N);
W = -180:0.001:180;

figure; 
nyquist(P_10, exp(W), 'r');
hold on;
grid on;
nyquist(P_0, 'b');
legend('P(s)', 'P_0 (s)');


%% Navrh regulatoru
%% II) Navrh regulatoru pomoci Control system tuner
K_p = 0.0708411782367862;
K_i = 0.0943211712628501;

C_cst = tf([K_p, K_i],[1, 0]);
L0_cst = C_cst * P;
S0 = 1/(1+L0_cst);
T0 = L0_cst/(1+L0_cst);