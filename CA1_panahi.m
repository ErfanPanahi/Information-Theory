%% CA1 - Digital Communication
% Erfan Panahi 810198369
clc
clear
fprintf("CA1 - Digital Communication - Erfan Panahi 810198369\n\n");

%% Problem 10. (Test cases)
ts_p1 = [3/4 1/4;1/4 3/4];
G1_p1 = entropy (ts_p1,1);
G1_p2 = entropy (ts_p1,2);
G1_p3 = entropy (ts_p1,3);
fprintf("Problem 10:\n  test case: ts = [3/4 1/4;1/4 3/4]\n  G1 = %f (ex:1.56) , G2 = %f (ex:1.28), G3 = %f (ex:1.136)\n\n",G1_p1,G1_p2,G1_p3);

%% Problem 11. 
chain_p2 = 'AABABB';
avglen1_p2 = average_length(chain_p2,1);
avglen2_p2 = average_length(chain_p2,2);
avglen3_p2 = average_length(chain_p2,3);
fprintf("Problem 11:\n  test case: chain = 'AABABB'\n  n1 = %f , n2 = %f , n3 = %f \n\n",avglen1_p2,avglen2_p2,avglen3_p2);

%% Problem 12.
ts = [0.5 0.5;0.8 0.2];
L = 1e4;
a = 0.5;
b = 0.8;
HX = (b*hb(a) + a*hb(b))/(a+b); 
fprintf("Problem 12:\n  H(X) = %f\n\n",HX);
%P1 = 8/13 , P2 = 5/13
P1 = 8/13; 
P2 = 5/13;
si_vec = rand > cumsum([0,P1,P2]);
for i = 1:L
   ns = [0,ts(si_vec(2)+1,:)];
   si_vec = rand > cumsum(ns);
   if si_vec(2) == 0
       chain(i) = 'A';
   else
       chain(i) = 'B';
   end
end
Gk = zeros(1,10);
avg = zeros(1,10);
eff = zeros(1,10);
for k = 1:10
    Gk(k) = entropy(ts,k);
    avg(k) = average_length(chain,k);
    H_hat = avg(k)/k;
    eff(k) = HX/H_hat;
end
K = 1:10;
figure(1);
subplot(2,1,1);
plot(K,Gk,'b',K,avg,'r');
title('Problem 12','Interpreter','latex');
legend('G_k','average length');
subplot(2,1,2);
plot(K,HX*ones(1,10),'r',K,Gk,'b',K,eff,'g');
legend('H(X)','efficiency');

%% Problem 13
L = 1e3;
HX = zeros(1,3);

% X1
prob1 = [0.01,0.29,0.7];
HX(1) = -sum(prob1.*log2(prob1));
ts1 = ones(3,1) * prob1;
alphabet1 = 1:3;
chain_X1 = char(64+randsrc(L,1,[alphabet1;prob1]));

% X2
prob2 = kron(prob1,prob1);
HX(2) = -sum(prob2.*log2(prob2));
ts2 = ones(9,1) * prob2;
alphabet2 = 1:9;
chain_X2 = char(64+randsrc(L,1,[alphabet2;prob2]));

% X3
prob3 = kron(prob1,prob2);
HX(3) = -sum(prob3.*log2(prob3));
ts3 = ones(27,1) * prob3;
alphabet3 = 1:27;
chain_X3 = char(64+randsrc(L,1,[alphabet3;prob3]));
chain = [chain_X1';chain_X2';chain_X3'];

Gk = zeros(3,10);
avg = zeros(3,10);
eff = zeros(3,10);
for k = 1:10
    Gk(1,k) = entropy(ts1,k);
    Gk(2,k) = entropy(ts2,k);
    Gk(3,k) = entropy(ts3,k);
    for i = 1:3
        avg(i,k) = average_length(chain(i,:),k);
        eff(i,k) = k*HX(i)/avg(i,k);
    end
end

K = 1:10;
figure(2);
subplot(2,1,1);
plot(K,Gk(1,:),'r',K,avg(1,:),'b');
title('Problem 13 (X)','Interpreter','latex');
legend('G_k (X)','average length (X)');
xlim([1,10]);
subplot(2,1,2);
plot(K,HX(1).'*ones(1,10),'b',K,Gk(1,:),'r',K,eff(1,:),'g');
xlim([1,10]);
legend('H(X)','G_k (X)','efficiency (X)');

figure(3);
subplot(2,1,1);
plot(K,Gk(2,:),'r',K,avg(2,:),'b');
title('Problem 13 ($X^2$)','Interpreter','latex');
legend('G_k (X^2)','average length (X^2)');
xlim([1,10]);
subplot(2,1,2);
plot(K,HX(2).'*ones(1,10),'b',K,Gk(2,:),'r',K,eff(2,:),'g');
xlim([1,10]);
legend('H((X^2))','G_k (X^2)','efficiency (X^2)');

figure(4);
subplot(2,1,1);
plot(K,Gk(3,:),'r',K,avg(3,:),'b');
title('Problem 13 ($X^3$)','Interpreter','latex');
legend('G_k (X^3)','average length (X^3)');
xlim([1,10]);
subplot(2,1,2);
plot(K,HX(3).'*ones(1,10),'b',K,Gk(3,:),'r',K,eff(3,:),'g');
xlim([1,10]);
legend('H(X^3)','G_k (X^3)','efficiency (X)');

fprintf("Problem 13:\n  H(X) = %f\n  H(X2) = %f\n  H(X3) = %f\n\n",HX(1),HX(2),HX(3));

% figure(3);
% subplot(1,3,1);
% plot(K(1:5),Gk(:,1:5),K(1:5),avg(:,1:5));
% legend('G_k (X)','G_k (X^2)','G_k (X^3)','average length (X)','average length (X^2)','average length (X^3)');
% subplot(1,3,2);
% plot(K(1:5),HX.'*ones(1,5),K(1:5),eff(:,1:5));
% title('Problem 13 (k = 1:5)','Interpreter','latex');
% legend('H(X)','H(X^2)','H(X^3)','efficiency (X)','efficiency (X^2)','efficiency (X^3)');
% subplot(1,3,3);
% plot(K(1:5),Gk(:,1:5),K(1:5),avg(:,1:5),K(1:5),eff(:,1:5));
% legend('G_k (X)','G_k (X^2)','G_k (X^3)','average length (X)','average length (X^2)','average length (X^3)','efficiency (X)','efficiency (X^2)','efficiency (X^3)');

%% Functions

% Problem 10.
function [Gk] = entropy (ts,k)
    [V,D] = eig(ts.');
    [m,i] = max(sum(D));
    P = V(:,i);
    Pi = P./sum(P);
    H1 = -sum(Pi.*log2(Pi));
    Hi = sum(-ts .* log2(ts),2);
    H_cond = sum(Pi .* Hi);
    Gk = (H1 + k * H_cond)/k;
end

%Problem 11.
function avglen = average_length(chain,k)
    n = length(chain);
    kSymbols = [];
    for i = 1:n-k+1
        kSymbols = [kSymbols chain(i:k+i-1)];
    end
    kSymbols = reshape(kSymbols,[k,n-k+1]);
    [C,ia,ic] = unique(kSymbols.','rows','stable');
    p = zeros(1,length(ia));
    for i = 1:length(ia)
        x = find(ic == i);
        p(i) = length(x)/size(kSymbols,2);
    end
    [dict,avglen] = huffmandict(cellstr(C),p);
end

function H = hb(p)
    H = -p*log2(p) - (1-p)*log2(1-p);
end

