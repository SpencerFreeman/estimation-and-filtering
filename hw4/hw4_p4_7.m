clear;clc;close all


I = eye(3);

T = 3*I;

F = 5*I;

n = 5;

P = zeros(3);
G=P;
for i = 1:n

    j = i - 1;

    P = P + F^j * T * F'^j;

    G = G + F^j;
end

inv(I - F*T*F')*(I - (F*T*F')^n)
P
inv(I - F')*(I - F^n)
G















