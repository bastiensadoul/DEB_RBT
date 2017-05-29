clear all; close all; clc
t = linspace(0,100,100)';
t_f = 50;

p_M_Q = 1000;
p_M = 500;

if p_M_Q < p_M
p_M_t = min(p_M, p_M_Q +(p_M - p_M_Q)/ t_f * t);
else
p_M_t = max(p_M, p_M_Q +(p_M - p_M_Q)/ t_f * t);
end 

plot(t, p_M_t, 'r')
