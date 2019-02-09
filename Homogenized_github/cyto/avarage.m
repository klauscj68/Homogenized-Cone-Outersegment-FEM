clear all
close all
clc

load cyto

size=500;
n_campioni=200;

index=randi(size,1,n_campioni);


ava=sum(drop(index,:),1)/n_campioni;


figure(1)
a=plot(time_downsample,ava);
set(a,'linewidth',2);
aa=xlabel('time [s]');
bb=ylabel('drop [%]');
set(aa,'fontsize',20);
set(bb,'fontsize',20);