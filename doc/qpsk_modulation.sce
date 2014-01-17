clc
clear
close
b=[0 1 1 0 1 1 1 0]
 
one=ones(1,100)
zero=-ones(1,100)
 
//generation of bodd signal
bodd=[]
for i=1:2:length(b)
    if(b(i)==1)
        bodd=[bodd one]
        bodd=[bodd one]
    else
        bodd=[bodd zero]
        bodd=[bodd zero]
    end
end
 
//generation of beven signal
beven=[]
for i=2:2:length(b)
    if(b(i)==1)  then
        beven=[beven one]
        beven=[beven one]
    else
        beven=[beven zero]
        beven=[beven zero]
    end
end
beven=[zeros(1,100) beven]
 
//resizing bodd n beven signals
bodd=bodd(1:length(b)*100)
beven=beven(1:length(b)*100)
 
 
//sin and cos carriers
f=1
t=0:0.01:(length(b)-0.01)
coscar=cos(2*%pi*f*t)
sincar=sin(2*%pi*f*t)
 
 
//Idata & Qdata
idata=beven.*coscar
qdata=bodd.*sincar
 
//Qpsk output
qpsk=idata+qdata
 
 
//plot 1 bodd
subplot(7,1,1)
a=gca()
a.data_bounds=[0,-2;length(b),2]
a.grid=[1 -1]
t=0:0.01:(length(b)-0.01)
plot(t,bodd)
title('bodd','color','blue','edgecolor','red')
 
 
//plot 2 beven
subplot(7,1,2)
a=gca()
a.data_bounds=[0,-2;length(b),2]
a.grid=[1 -1]
t=0:0.01:(length(b)-0.01)
plot(t,beven)
title('beven','color','blue','edgecolor','red')
 
 
//coscarrier
subplot(7,1,3)
a=gca()
a.data_bounds=[0,-2;length(b),2]
a.grid=[1 -1]
t=0:0.01:(length(b)-0.01)
plot(t,coscar)
title('coscar','color','blue','edgecolor','red')
 
 
 
//sincarrier
subplot(7,1,4)
a=gca()
a.data_bounds=[0,-2;length(b),2]
a.grid=[1 -1]
t=0:0.01:(length(b)-0.01)
plot(t,sincar)
title('sincar','color','blue','edgecolor','red')
 
 
 
//Idata
subplot(7,1,5)
a=gca()
a.data_bounds=[0,-2;length(b),2]
a.grid=[1 -1]
t=0:0.01:(length(b)-0.01)
plot(t,idata)
title('idata','color','blue','edgecolor','red')
 
 
 
 
//Qdata
subplot(7,1,6)
a=gca()
a.data_bounds=[0,-2;length(b),2]
a.grid=[1 -1]
t=0:0.01:(length(b)-0.01)
plot(t,qdata)
title('qdata','color','blue','edgecolor','red')
 
 
//QPSK
subplot(7,1,7)
a=gca()
a.data_bounds=[0,-2;length(b),2]
a.grid=[1 -1]
t=0:0.01:(length(b)-0.01)
plot(t,qpsk)
title('qpsk','color','blue','edgecolor','red')
end