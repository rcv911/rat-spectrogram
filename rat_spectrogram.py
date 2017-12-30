from numpy import * 
from matplotlib.pyplot import *
from numpy.fft import rfft
from matplotlib import rcParams
from struct import unpack

##parameters

dt=0.00195075 # sample rate
channel=4
dlbite=2

sn=20000  # start time
se=25000  # stop time

sd=abs(se-sn)
start=int((1/dt)*sn)
end=int((1/dt)*se)
duration=int((1/dt)*sd)

lx=10 #figsize
ly=5  #figsize

Hz=50 #we cut 50 Herz from analysis 

##file

f=open('rat22agonist 12mg.wdq','rb')
f.read(5296+channel*dlbite*start)
x=f.read(channel*dlbite*duration)
bite=(channel*dlbite*duration)//2 

l1=unpack(bite*'h',x)

data=l1[2::4]
data=array(data, dtype=float64)
data=(data/2**15)*10

L=int((1/dt)*1)*6  # L=int((1/dt)*1) = 512 => 
# => window length = 1 sec (it because is sample frequency = 512)
S=int((1/dt)*0.1)

k=(len(data)-L+S)//S  # number of windows

Fn=1/(2*dt)
Ln=int((L/2+1)*Hz/Fn)
num=int(L/(2*S))

spec=zeros((k +(2*num), L//2 + 1))
print(L, S, k)

##spectrogram
for el in range(0,k):
    #print(data.shape)
    mas=data[el*S:L+el*S] 
    spec[el+num]=abs(rfft(mas))/(len(mas)/2)


##draw
time=linspace(sn, se, len(data))
"""
time=zeros(len(data))
for el in range(len(data)-1):
    time[el+1]=time[el]+dt
"""
figure('spectrogram', figsize=(lx,ly))
subplot(2,1,1)
xlim(sn, se)
plot(time,data)
grid()
subplot(2,1,2)
imshow(spec.T[:Ln,:], cmap='nipy_spectral' , origin='lower', extent=(time[S],time[-S],0,Hz), aspect=((time[-L//2]-time[L//2])/(Hz-0))/((lx/(0.5*ly)*1/0.775)), )
show()