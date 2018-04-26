# -*- coding: utf-8 -*-
"""
Created on Sat Apr 14 20:47:31 2018

@author: 脐橙
"""
def taylor1(d,m):
    ux_h=Decimal('0')
    for i in d:
        ux_h=ux_h+i
    return ux_h
def taylor2(h,d,m):
    dux_h=d[1]/h
    for i in range(2,m):
        d[i]=d[i]*Decimal(str(i))/h
        dux_h=dux_h+d[i]
    return dux_h
def f(theta,x):
    z1=n
    z=(z1*(z1+1)/(1-x*x))**0.5-x*sin(2*theta)/(2*(1-x*x))
    z=-1.0/z
    return z
def RK(theta,value_at_theta,emm):
    h=-3.14159/100
    if emm==0:
        nn=100
    else:
        nn=50
    x=value_at_theta
    for i in range(0,nn):
        k1=f(theta,x)
        k2=f(theta+h/2,x+h*k1/2)
        k3=f(theta+h/2,x+h*k2/2)
        k4=f(theta+h,x+h*k3)
        x=x+h*(k1+2*k2+2*k3+k4)/6
        theta=theta+h
    return x
def m_ders(point,value_at_point,der_at_point,m,h):
    dd=[Decimal('0')]
    point=Decimal(str(point))
    if m<=n:
        t=m-2
    else:
        t=n-2
    dd*=(m+2)
    dd[0]=value_at_point
    dd[1]=der_at_point*h
    for k in range(0,t+1):
        k1=k
        dd[k+2]=((Decimal('2')*h*Decimal(str(k1+1)))*point/Decimal(str(k1+2)))*dd[k+1]
        dd[k+2]=dd[k+2]+((h*h*Decimal(str(k1-n))*Decimal(str(k1+n+1)))/(Decimal(str(k1+1))*Decimal(str(k1+2))))*dd[k]
        dd[k+2]=dd[k+2]/(Decimal('1')-point*point)
    return dd
def newton(point,value_at_prepoint,h,error,index,m):
    global roots
    global ders
    d=[Decimal('0')]
    d*=(m+1)
    d=m_ders(roots[index-1],value_at_prepoint,ders[index-1],m,h)
    t=Decimal(str(point))
    point=Decimal(str(point))-taylor1(d,m)/taylor2(h,d,m)
    j=1
    while ((point-t).copy_abs()>error) and (j<=10):
        h=h+point-t
        d=m_ders(roots[index-1],value_at_prepoint,ders[index-1],m,h)
        t=point
        point=point-taylor1(d,m)/taylor2(h,d,m)
        j=j+1
    roots[index]=point
    h=h+point-t
    d=m_ders(roots[index-1],value_at_prepoint,ders[index-1],m,h)
    ders[index]=taylor2(h,d,m)
def find_first_root(xe,uxe,erro,m):
    x=RK(0,xe,1)
    newton(Decimal(str(x)),Decimal(str(uxe)),Decimal(str(x-xe)),10**(-erro),1,m)
from decimal import *
from math import sin
fp=open("input.txt","r")
n=fp.readline().strip('\n')
n=int(n)
prec=fp.readline().strip('\n')
prec=int(prec)
fp.close()
getcontext().prec = prec
r=n*(n+1)
m=30
value_of_zero1=Decimal('1')
value_of_zero2=Decimal('0')
der_of_zero1=Decimal('0')
der_of_zero2=Decimal('1')
if n%2==0:
    for i in range(0,n-1,2):
        value_of_zero1=Decimal(str(-(i+1)))/Decimal(str(i+2))*value_of_zero1
else:
    for i in range(1,n-1):
        value_of_zero1=Decimal(str(-i))/Decimal(str(i+1))*value_of_zero1
        der_of_zero2=Decimal(str(2*i+3))/Decimal(str(i+2))*value_of_zero1-Decimal(str(i+1))/Decimal(str(i+2))*der_of_zero2
roots=[Decimal('0')]
ders=[Decimal('0')]
if n%2==0:
    roots[0]=Decimal('0')
    ders[0]=Decimal('0')
    nn=n//2
    roots*=(nn+1)
    ders*=(nn+1)
    find_first_root(0,float(value_of_zero1),10**(-prec),m)
else:
    nn=n//2+1
    roots*=(nn+1)
    ders*=(nn+1)
    roots[1]=0
    ders[1]=der_of_zero2
for i in range(1,nn):
    roots[i+1]=RK(3.14159/2,float(roots[i]),0)
    roots[i+1]=Decimal(str(roots[i+1]))
    newton(Decimal(str(roots[i+1])),Decimal(str('0')),Decimal(str(roots[i+1]-roots[i])),10**(-prec),i+1,m)

#########
w=[]
tol=Decimal(str(10**(-prec)))
for i in range(1,nn+1):
    p0=Decimal('1')
    p1=Decimal(str(roots[i]))
    for j in range(2,n):
        p2=Decimal(str(2*j-1))*Decimal(str(roots[i]))*p1-Decimal(str(j-1))*p0
        p2=p2/Decimal(str(j))
        p0=p1
        p1=p2
    p2=Decimal(str(n))*p2/Decimal(str(n+1))
    pp1=Decimal('2')*(Decimal('1')-Decimal(str(roots[i]))**Decimal('2'))
    pp2=Decimal(str(n+1))**Decimal('2')*p2**Decimal('2')
    w.append(str((pp1/pp2).quantize(tol)))
fq=open("points_and_coefs.txt","w")
fq1=open("points.txt","w")
fq2=open("coefs.txt","w")
for i in range(1,len(roots)):
    roots[i]=str(roots[i])
    fq.write("x{} = {}".format(i,roots[i][0:17]))
    fq.write('\t\t')
    fq.write("A{} = {}".format(i,w[i-1][0:17]))
    fq.write('\n')
    fq1.write(roots[i])
    fq1.write('\n')
    fq2.write(w[i-1])
    fq2.write('\n')
fq.close()
fq1.close()
fq2.close()