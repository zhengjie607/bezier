
from numpy import *
import sympy

# 计算bezier曲线在t处的坐标
def q(ctrlPoly, t):
    if len(ctrlPoly)==4:
        return (1.0-t)**3 * ctrlPoly[0] + 3*(1.0-t)**2 * t * ctrlPoly[1] + 3*(1.0-t)* t**2 * ctrlPoly[2] + t**3 * ctrlPoly[3]
    elif len(ctrlPoly)==2:
        return (1.0-t)*ctrlPoly[0]+t*ctrlPoly[1]
    
# 计算bezier曲线在t处的一阶导数
def qprime(ctrlPoly, t):
    if len(ctrlPoly)==4:
        return 3*(1.0-t)**2 * (ctrlPoly[1]-ctrlPoly[0]) + 6*(1.0-t) * t * (ctrlPoly[2]-ctrlPoly[1]) + 3*t**2 * (ctrlPoly[3]-ctrlPoly[2])
    elif len(ctrlPoly)==2:
        return ctrlPoly[1]-ctrlPoly[0]

# 计算bezier曲线在t处的二阶导数
def qprimeprime(ctrlPoly, t):
    if len(ctrlPoly)==4:
        return 6*(1.0-t) * (ctrlPoly[2]-2*ctrlPoly[1]+ctrlPoly[0]) + 6*(t) * (ctrlPoly[3]-2*ctrlPoly[2]+ctrlPoly[1])
    elif len(ctrlPoly)==2:
        return (0,0)

#计算bezier曲线在x点处的y值
def x_to_y(ctrlPoly,x):
    if x<ctrlPoly.min(axis=0)[0] or x>ctrlPoly.max(axis=0)[0]:
        print('x value must between ['+str(ctrlPoly.min(axis=0)[0])+','+str(ctrlPoly.max(axis=0)[0])+']')
        return 
    t=sympy.Symbol('t')
    if len(ctrlPoly)==4:
        C=array([[-1,3,-3,1],[3,-6,3,0],[-3,3,0,0],[1,0,0,0]])
        X=array([[ctrlPoly[0][0]],[ctrlPoly[1][0]],[ctrlPoly[2][0]],[ctrlPoly[3][0]]])
        Y=array([[ctrlPoly[0][1]],[ctrlPoly[1][1]],[ctrlPoly[2][1]],[ctrlPoly[3][1]]])
        CX=C.dot(X)
        tarr=sympy.solve(CX[0][0]*t**3+CX[1][0]*t**2+CX[2][0]*t+CX[3][0]-x,t)
        t=tarr[1].as_real_imag()[0]
        T=array([[t**3,t**2,t,1]])
        return T.dot(C.dot(Y))[0][0]
    elif len(ctrlPoly)==2:
        t=float(sympy.solve((1.0-t)*ctrlPoly[0][0]+t*ctrlPoly[1][0]-x,t)[0])
        return (1-t)*ctrlPoly[0][1]+t*ctrlPoly[1][1]
    
    
def evaluateRange(bezier1,bezier2):
    min1=bezier1.min(axis=0)[0]
    min2=bezier2.min(axis=0)[0]
    max1=bezier1.max(axis=0)[0]
    max2=bezier2.max(axis=0)[0]
    if min1<min2:
        begin=min2
    else:
        begin=min1
    if max1<max2:
        end=max1
    else:
        end=max2
    begin=int(begin)
    end=int(end)+1
    return begin,end
def circulation(begin,end,step,bezier1,bezier2,cof,timer):
    co=cof
    tim=timer
    delta=x_to_y(bezier2,begin/co)-x_to_y(bezier1,begin/co)
    for _ in range(begin,end,step):
        deltay= x_to_y(bezier2,_/co)-x_to_y(bezier1,_/co)
        if delta*deltay==0:
            x=_/co
            y=x_to_y(bezier2,x)
            return x,y
        if delta*deltay<0:
            if tim==0:
                x=(2*_-step)/(2*co)
                y=x_to_y(bezier2,x)
                return x,y
            tim-=1
            co*=10
            begin=(_-step)*10
            end=_*10
            return circulation(begin,end,step,bezier1,bezier2,co,tim)
        delta=deltay
    return None

#计算两条bezier曲线的交点，目前支持一次和三次bezier曲线
def getNode(bezier1,bezier2,count=3,step=1):
    begin,end=evaluateRange(bezier1,bezier2)
    print(begin)
    print(end)
    return circulation(begin,end,step,bezier1,bezier2,1,count)
        
