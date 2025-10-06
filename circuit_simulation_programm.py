# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 14:36:24 2024

@author: sarah
"""

import matplotlib.pyplot as plt
import numpy as np
    
T_voiture = {'Dodge' : {'masse': 1760,'acceleration' : 5.1,'longueur' : 5.28, 'largeur' : 1.95,'hauteur' : 1.35,'Cz' : 0.38,'Cx': 0.3,'mu': 0.1},
      'Toyota' : {'masse': 1615,'acceleration' : 5,'longueur' : 4.51,'largeur' : 1.81,'hauteur' : 1.27,'Cz': 0.29,'Cx': 0.3,'mu': 0.1},
      'Chevrolet' : {'masse': 1498,'acceleration' : 5.3,'longueur' : 4.72,'largeur' : 1.88,'hauteur' : 1.3,'Cz': 0.35,'Cx': 0.3,'mu': 0.1},
      'Mazda' : {'masse': 1385, 'acceleration' : 5.2,'longueur' : 4.3, 'largeur' : 1.75,'hauteur' : 1.23, 'Cz': 0.28,'Cx': 0.3, 'mu': 0.1},
      'Nissan' : {'masse': 1540, 'acceleration' : 5.8,'longueur' : 4.6,'largeur' : 1.79,'hauteur' : 1.36,'Cz': 0.34, 'Cx': 0.3, 'mu': 0.1},
      'Mitsubishi' : {'masse': 1600, 'acceleration' : 5, 'longueur' : 4.51, 'largeur' : 1.81,'hauteur' : 1.48,'Cz': 0.28, 'Cx': 0.3, 'mu': 0.1}}

LST_Voiture = T_voiture.keys()
print(LST_Voiture)
voiture = input("Choisir la marque de la voiture voulu parmi la liste ci dessus en rentrant son nom: ")

dt = 0.001
g = 9.81
r = 6
k = 1/2 * 1.3 * T_voiture[voiture]['largeur'] * T_voiture[voiture]['hauteur'] * T_voiture[voiture]['Cx']
k2 = 1/2 * 1.3 * T_voiture[voiture]['largeur'] * T_voiture[voiture]['longueur'] * T_voiture[voiture]['Cz'] 
    
def rectiligne(v0,alpha,d,nitro):
    
    if nitro == True : 
        acceleration = T_voiture[voiture]['acceleration']*1.3
    else :
        acceleration = T_voiture[voiture]['acceleration']
    it1 = 0
    v = np.arange(0,4500,1,float)
    x = np.arange(0,4500,1,float)

    for i in range(0,4500):
        if i == 0:
            v[i] = v0
            x[i]= 0
        else :
            v[i] = (g*np.sin(alpha) +  acceleration - T_voiture[voiture]['mu']*g*np.cos(alpha) - k/T_voiture[voiture]['masse'] * (v[i-1])**2) * dt + v[i-1]
            K = g*np.sin(alpha) +  acceleration - T_voiture[voiture]['mu']*g*np.cos(alpha)
            x[i] = (K - k/T_voiture[voiture]['masse'] *(v[i-1])**2)* (dt**2)/2 + v[i-1]*dt + x[i-1]
            if x[i] > (d-0.05) and x[i] < (d+0.05):
                it1 = i
    return (it1*dt, v[it1])


def looping(v0,nitro):
    
    if nitro == True : 
        acceleration = T_voiture[voiture]['acceleration']*1.3
    else :
        acceleration = T_voiture[voiture]['acceleration']
    f = 0
    fx = 0
    it2 = 0
    t = np.arange(0,4000)
    w = np.arange(0,4000,1,float)
    theta = np.arange(0,4000,1,float)
    
    for i in range(0,4000):
        if i == 0:
            w[i] = v0 / r
            theta[i]= 0
        else :
            w[i] = (acceleration - g*np.sin(theta[i-1]) - T_voiture[voiture]['mu']*g*np.cos(theta[i-1]) -(T_voiture[voiture]['mu']*r + k/T_voiture[voiture]['masse'])*w[i-1]**2)/r *dt + w[i-1]
            K = (acceleration - g*np.sin(theta[i-1]) - T_voiture[voiture]['mu']*g*np.cos(theta[i-1]) - (T_voiture[voiture]['mu']*r + k/T_voiture[voiture]['masse'])*w[i-1]**2)/r
            theta[i] = K/2 * dt**2 + w[i-1]*dt + theta[i-1]
            if theta[i] > (2*np.pi - 0.05) and theta[i] < (2*np.pi + 0.05):
                it2 = i
                break;
        if theta[i] < np.pi and theta[i] > 0 :
            f = f - T_voiture[voiture]['mu']*(r*w[i-1]**2 + g*np.cos(theta[i-1]))
            fx = fx - k*w[i-1]**2
   
    wf = -np.pi*r*f
    wfx = -np.pi*r*fx
    w = w*r
    
    vmin = np.sqrt((0.5*T_voiture[voiture]['masse']*2*g*r + T_voiture[voiture]['masse']*g*2*r - wf - wfx)/(0.5*T_voiture[voiture]['masse']))
    if v0 < vmin:
        return 0,0
    else:
        plt.plot(t,w)
        plt.xlim(0,it2)
        plt.ylim(0,22)
        plt.xlabel("t en ms")
        plt.ylabel("vitesse (m/s)")
        plt.title("Evolution de la vitesse de la voiture dans le looping en fonction du temps")
        plt.show()
        return (it2*dt,w[it2])


def saut(v0):
    
    it3 = 0
    vx = np.arange(0,1500,1,float)
    x = np.arange(0,1500,1,float)
    vy = np.arange(0,1500,1,float)
    y = np.arange(0,1500,1,float)
    
    for i in range(0,1500):
        if i == 0:
            vx[i] = v0
            x[i] = 0
            vy[i] = 0
            y[i] = 1
            
        else : 
            vx[i] = (-k/T_voiture[voiture]['masse'])*((vx[i-1])**2)*dt+vx[i-1]
            K1 = (-k/T_voiture[voiture]['masse'])*((vx[i-1])**2)
            x[i] =  K1/2*dt**2 + vx[i-1]*dt + x[i-1]
            vy[i] = (-g + (k2/T_voiture[voiture]['masse'])*((vy[i-1])**2))*dt+vy[i-1]
            K2 = (-g + (k2/T_voiture[voiture]['masse'])*((vy[i-1])**2))
            y[i] = K2/2 *dt**2 + vy[i-1]*dt + y[i-1]
            if y[i] > -0.05 and y[i] < 0.05:
                it3 = i
    
    if x[it3] < 9:
        return 0,0
    else:
        plt.figure()
        plt.plot(x,y)
        plt.xlim(0,11)
        plt.ylim(0,4)
        plt.xlabel("x (m)")
        plt.ylabel("y (m)")
        plt.title("Trajectoire de la voiture dans le ravin")
        plt.show()
        return (it3*dt,vx[it3])

a = rectiligne(0,np.radians(3.699),31, True)
b = looping(a[1],False)
c = saut(b[1])
d = rectiligne(c[1],0,10,False)
t_total= a[0]+b[0]+c[0]+d[0]

def resultat():
    print (f"En utilisant la nitro sur la pente, la {voiture} fait un temps total de {t_total} s.")
    print (f"La piste d'élan : La {voiture} fait un temps de {a[0]} s et sort de  avec une vitesse de {a[1]} m.s-1.")
    print (f"Le looping : La {voiture} fait un temps de {b[0]} s et sort de  avec une vitesse de {b[1]} m.s-1.")
    print (f"La ravin : La {voiture} fait un temps de {c[0]} s et sort de  avec une vitesse de {c[1]} m.s-1.")
    print (f"Le sprint final : La {voiture} fait un temps de {d[0]} s et sort de  avec une vitesse de {d[1]} m.s-1.")

def circuit():
    if b[1] == 0:
        return "La voiture ne peut pas passer le looping"
    elif c[1] == 0:
        return "La voiture ne peut pas passer le ravin"
    else:
        if t_total < 8:
            resultat()
            return f"Avec le temps total de {t_total} s, Dom torreto peut battre Owen Shaw."
        else :
            resultat()
            return f"Avec le temps total de {t_total} s, Dom torreto ne peut pas battre Owen Shaw."
           

print(circuit())

