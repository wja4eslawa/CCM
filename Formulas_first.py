"""
This file is designed to make it easy to get the numbers from the formulas received in continuous chain model in

--gaussian approach for an ideal chain 

-- real chains or chains in good solution using Douglas-Freed approximation described in (J. Douglas and K. F. Freed, “Renormalization and the two-parameter theory,”
Macromolecules, vol. 17, no.11, pp. 2344–2354, 1984)

It provides results for 

-- size rartio R_g^2_complex/R_g^2_chain for both gaussian and real chains

-- size rartio R_g/R_h for gaussian and chains

-- R_g and R_h for ideal  chains

-- R_g for real chains
"""

import tkinter as tk
from tkinter import ttk
import numpy as np

###Dictioneries that contain lambda functions for parameters in gaaussian case

RG={"rosette":lambda fc,fr:3*(2*fc*(3*fc-2)+fr*(2*fr-1)+8*fc*fr)/(12*(fc+fr)**2),
        "pom-pom":lambda f1,f2,a:3*((f1+f2)*(3*f1+3*f2-2+3*a**2)+3*a*(2*f1*f2+f1+f2)+a**3)/(6*(f1+f2+a)**2),
        "db":lambda a:3*((a+1)*(a**2+5*a+3))/(6*(2+a)**2),
        "snow":lambda f,fc:3*((3*f-2)*(3*fc-2))/(6*f*fc),
        "bb":lambda f,n,a:3*(a**3*n**3+(n-1)*((2*a**2*n**2-a**2*n+3*a*n-2)*f+(a*n**2-2*a*n+3*n-3)*f**2))/(6*(a*n+f*(n-1))**2)
        }
nn=lambda n:np.arange(1,n-1)
RH={"rosette":lambda fc,fr:2*np.sqrt(2/np.pi)*(4*fc**2*(np.sqrt(2)-1)/3+4*fc*(2-np.sqrt(2))/3
                                               +np.pi*fr*(np.sqrt(2)-1)/2+np.pi*fr**2*(2-np.sqrt(2))/2
                                               +fc*fr*(5*np.arctan(1/2)/2+1-np.pi/4))/(fc+fr)**2,
        "pom-pom":lambda f1,f2,a:8*np.sqrt(2/np.pi)/3*((f1+f2)*((a+1)**(3/2)-a**(3/2)-2**(1/2)+1)+a**(3/2)
                                                       +(f1**2+f2**2)*(2**(1/2)-1)+f1*f2*((a+2)**(3/2)-2*(a+1)**(3/2)+(a)**(3/2)))/(f1+f2+a)**2,
        "db":lambda a:2*np.sqrt(2/np.pi)*(np.pi+4/3*a**(3/2)+(a+(a**(3/2)/2+2*a**(1/2))*np.arcsin(np.sqrt(a/(4+a)))-a**(3/2)*np.pi/4)*2
                                          +((np.arctan((4*a+1+np.sqrt(4*a+2))/(2*np.sqrt(a)))-np.arctan((4*a+1-np.sqrt(4*a+2))/(2*np.sqrt(a))))*np.sqrt(4*a+2)-
                                                2*(np.arcsin(1/(2*np.sqrt(4*a+1)))+np.arctan(1/(2*np.sqrt(a))))))/((2+a)**2),

        
        "snow":lambda f,fs:2*np.sqrt(2/np.pi)*(-(4/141*(3*np.sqrt(3)-np.sqrt(2)-4))*(-110+149*f-9*np.sqrt(3)*fs*f+15*np.sqrt(3)*f+21*np.sqrt(3)*f**2+110*fs-30*np.sqrt(3)+38*np.sqrt(2)
                                                                                     +30*np.sqrt(3)*fs+3*np.sqrt(3)*np.sqrt(2)*f**2-18*np.sqrt(3)*np.sqrt(2)*f-36*np.sqrt(3)*np.sqrt(2)*fs
                                                                                     +36*np.sqrt(3)*np.sqrt(2)+39*np.sqrt(3)*np.sqrt(2)*f*fs+47*fs*f**2-17*f**2-127*fs*f+11*f**2*np.sqrt(2)
                                                                                     -19*f*np.sqrt(2)-38*fs*np.sqrt(2)+49*fs*f*np.sqrt(2))*fs)/(f*fs)**2,
        "bb":lambda f,n,a:2*np.sqrt(2/np.pi)*((4/3*(n-1))*f+(4/3)*a**(3/2)*n+(1/2*(n-1))*f*(f-1)*(-8/3+(8/3)*np.sqrt(2))+2*f*(n-1)*(-4/3-(4/3)*a**(3/2)+(4/3)*(1+a)**(3/2))
                                              +(n-1)*(-(8/3)*a**(3/2)+(8/3)*np.sqrt(2)*a**(3/2))+2*f*(-4/3*(1+a)**(3/2)*n+4/3*(1+a)**(3/2)
                                        +4/3*a**(3/2)*n-4/3*a**(3/2)+4/3*(a*(n-2)+a+1)**(3/2)*n-4/3*(a*(n-2)+a+1)**(3/2)-4/3*(a*(n-2)+a)**(3/2)*n+4/3*(a*(n-2)+a)**(3/2))+
                                        np.sum((4/3)*f**2*(a*nn(n)+2)**(3/2)*(n-nn(n)-1)-(8/3)*f*(a*nn(n)+a+1)**(3/2)*nn(n)+(8/3)*(a*(nn(n)+1))**(3/2)*(f*nn(n)-n+nn(n)+1)
                                               +(4/3)*(a*(nn(n)+2))**(3/2)*(n-nn(n)-1)-(8/3)*(a*nn(n)+1)**(3/2)*f*(f*n-f*nn(n)-f-nn(n))
                                               +(4/3)*(a*nn(n))**(3/2)*(f**2*n-f**2*nn(n)-f**2-2*f*nn(n)+n-nn(n)-1))
                                        )/(a*n+f*(n-1))**2
        }
ratio={"rosette":lambda fc,fr:(2*fc*(3*fc-2)+fr*(2*fr-1)+8*fc*fr)/(2*(fc+fr)**3),
        "pom-pom":lambda f1,f2,a:((f1+f2)*(3*f1+3*f2-2+3*a**2)+3*a*(2*f1*f2+f1+f2)+a**3)/((f1+f2+a)**3),
        "db":lambda a:((a+1)*(a**2+5*a+3))/((a+2)**3),
        "snow":lambda f,fc:((3*f-2)*(3*fc-2))/((f*fc)**2),
        "bb":lambda f,n,a:(a**3*n**3+(n-1)*((2*a**2*n**2-a**2*n+3*a*n-2)*f+(a*n**2-2*a*n+3*n-3)*f**2))/((a*n+f*(n-1))**3)
        }

"""For real polymers gyration radius is given as R_g^2=R_g^2_gauss(1+C*u) where u is a coupling constant (see related papers) and 
C is the contribtion from the diagrams. It is provided in the dictionary with some support functions"""

### Support functions for C

C_rosette =lambda fc,fr: (-(2/71295*(101*np.sqrt(2)-138))*fc*(30217*fc*np.sqrt(2)+19012*fc**2-23450*np.sqrt(2)-12146*fc+2380)+(1/6)*np.pi*fr*(14*fr**2*np.sqrt(2)
    -26*fr*np.sqrt(2)-24*fr**2+12*np.sqrt(2)+36*fr-15)-(1/24)*np.arcsin((1/5)*np.sqrt(5))*fc*fr*(791*fr+656*fc-939)+(1/48)*np.pi*fc*fr*(112*np.sqrt(2)*(fr-1)
    +464*fc+599*fr-843)-(808/15)*fc*fr*np.sqrt(2)*(-1+fc)+(1/420)*fc*fr*(-17881*fr+20944*fc-13339))/(6*fc**2+8*fc*fr+2*fr**2-4*fc-fr)

C_pp=lambda f1,f2,a:(2/315)*(-18312*f1*f2+8955*f2+9156*a*f1**2-4242*np.sqrt(2)*f1**3-7680*np.sqrt(2)*f2-4242*np.sqrt(2)*f2**3
+11922*np.sqrt(2)*f2**2+4242*a*np.sqrt(2)*f2+4242*a*np.sqrt(2)*f1+8955*f1+9156*a*f2**2-9156*f2*a+9156*f1*f2**2
+9156*f1**2*f2-9156*f1*a-4242*a*np.sqrt(2)*f1**2-4242*np.sqrt(2)*f1*f2**2+8484*np.sqrt(2)*f1*f2-14751*f1**2-14751*f2**2
+5796*f1**3+5796*f2**3-4242*np.sqrt(2)*f1**2*f2-4242*a*np.sqrt(2)*f2**2-201*a**5*f1/(a+1)**(3/2)-201*a**5*f2/(a+1)**(3/2)
-201*a**(7/2)*f1*f2+5460*f1*f2*a**(3/2)+4872*f1*f2*a**(5/2)-2100*a**(3/2)*f1*f2**2-2100*a**(3/2)*f1**2*f2-1260*a**(3/2)*f1**2*f2**2
-8484*a*f1**2/(a+1)**(3/2)-8484*a*f2**2/(a+1)**(3/2)-8484*a**3*f1**2/(a+1)**(3/2)-8484*a**3*f2**2/(a+1)**(3/2)-2436*a**4*f1**2/(a+1)**(3/2)
-2436*a**4*f2**2/(a+1)**(3/2)+10086*a**2*f1/(a+1)**(3/2)+10086*a**2*f2/(a+1)**(3/2)+6474*a**3*f1/(a+1)**(3/2)+6474*a**3*f2/(a+1)**(3/2)
-16062*f1*f2/(a+1)**(3/2)+7479*a*f1/(a+1)**(3/2)+7479*a*f2/(a+1)**(3/2)-12096*a**2*f1**2/(a+1)**(3/2)-12096*a**2*f2**2/(a+1)**(3/2)
+1431*a**4*f1/(a+1)**(3/2)+1431*a**4*f2/(a+1)**(3/2)-33936*f1*f2**2/(2+a)**(3/2)-33936*f1**2*f2/(2+a)**(3/2)
+8232*f1**2*f2/(a+1)**(3/2)+8232*f1*f2**2/(a+1)**(3/2)+3360*a**(3/2)*f1**2+3360*a**(3/2)*f2**2+61440*f1*f2/(2+a)**(3/2)
-201*a**(7/2)+85728*f1*f2*a**2/(2+a)**(3/2)+28308*f1*f2**2*a/(a+1)**(3/2)+19068*f1**2*f2*a**3/(a+1)**(3/2)+19068*f1*f2**2*a**3/(a+1)**(3/2)
+402*a**5*f1*f2/(a+1)**(3/2)-201*a**5*f1*f2/(2+a)**(3/2)+119664*f1*f2*a/(2+a)**(3/2)+4872*f1**2*f2*a**4/(a+1)**(3/2)+28308*f1**2*f2*a/(a+1)**(3/2)
+4872*f1*f2**2*a**4/(a+1)**(3/2)-2436*f1**2*f2*a**(5/2)-2436*f1*f2**2*a**(5/2)-67044*f1*f2*a**2/(a+1)**(3/2)-36636*f1*f2*a**3/(a+1)**(3/2)
-54606*f1*f2*a/(a+1)**(3/2)-49644*f1**2*f2*a**2/(2+a)**(3/2)-16968*f1**2*f2*a**3/(2+a)**(3/2)-2436*f1**2*f2*a**4/(2+a)**(3/2)
-67872*f1*f2**2*a/(2+a)**(3/2)+34272*f1*f2**2*a**2/(a+1)**(3/2)-49644*f1*f2**2*a**2/(2+a)**(3/2)-16968*f1*f2**2*a**3/(2+a)**(3/2)
-2436*f1*f2**2*a**4/(2+a)**(3/2)-67872*f1**2*f2*a/(2+a)**(3/2)+34272*f1**2*f2*a**2/(a+1)**(3/2)+2520*a**2*f1**2*f2**2/(a+1)**(3/2)
-2520*a**2*f1**2*f2**2/(2+a)**(3/2)+2520*a**3*f1**2*f2**2/(a+1)**(3/2)-1260*a**3*f1**2*f2**2/(2+a)**(3/2)+27156*f1*f2*a**3/(2+a)**(3/2)
+2862*f1*f2*a**4/(2+a)**(3/2)-7734*f1*f2*a**4/(a+1)**(3/2)+11922*np.sqrt(2)*f1**2-7680*np.sqrt(2)*f1+2235*f1/(a+1)**(3/2)+201*a**(7/2)*f2+201*a**(7/2)*f1
+2436*a**(5/2)*f2**2+2436*a**(5/2)*f1**2-3360*a**(3/2)*f1-3360*a**(3/2)*f2-2436*a**(5/2)*f1-2436*a**(5/2)*f2-2436*f2**2/(a+1)**(3/2)+2235*f2/(a+1)**(3/2)
-2436*f1**2/(a+1)**(3/2))/(a**3+3*a**2*f1+3*a**2*f2+6*a*f1*f2+3*a*f1+3*a*f2+3*f1**2+6*f1*f2+3*f2**2-2*f1-2*f2)

C_db=lambda a:(-(1/144)*(np.arctan((1/2)*(4*a+1+np.sqrt(4*a+2))/np.sqrt(a))-np.arctan((1/2)*(4*a+1-np.sqrt(4*a+2))/np.sqrt(a)))*(32*a**4-1088*a**3-2128*a**2-1684*a-425)/
((2*a+1)*np.sqrt(4*a+2))+(1/72)*np.arcsin(1/np.sqrt(4*a+1))*(288*a**2-964*a-85)-(1/630)*(1072*a**6+8308*a**5+582*a**4+42588*a**3+36505*a**2+7630*a+420)/
(np.sqrt(a)*(4*a+1)*(2*a+1))-(1/72)*np.pi*(192*a**2-332*a+47)-(1/3)*np.arctan(1/(2*np.sqrt(a)))-(4/9)*a**(7/2)/((2*a+1)*(4*a+1)))/(2+a)**2

C_snow=lambda f,fs:(1/105)*(-27840+62573*np.sqrt(2)*f**2*fs**2-57183*np.sqrt(2)*f*fs**2+68495*np.sqrt(2)*f**3*fs
-12400*f-130155*np.sqrt(2)*f**2*fs+106512*np.sqrt(3)*f*fs**2-141008*np.sqrt(3)*f**2*fs**2
+74413*np.sqrt(2)*f*fs+34720*np.sqrt(3)*f**4*fs-133662*np.sqrt(3)*f*fs-11200*np.sqrt(3)*f**4*fs**2
-17535*np.sqrt(2)*f**4*fs-25753*np.sqrt(2)*f**3*fs**2+291178*np.sqrt(3)*f**2*fs-180096*np.sqrt(3)*f**3*fs
+5565*np.sqrt(2)*f**4*fs**2+69216*np.sqrt(3)*f**3*fs**2+11970*np.sqrt(2)*f**4-23520*np.sqrt(3)*f**4-45570*np.sqrt(2)*f**3
+110880*np.sqrt(3)*f**3+75530*np.sqrt(2)*f**2-150170*np.sqrt(3)*f**2-22350*np.sqrt(2)*f+27150*np.sqrt(3)*f
+11970*np.sqrt(2)*fs**2-23520*np.sqrt(3)*fs**2+12730*np.sqrt(2)*fs-12140*np.sqrt(3)*fs+35660*np.sqrt(3)
-24700*np.sqrt(2)-34545*f**4*fs+211113*f**3*fs+11235*f**4*fs**2-82299*f**3*fs**2+121661*fs*f-314217*f**2*fs
-101801*f*fs**2+153419*f**2*fs**2+149340*f**2+23310*fs**2+23310*f**4-124950*f**3+4530*fs)/(f*(3*fs-2)*(3*f-2))

## Main dictionary for C paramether
C={"rosette":lambda fc,fr:C_rosette(fc,fr)/(fc+fr)**(1/2),
        "pom-pom":lambda f1,f2,a:C_pp(f1,f2,a)/(f1+f2+a)**(1/2),
        "db":lambda a:C_db(a)/(RG["db"](a)/3*(2+a)**0.5),
        "snow":lambda f,fs:C_snow(f,fs)/(f*fs)**0.5,
        "bb":1#lambda f,n,a:C_bottbrush(f,n,a)/(RG["bb"](f,n,a)/3*(a*n+f*(n-1))**(5/2))
        }

### Dictionary with citations and comments
citations={"rosette": 'Rosette polymer\n\n Gyration radius and size ratio of an idea chain:\n\n V. Blavatska and R. Metzler,\nJ. Phys. A: Math. Theor., vol. 48, no. 13, p. 135001, 2015.\n'+
           'Hydrodynamic radius and size ratio of an idea chain:\n\n K. Haydukivska, V. Blavatska, and J.  Paturej,\n Sci. Rep., vol. 10, no. 1, p. 14127, 2020.\n'+
                'Gyration radius and size ratio of a real chain:\n\nK. Haydukivska, V. Blavatska, J. S. Klos, and J. Paturej,\n PRE, vol. 105, no. 3, p. 034502, 2022.\n'+
                "\n* Other citations regarding the architecture\n can be found inside this papers\n"+
                '** Results for radii have to be multiplyed\n by a corresponding Khun length in order\n to be compaired with experiments or simulation\n',

           "pom-pom":'Pom-pom polymer\nGyration radius and size ratio:\n\n K. Haydukivska, O. Kalyuzhnyi, V. Blavatska, and J. Ilnytskyi,\n   J. Mol. Liq., vol. 328, p. 115456, 2021.\n '+
                                                                         'K. Haydukivska, O. Kalyuzhnyi, V. Blavatska, and J. Ilnytskyi,\n  Condens. Matter Phys., vol. 25, no. 2, p. 23302, 2022.\n'+
                                                                         'K. Haydukivska and Blavatska,\n  Condens. Matter Phys., vol. 26, no. 2, p. 23301, 2023.\n'+
                "\n* Other citations regarding the architecture\n can be found inside this papers\n"+
                '** Results for radii have to be multiplyed\n by a corresponding Khun length in order\n to be compaired with experiments or simulation\n',
        
           "db":'Dumbbell polymer\n\nGyration radius and size ratio:\n\nK. Haydukivska, V. Blavatska, and J. Paturej, \n PRE, vol. 108, no. 3, p. 034502, 2023.\n'+
                "\n* Other citations regarding the architecture\n can be found inside this papers\n"+
                '** Results for radii have to be multiplyed\n by a corresponding Khun length in order\n to be compaired with experiments or simulation\n'+
                "*** For hydrodynamic radius and\n corresponding ratio please cite:\n K. Haydukivska and V. Blavatska,\n J. Phys. A: Math. Theor., vol. 55, no. 14, p. 145001, 2022.\n as rusult is derived from diagrams calculated there",
        
           "snow":'Snowflake polymer\n\nGyration radius and size ratio:\n\nK. Haydukivska, V. Blavatska, and J. Paturej,\n J. Mol. Liq., vol. 392, p. 123430, 2023.\n'+
                "\n* Other citations regarding the architecture\n can be found inside this papers\n"+
                '** Results for radii have to be multiplyed\n by a corresponding Khun length in order\n to be compaired with experiments or simulation\n'+
                "*** For hydrodynamic radius and \ncorresponding ratio please cite:\n K. Haydukivska and V. Blavatska,\n J. Phys. A: Math. Theor., vol. 55, no. 14, p. 145001, 2022.\n as rusult is derived from diagrams calculated there",
        
            "bb":'Bottlebrush polymer\n\nResults for an ideal chains :\n K. Haydukivska and V. Blavatska,\n J. Phys. A: Math. Theor., vol. 55, no. 14, p. 145001, 2022.\n'+
                "\n* Other citations regarding the architecture\n can be found inside this papers\n"+
                "\n Data for excluded volume interactions will be added\n on a later stage after publication of the results\n\n"+
                '** Results for radii have to be multiplyed\n by a corresponding Khun length in order\n to be compaired with experiments or simulation\n'
        }

### Functions that receive the parameters from input and provide an output

def parameters(F):
    ## input tuple
    arch_par=[i.get() for i in F]
    ##message for errors string
    message=""
    for i in range(len(arch_par)):
        try:
            arch_par[i]=float(arch_par[i])
        except:
            arch_par=[]
            message="Input has to be a possitive number or 0"
            return arch_par,message
    if any([i<0 for i in arch_par]): message="Input has to be a possitive number or 0"
    if sum(arch_par)==0:message="Atleast one parameter must be non-zero"
    return arch_par,message

def gyration_raduis(topology,out,cite,*F):
    global RG
    arch_par,message=parameters(F)
    if message!="":
        cite.set(message)
    else:
        Rg=RG[topology](*arch_par[:-1])
        out.set(round((Rg*arch_par[-1]**1)**0.5,3))
        cc=citations[topology]
        cite.set(cc)
        
    
    return Rg

def gyration_raduis_DF(topology,out,cite,*F):
    global RG,C
    arch_par,message=parameters(F)
    if message!="":
        cite.set(message)
    else:
        if topology =="bb":
            out.set(np.NaN)
            cc=citations[topology]
            cite.set(cc)
        else:
            a=(3/32)*C[topology](*arch_par[:-1])/RG[topology](*arch_par[:-1])
            Rg=(1+a)*RG[topology](*arch_par[:-1])*arch_par[-1]**1.17518
            out.set(round((Rg)**0.5,3))
            cc=citations[topology]
            cite.set(cc)
          
def hydrodynamic_raduis(topology,out,cite,*F):
    global RH
    arch_par,message=parameters(F)
    if message!="":
        cite.set(message)
    else:
        Rh=(RH[topology](*arch_par[:-1]))**(-1)
        out.set(round((Rh*arch_par[-1]**0.5),3))
        cc=citations[topology]
        cite.set(cc)
        return Rh
          
def size_ratio_g_c(topology,out,cite,*F):
    global ratio
    arch_par,message=parameters(F)
    if message!="":
        cite.set(message)
    else:
        gc=ratio[topology](*arch_par)
        out.set(round(gc,3))
        cc=citations[topology]
        cite.set(cc)

def size_ratio_rho(topology,out,cite,*F):
    global RG,RH
    arch_par,message=parameters(F)
    if message!="":
        cite.set(message)
    else:
        Rg=RG[topology](*arch_par)**0.5
        Rh=(RH[topology](*arch_par))**(-1)
        out.set(round(Rg/Rh,3))
        cc=citations[topology]
        cite.set(cc)

def size_ratio_g_c_DF(topology,out,cite,*F):
    global ratio,C,RG
    arch_par,message=parameters(F)
    if message!="":
        cite.set(message)
    else:
        if topology == "bb":
            out.set(np.NaN)
            cc=citations[topology]
            cite.set(cc)
        else:
            a=(3/32)*C[topology](*arch_par)+1/4
            gc=ratio[topology](*arch_par)*(1-a)/(1-73/560)
            out.set(round(gc,3))
            cc=cite.get()
            cc=citations[topology]
            cite.set(cc)

### Code for the interface

window=tk.Tk()
window.title("SChoPS")
window.geometry('1200x1000')
window.resizable(0, 0)
title_frame=ttk.Frame(window)
title=tk.Label(title_frame,
                text="Size characteristics of Polymer structures\n as predicted by continuous chain model and Douglas-Freed approximation",
                font="Calibri 14 bold",anchor=tk.CENTER)
title_frame.place(relx = 0.5, 
                   rely = 0.04,
                   anchor = 'center')

title.pack(fill=tk.X)


##Comments and references
comments_frame=ttk.Frame(window)
comments_frame.place(relx=0.66,rely=0.55,relwidth=0.33,relheight=0.4)

#variables 
cite=tk.StringVar(value="")
comments_gc_gauss=tk.StringVar()
comments_gc_DF=tk.StringVar()
comments_rho_gauss=tk.StringVar()
comments_Rg_gauss=tk.StringVar()
comments_Rh_gauss=tk.StringVar()
comments_Rg_DF=tk.StringVar()

tk.Label(comments_frame,text="References and comments",font="12").pack()
tk.Label(comments_frame,textvariable=cite,font="Arial 10").pack()

#Button labels
name_rg="Gyration radius \n of an ideal chain"
name_rh="Hydrodynamic radius \n of an ideal chain"
name_rg_DF="Gyration radius \n of a real chain"
name_gc="Size ratio g \n of an ideal chain"
name_rho="Size ratio rho \n of an ideal chain"
name_gc_DF="Size ratio g \n of a real chain"

"""A block for displaying input and results for rosrette polymers"""
#frame
rosette_frame=ttk.Frame(window)
rosette_frame.place(relx=0,rely=0.08,relwidth=0.33,relheight=0.45)

rosette_frame_input=ttk.Frame(rosette_frame)
#variables
rosette_fc=tk.StringVar(value="0")
rosette_fr=tk.StringVar(value="0")
rosette_N=tk.StringVar(value="0")

rosette_gc_gauss=tk.StringVar(rosette_frame)
rosette_gc_DF=tk.StringVar(rosette_frame)
rosette_rho_gauss=tk.StringVar(rosette_frame)
rosette_Rg_gauss=tk.StringVar(rosette_frame)
rosette_Rh_gauss=tk.StringVar(rosette_frame)
rosette_Rg_DF=tk.StringVar(rosette_frame)

#Labels and inputs
rosette=tk.Label(rosette_frame,text="Rosette polymer",font="12")

rosette_fc_text=tk.Label(rosette_frame_input,text="Number of chains").pack()
rosette_fc_entry=ttk.Entry(rosette_frame_input,textvariable=rosette_fc,width=8).pack()
rosette_fr_text=tk.Label(rosette_frame_input,text="Number of rings").pack()
rosette_fr_entry=ttk.Entry(rosette_frame_input,textvariable=rosette_fr,width=8).pack()
rosette_N_text=tk.Label(rosette_frame_input,text="Polymerization degree \n of an arm").pack()
rosette_N_entry=ttk.Entry(rosette_frame_input,textvariable=rosette_N,width=8).pack()

#Buttons and outputs
rosette_button_gc_gauss=ttk.Button(rosette_frame,command=lambda:size_ratio_g_c("rosette",rosette_gc_gauss,cite,rosette_fc,rosette_fr),text=name_gc)
rosette_gc_gauss_out=tk.Label(rosette_frame,textvariable=rosette_gc_gauss,font="Calibri 12 bold")

rosette_button_gc_DF=ttk.Button(rosette_frame,command=lambda:size_ratio_rho("rosette",rosette_rho_gauss,cite,rosette_fc,rosette_fr),text=name_rho)
rosette_gc_out_DF=tk.Label(rosette_frame,textvariable=rosette_rho_gauss,font="Calibri 12 bold")

rosette_button_rho_gauss=ttk.Button(rosette_frame,command=lambda:size_ratio_g_c_DF("rosette",rosette_gc_DF,cite,rosette_fc,rosette_fr),text=name_gc_DF)
rosette_rho_gauss_out=tk.Label(rosette_frame,textvariable=rosette_gc_DF,font="Calibri 12 bold")

rosette_button_Rg_gauss=ttk.Button(rosette_frame,command=lambda:gyration_raduis("rosette",rosette_Rg_gauss,cite,rosette_fc,rosette_fr,rosette_N),text=name_rg)
rosette_Rg_gauss_out=tk.Label(rosette_frame,textvariable=rosette_Rg_gauss,font="Calibri 12 bold")

rosette_button_Rh_gauss=ttk.Button(rosette_frame,command=lambda:hydrodynamic_raduis("rosette",rosette_Rh_gauss,cite,rosette_fc,rosette_fr,rosette_N),text=name_rh)
rosette_Rh_gauss_out=tk.Label(rosette_frame,textvariable=rosette_Rh_gauss,font="Calibri 12 bold")

rosette_button_Rg_DF=ttk.Button(rosette_frame,command=lambda:gyration_raduis_DF("rosette",rosette_Rg_DF,cite,rosette_fc,rosette_fr,rosette_N),text=name_rg_DF)
rosette_Rg_out_DF=tk.Label(rosette_frame,textvariable=rosette_Rg_DF,font="Calibri 12 bold")

#picture

rosette_c=tk.Canvas(rosette_frame,height=200,width=200)
rosette_pic=tk.PhotoImage(file="pics/rosette1.png",master = rosette_c)
rosette_c.create_image(0,0,image=rosette_pic,anchor="nw")

#grid
rosette_frame.columnconfigure((0,1,2),weight=1)
rosette_frame.rowconfigure((0,),weight=1)
rosette_frame.rowconfigure((1,),weight=1)
rosette_frame.rowconfigure((2,3),weight=3)


rosette.grid(row=0,column=0,columnspan=3)

rosette_c.grid(row=1,column=0,columnspan=2)
rosette_frame_input.grid(row=1,column=2)

rosette_button_gc_gauss.grid(row=2,column=0,sticky='n')
rosette_gc_gauss_out.grid(row=2,column=0,sticky='s')
rosette_button_gc_DF.grid(row=2,column=1,sticky='n')
rosette_gc_out_DF.grid(row=2,column=1,sticky='s')
rosette_button_rho_gauss.grid(row=2,column=2,sticky='n')
rosette_rho_gauss_out.grid(row=2,column=2,sticky='s')

rosette_button_Rg_gauss.grid(row=3,column=0,sticky='n')
rosette_Rg_gauss_out.grid(row=3,column=0,sticky='s')
rosette_button_Rh_gauss.grid(row=3,column=1,sticky='n')
rosette_Rh_gauss_out.grid(row=3,column=1,sticky='s')
rosette_button_Rg_DF.grid(row=3,column=2,sticky='n')
rosette_Rg_out_DF.grid(row=3,column=2,sticky='s')

"""A block for displaying input and results for pom-pom polymers"""
#frame
pom_pom_frame=ttk.Frame(window)
pom_pom_frame.place(relx=0.33,rely=0.08,relwidth=0.33,relheight=0.45)
#variables
pom_pom_f1=tk.StringVar(value="0")
pom_pom_f2=tk.StringVar(value="0")
pom_pom_a=tk.StringVar(value="0")
pom_pom_N=tk.StringVar(value="0")

pom_pom_gc_gauss=tk.StringVar(pom_pom_frame)
pom_pom_gc_DF=tk.StringVar()
pom_pom_rho_gauss=tk.StringVar(pom_pom_frame)
pom_pom_Rg_gauss=tk.StringVar(pom_pom_frame)
pom_pom_Rh_gauss=tk.StringVar(pom_pom_frame)
pom_pom_Rg_DF=tk.StringVar(pom_pom_frame)

#Labels and inputs
pom_pom=tk.Label(pom_pom_frame,text="Pom-pom polymer",font="12")
pom_pom_frame_input=ttk.Frame(pom_pom_frame)
pom_pom_f1_text=tk.Label(pom_pom_frame_input,text="Number of arms \n on each pom").pack()
pom_pom_f1_entry=ttk.Entry(pom_pom_frame_input,textvariable=pom_pom_f1,width=5).pack()
pom_pom_f2_text=tk.Label(pom_pom_frame_input,text="Number of arms \n on left pom")
pom_pom_f2_entry=ttk.Entry(pom_pom_frame_input,textvariable=pom_pom_f2,width=5).pack()
pom_pom_a_text=tk.Label(pom_pom_frame_input,text="Relative polimerization \n of the backbone").pack()
pom_pom_a_entry=ttk.Entry(pom_pom_frame_input,textvariable=pom_pom_a,width=8).pack()
pom_pom_N_text=tk.Label(pom_pom_frame_input,text="Polimerization \n degree of a side arm").pack()
pom_pom_N_entry=ttk.Entry(pom_pom_frame_input,textvariable=pom_pom_N,width=8).pack()

#Buttons and outputs
pom_pom_button_gc_gauss=ttk.Button(pom_pom_frame,command=lambda:size_ratio_g_c("pom-pom",pom_pom_gc_gauss,cite,pom_pom_f1,pom_pom_f2,pom_pom_a),text=name_gc)
pom_pom_gc_gauss_out=tk.Label(pom_pom_frame,textvariable=pom_pom_gc_gauss,font="Calibri 12 bold")

pom_pom_button_gc_DF=ttk.Button(pom_pom_frame,command=lambda:size_ratio_rho("pom-pom",pom_pom_rho_gauss,cite,pom_pom_f1,pom_pom_f2,pom_pom_a),text=name_rho)
pom_pom_gc_out_DF=tk.Label(pom_pom_frame,textvariable=pom_pom_rho_gauss,font="Calibri 12 bold")

pom_pom_button_rho_gauss=ttk.Button(pom_pom_frame,command=lambda:size_ratio_g_c_DF("pom-pom",pom_pom_gc_DF,cite,pom_pom_f1,pom_pom_f2,pom_pom_a),text=name_gc_DF)
pom_pom_rho_gauss_out=tk.Label(pom_pom_frame,textvariable=pom_pom_gc_DF,font="Calibri 12 bold")

pom_pom_button_Rg_gauss=ttk.Button(pom_pom_frame,command=lambda:gyration_raduis("pom-pom",pom_pom_Rg_gauss,cite,pom_pom_f1,pom_pom_f2,pom_pom_a,pom_pom_N),text=name_rg)
pom_pom_Rg_gauss_out=tk.Label(pom_pom_frame,textvariable=pom_pom_Rg_gauss,font="Calibri 12 bold")

pom_pom_button_Rh_gauss=ttk.Button(pom_pom_frame,command=lambda:hydrodynamic_raduis("pom-pom",pom_pom_Rh_gauss,cite,pom_pom_f1,pom_pom_f2,pom_pom_a,pom_pom_N),text=name_rh)
pom_pom_Rh_gauss_out=tk.Label(pom_pom_frame,textvariable=pom_pom_Rh_gauss,font="Calibri 12 bold")

pom_pom_button_Rg_DF=ttk.Button(pom_pom_frame,command=lambda:gyration_raduis_DF("pom-pom",pom_pom_Rg_DF,cite,pom_pom_f1,pom_pom_f2,pom_pom_a,pom_pom_N),text=name_rg_DF)
pom_pom_Rg_out_DF=tk.Label(pom_pom_frame,textvariable=pom_pom_Rg_DF,font="Calibri 12 bold")

#picture

pom_pom_c=tk.Canvas(pom_pom_frame,height=200,width=200)
pom_pom_pic=tk.PhotoImage(file="pics/pp1.png",master = pom_pom_c)
pom_pom_c.create_image(0,0,image=pom_pom_pic,anchor="nw")

#grid
pom_pom_frame.columnconfigure((0,1,2),weight=1)
pom_pom_frame.rowconfigure((0,),weight=1)
pom_pom_frame.rowconfigure((1,),weight=1)
pom_pom_frame.rowconfigure((2,3),weight=3)


pom_pom.grid(row=0,column=0,columnspan=3)

pom_pom_c.grid(row=1,column=0,columnspan=2)
pom_pom_frame_input.grid(row=1,column=2)

pom_pom_button_gc_gauss.grid(row=2,column=0,sticky='n')
pom_pom_gc_gauss_out.grid(row=2,column=0,sticky='s')
pom_pom_button_gc_DF.grid(row=2,column=1,sticky='n')
pom_pom_gc_out_DF.grid(row=2,column=1,sticky='s')
pom_pom_button_rho_gauss.grid(row=2,column=2,sticky='n')
pom_pom_rho_gauss_out.grid(row=2,column=2,sticky='s')

pom_pom_button_Rg_gauss.grid(row=3,column=0,sticky='n')
pom_pom_Rg_gauss_out.grid(row=3,column=0,sticky='s')
pom_pom_button_Rh_gauss.grid(row=3,column=1,sticky='n')
pom_pom_Rh_gauss_out.grid(row=3,column=1,sticky='s')
pom_pom_button_Rg_DF.grid(row=3,column=2,sticky='n')
pom_pom_Rg_out_DF.grid(row=3,column=2,sticky='s')

"""A block for displaying input and results for Dumbbell polymers"""
#frame
db_frame=ttk.Frame(window)
db_frame.place(relx=0.66,rely=0.08,relwidth=0.33,relheight=0.45)
#variables
db_a=tk.StringVar(value="0")
db_N=tk.StringVar(value="0")

db_gc_gauss=tk.StringVar(db_frame)
db_gc_DF=tk.StringVar(db_frame)
db_rho_gauss=tk.StringVar(db_frame)
db_Rg_gauss=tk.StringVar(db_frame)
db_Rh_gauss=tk.StringVar(db_frame)
db_Rg_DF=tk.StringVar(db_frame)

#Labels and inputs
db=tk.Label(db_frame,text="Dumbbell polymer",font="12 ")
db_frame_input=ttk.Frame(db_frame)
db_a_text=tk.Label(db_frame_input,text="Relative polimerization \n of the backbone").pack()
db_a_entry=ttk.Entry(db_frame_input,textvariable=db_a,width=8).pack()
db_N_text=tk.Label(db_frame_input,text="Polimerization degree\n of a ring").pack()
db_N_entry=ttk.Entry(db_frame_input,textvariable=db_N,width=8).pack()

#Buttons and outputs
db_button_gc_gauss=ttk.Button(db_frame,command=lambda:size_ratio_g_c("db",db_gc_gauss,cite,db_a),text=name_gc)
db_gc_gauss_out=tk.Label(db_frame,textvariable=db_gc_gauss,font="Calibri 12 bold")

db_button_gc_DF=ttk.Button(db_frame,command=lambda:size_ratio_rho("db",db_rho_gauss,cite,db_a),text=name_rho)
db_gc_out_DF=tk.Label(db_frame,textvariable=db_rho_gauss,font="Calibri 12 bold")

db_button_rho_gauss=ttk.Button(db_frame,command=lambda:size_ratio_g_c_DF("db",db_gc_DF,cite,db_a),text=name_gc_DF)
db_rho_gauss_out=tk.Label(db_frame,textvariable=db_gc_DF,font="Calibri 12 bold")

db_button_Rg_gauss=ttk.Button(db_frame,command=lambda:gyration_raduis("db",db_Rg_gauss,cite,db_a,db_N),text=name_rg)
db_Rg_gauss_out=tk.Label(db_frame,textvariable=db_Rg_gauss,font="Calibri 12 bold")

db_button_Rh_gauss=ttk.Button(db_frame,command=lambda:hydrodynamic_raduis("db",db_Rh_gauss,cite,db_a,db_N),text=name_rh)
db_Rh_gauss_out=tk.Label(db_frame,textvariable=db_Rh_gauss,font="Calibri 12 bold")

db_button_Rg_DF=ttk.Button(db_frame,command=lambda:gyration_raduis_DF("db",db_Rg_DF,cite,db_a,db_N),text=name_rg_DF)
db_Rg_out_DF=tk.Label(db_frame,textvariable=db_Rg_DF,font="Calibri 12 bold")

#picture

db_c=tk.Canvas(db_frame,height=200,width=200)
db_pic=tk.PhotoImage(file="pics/db1.png",master = db_c)
db_c.create_image(0,0,image=db_pic,anchor="nw")

#grid
db_frame.columnconfigure((0,1,2),weight=1)
db_frame.rowconfigure((0,),weight=1)
db_frame.rowconfigure((1,),weight=1)
db_frame.rowconfigure((2,3),weight=3)


db.grid(row=0,column=0,columnspan=3)

db_c.grid(row=1,column=0,columnspan=2)
db_frame_input.grid(row=1,column=2)

db_button_gc_gauss.grid(row=2,column=0,sticky='n')
db_gc_gauss_out.grid(row=2,column=0,sticky='s')
db_button_gc_DF.grid(row=2,column=1,sticky='n')
db_gc_out_DF.grid(row=2,column=1,sticky='s')
db_button_rho_gauss.grid(row=2,column=2,sticky='n')
db_rho_gauss_out.grid(row=2,column=2,sticky='s')

db_button_Rg_gauss.grid(row=3,column=0,sticky='n')
db_Rg_gauss_out.grid(row=3,column=0,sticky='s')
db_button_Rh_gauss.grid(row=3,column=1,sticky='n')
db_Rh_gauss_out.grid(row=3,column=1,sticky='s')
db_button_Rg_DF.grid(row=3,column=2,sticky='n')
db_Rg_out_DF.grid(row=3,column=2,sticky='s')

"""A block for displaying input and results for Snowflake polymers"""
#frame
snow_frame=ttk.Frame(window)
snow_frame.place(relx=0.0,rely=0.55,relwidth=0.33,relheight=0.4)
#variables
snow_fc=tk.StringVar(value="0")
snow_f=tk.StringVar(value="0")
snow_N=tk.StringVar(value="0")

snow_gc_gauss=tk.StringVar(snow_frame)
snow_gc_DF=tk.StringVar()
snow_rho_gauss=tk.StringVar(snow_frame)
snow_Rg_gauss=tk.StringVar(snow_frame)
snow_Rh_gauss=tk.StringVar(snow_frame)
snow_Rg_DF=tk.StringVar(snow_frame)

#Labels and inputs
snow_frame_input=ttk.Frame(snow_frame)
snow=tk.Label(snow_frame,text="Snowflake polymer",font="12")
snow_f_text=tk.Label(snow_frame_input,text="Number of inner arms").pack()
snow_f_entry=ttk.Entry(snow_frame_input,textvariable=snow_f,width=8).pack()
snow_fc_text=tk.Label(snow_frame_input,text="Number of arms per\n outer branching points").pack()

snow_fc_entry=ttk.Entry(snow_frame_input,textvariable=snow_fc,width=8).pack()
snow_N_text=tk.Label(snow_frame_input,text="Polimerization \n degree of an arm").pack()
snow_N_entry=ttk.Entry(snow_frame_input,textvariable=snow_N,width=8).pack()

#Buttons and outputs
snow_button_gc_gauss=ttk.Button(snow_frame,command=lambda:size_ratio_g_c("snow",snow_gc_gauss,cite,snow_f,snow_fc),text=name_gc)
snow_gc_gauss_out=tk.Label(snow_frame,textvariable=snow_gc_gauss,font="Calibri 12 bold")

snow_button_gc_DF=ttk.Button(snow_frame,command=lambda:size_ratio_rho("snow",snow_rho_gauss,cite,snow_f,snow_fc),text=name_rho)
snow_gc_out_DF=tk.Label(snow_frame,textvariable=snow_rho_gauss,font="Calibri 12 bold")

snow_button_rho_gauss=ttk.Button(snow_frame,command=lambda:size_ratio_g_c_DF("snow",snow_gc_DF,cite,snow_f,snow_fc),text=name_gc_DF)
snow_rho_gauss_out=tk.Label(snow_frame,textvariable=snow_gc_DF,font="Calibri 12 bold")

snow_button_Rg_gauss=ttk.Button(snow_frame,command=lambda:gyration_raduis("snow",snow_Rg_gauss,cite,snow_f,snow_fc,snow_N),text=name_rg)
snow_Rg_gauss_out=tk.Label(snow_frame,textvariable=snow_Rg_gauss,font="Calibri 12 bold")

snow_button_Rh_gauss=ttk.Button(snow_frame,command=lambda:hydrodynamic_raduis("snow",snow_Rh_gauss,cite,snow_f,snow_fc,snow_N),text=name_rh)
snow_Rh_gauss_out=tk.Label(snow_frame,textvariable=snow_Rh_gauss,font="Calibri 12 bold")

snow_button_Rg_DF=ttk.Button(snow_frame,command=lambda:gyration_raduis_DF("snow",snow_Rg_DF,cite,snow_f,snow_fc,snow_N),text=name_rg_DF)
snow_Rg_out_DF=tk.Label(snow_frame,textvariable=snow_Rg_DF,font="Calibri 12 bold")

#picture

snow_c=tk.Canvas(snow_frame,height=200,width=200)
snow_pic=tk.PhotoImage(file="pics/snow1.png",master = snow_c)
snow_c.create_image(0,0,image=snow_pic,anchor="nw")

#grid
snow_frame.columnconfigure((0,1,2),weight=1)
snow_frame.rowconfigure((0,),weight=1)
snow_frame.rowconfigure((2,3),weight=3)
snow.grid(row=0,column=0,columnspan=3)


snow_c.grid(row=1,column=0,columnspan=2)
snow_frame_input.grid(row=1,column=2)

snow_button_gc_gauss.grid(row=2,column=0,sticky='n')
snow_gc_gauss_out.grid(row=2,column=0,sticky='s')
snow_button_gc_DF.grid(row=2,column=1,sticky='n')
snow_gc_out_DF.grid(row=2,column=1,sticky='s')
snow_button_rho_gauss.grid(row=2,column=2,sticky='n')
snow_rho_gauss_out.grid(row=2,column=2,sticky='s')

snow_button_Rg_gauss.grid(row=3,column=0,sticky='n')
snow_Rg_gauss_out.grid(row=3,column=0,sticky='s')
snow_button_Rh_gauss.grid(row=3,column=1,sticky='n')
snow_Rh_gauss_out.grid(row=3,column=1,sticky='s')
snow_button_Rg_DF.grid(row=3,column=2,sticky='n')
snow_Rg_out_DF.grid(row=3,column=2,sticky='s')

"""A block for displaying input and results for bottlebrush polymers"""
#frame
bb_frame=ttk.Frame(window)
bb_frame.place(relx=0.33,rely=0.55,relwidth=0.33,relheight=0.4)
#variables
bb_f=tk.StringVar(value="0")
bb_n=tk.StringVar(value="0")
bb_a=tk.StringVar(value="0")
bb_N=tk.StringVar(value="0")

bb_gc_gauss=tk.StringVar(bb_frame)
bb_gc_DF=tk.StringVar()
bb_rho_gauss=tk.StringVar(bb_frame)
bb_Rg_gauss=tk.StringVar(bb_frame)
bb_Rh_gauss=tk.StringVar(bb_frame)
bb_Rg_DF=tk.StringVar(bb_frame)

#Labels and inputs
bb_frame_input=ttk.Frame(bb_frame)
bb=tk.Label(bb_frame,text="Weakly grafted Bottlebrush",font="12")
bb_n_text=tk.Label(bb_frame_input,text="Number of branching point \n and arms per each").pack()
bb_n_entry=ttk.Entry(bb_frame_input,textvariable=bb_n,width=8).pack()
bb_f_text=tk.Label(bb_frame_input,text="Number of arms per \n branching points")
bb_f_entry=ttk.Entry(bb_frame_input,textvariable=bb_f,width=8).pack()
bb_a_text=tk.Label(bb_frame_input,text="Relative polimerization\n degree of a spacer").pack()
bb_a_entry=ttk.Entry(bb_frame_input,textvariable=bb_a,width=8).pack()
bb_N_text=tk.Label(bb_frame_input,text="Polimerization \n degree of a side arm").pack()
bb_N_entry=ttk.Entry(bb_frame_input,textvariable=bb_N,width=8).pack()


#Buttons and outputs
bb_button_gc_gauss=ttk.Button(bb_frame,command=lambda:size_ratio_g_c("bb",bb_gc_gauss,cite,bb_f,bb_n,bb_a),text=name_gc)
bb_gc_gauss_out=tk.Label(bb_frame,textvariable=bb_gc_gauss,font="Calibri 12 bold")

bb_button_gc_DF=ttk.Button(bb_frame,command=lambda:size_ratio_rho("bb",bb_rho_gauss,cite,bb_f,bb_n,bb_a),text=name_rho)
bb_gc_out_DF=tk.Label(bb_frame,textvariable=bb_rho_gauss,font="Calibri 12 bold")

bb_button_rho_gauss=ttk.Button(bb_frame,command=lambda:size_ratio_g_c_DF("bb",bb_gc_DF,cite,bb_f,bb_n,bb_a),text=name_gc_DF)
bb_rho_gauss_out=tk.Label(bb_frame,textvariable=bb_gc_DF,font="Calibri 12 bold")

bb_button_Rg_gauss=ttk.Button(bb_frame,command=lambda:gyration_raduis("bb",bb_Rg_gauss,cite,bb_f,bb_n,bb_a,bb_N),text=name_rg)
bb_Rg_gauss_out=tk.Label(bb_frame,textvariable=bb_Rg_gauss,font="Calibri 12 bold")

bb_button_Rh_gauss=ttk.Button(bb_frame,command=lambda:hydrodynamic_raduis("bb",bb_Rh_gauss,cite,bb_f,bb_n,bb_a,bb_N),text=name_rh)
bb_Rh_gauss_out=tk.Label(bb_frame,textvariable=bb_Rh_gauss,font="Calibri 12 bold")

bb_button_Rg_DF=ttk.Button(bb_frame,command=lambda:gyration_raduis_DF("bb",bb_Rg_DF,cite,bb_f,bb_n,bb_a,bb_N),text=name_rg_DF)
bb_Rg_out_DF=tk.Label(bb_frame,textvariable=bb_Rg_DF,font="Calibri 12 bold")

#picture
bb_c=tk.Canvas(bb_frame,height=200,width=200)
bb_pic=tk.PhotoImage(file="pics/bb1.png",master = bb_c)
bb_c.create_image(0,0,image=bb_pic,anchor="nw")
#grid

bb_frame.columnconfigure((0,1,2),weight=1)
bb_frame.rowconfigure((0,),weight=1)
bb_frame.rowconfigure((1,),weight=1)
bb_frame.rowconfigure((2,3),weight=3)


bb.grid(row=0,column=0,columnspan=3)

bb_c.grid(row=1,column=0,columnspan=2)
bb_frame_input.grid(row=1,column=2)


bb_button_gc_gauss.grid(row=2,column=0,sticky='n')
bb_gc_gauss_out.grid(row=2,column=0,sticky='s')
bb_button_gc_DF.grid(row=2,column=1,sticky='n')
bb_gc_out_DF.grid(row=2,column=1,sticky='s')
bb_button_rho_gauss.grid(row=2,column=2,sticky='n')
bb_rho_gauss_out.grid(row=2,column=2,sticky='s')

bb_button_Rg_gauss.grid(row=3,column=0,sticky='n')
bb_Rg_gauss_out.grid(row=3,column=0,sticky='s')
bb_button_Rh_gauss.grid(row=3,column=1,sticky='n')
bb_Rh_gauss_out.grid(row=3,column=1,sticky='s')
bb_button_Rg_DF.grid(row=3,column=2,sticky='n')
bb_Rg_out_DF.grid(row=3,column=2,sticky='s')




window.mainloop()
