
from traceback import print_tb
import numpy as np

carga_inicial = 20
diferencial_carga = 0.25
carga_final = carga_inicial + 18 * diferencial_carga

distancia_AD = 1
distancia_AB = 2
distancia_BC = 2
distancia_AC = distancia_AB+ distancia_BC

ángulo_barra_inclinada = np.arctan(distancia_AD / distancia_BC)
ángulo_barra_inclinada_AE = np.arctan(distancia_AD / distancia_AB)

rango_de_ángulos = np.arange(0, 190, 10)
rango_de_cargas = np.arange(carga_inicial, carga_final + diferencial_carga, diferencial_carga)
rango_de_cargas
print(rango_de_cargas)

def cálculo_reacciones(carga, ángulo, distacia_AC):
    ángulo_radianes = (ángulo if ángulo <= 90 else 180 - ángulo) * np.pi / 180
    Ay = (-1 if ángulo <= 90 else 1) * carga * np.cos(ángulo_radianes)
    Dx = (-1 if ángulo <= 90 else 1) * distancia_AC * carga * np.cos(ángulo_radianes)
    Ax = - Dx - carga * np.sin(ángulo_radianes)
    return Ax, Ay, Dx

def  cálculo_fuerzas_nodo_C(carga, ángulo):
    ángulo_radianes =(ángulo if ángulo <= 90 else 180 - ángulo) * np.pi / 180
    CE = (1 if ángulo <= 90 else -1) * carga * np.cos(ángulo_radianes) / np.sin(ángulo_barra_inclinada)
    BC = carga * np.sin(ángulo_radianes) - CE * np.cos(ángulo_barra_inclinada)
    return CE, BC

def cálculo_fuerzas_nodo_B(fuerza_BC):
    BA = fuerza_BC
    BE = 0
    return BA, BE

def cálculo_fuerzas_nodo_E(fuerza_CE, ángulo_barra_inclinada_AE):
    EA = fuerza_CE * np.sin(ángulo_barra_inclinada) / np.sin(ángulo_barra_inclinada_AE)
    ED = EA * np.cos(ángulo_barra_inclinada_AE) + fuerza_CE * np.cos(ángulo_barra_inclinada)
    return EA, ED

def cálculo_fuerzas_nodo_D(fuerza_Dx):
    ED = - fuerza_Dx
    DA = 0
    return ED, DA

def cálculo_fuerzas_nodo_A(Ax, Ay, ángulo_barra_inclinada_AE):
    AE = Ay / np.sin(ángulo_barra_inclinada_AE)
    AB = Ax / np.cos(ángulo_barra_inclinada_AE)
    return AE, AB

fuerza_maxima = np.array([0, 0, 0, 0, 0, 0, 0])
for ángulo, carga in zip(rango_de_ángulos, rango_de_cargas):
    Ax, Ay, Dx = cálculo_reacciones(carga, ángulo, distancia_AC)
    CE, BC = cálculo_fuerzas_nodo_C(carga, ángulo)
    BA, BE = cálculo_fuerzas_nodo_B(BC)
    EA, ED = cálculo_fuerzas_nodo_E(CE, ángulo_barra_inclinada_AE)
    DE, DA = cálculo_fuerzas_nodo_D(Dx)
    AE, AB = cálculo_fuerzas_nodo_A(Ax, Ay, ángulo_barra_inclinada_AE)
    print(f"Reacciones ángulo", ángulo, "grados: ")
    print(Ax, Ay, Dx)
    print("Fuerzas internas: BA, EA, DA, BC, BE, CE, ED")
    print(BA, EA, DA, BC, BE, CE, ED)
    print("AE == EA,  AB == BA,  DE == ED")
    print(f"{AE} == {EA}, {AB} == {BA}, {DE} == {ED}")
    print()
    fuerza_maxima[0] = BA if np.absolute(BA) > fuerza_maxima[0] else fuerza_maxima[0]
    fuerza_maxima[1] = BC if np.absolute(BC) > fuerza_maxima[1] else fuerza_maxima[1]
    fuerza_maxima[2] = CE if np.absolute(CE) > fuerza_maxima[2] else fuerza_maxima[2]
    fuerza_maxima[3] = ED if np.absolute(ED) > fuerza_maxima[3] else fuerza_maxima[3]
    fuerza_maxima[4] = EA if np.absolute(EA) > fuerza_maxima[4] else fuerza_maxima[4]
    fuerza_maxima[5] = DA if np.absolute(DA) > fuerza_maxima[5] else fuerza_maxima[5]
    fuerza_maxima[6] = BE if np.absolute(BE) > fuerza_maxima[6] else fuerza_maxima[6]

print(fuerza_maxima)

######
## Esfuerzo admisible en MPa
esfuerzo_admisible = 21
## Área de la sección en mm2
área = 1200
## Fuerza máxima en KN
fuerza_maxima = esfuerzo_admisible * área / 1000
fuerza_maxima
print("Fuerza máxima = ", fuerza_maxima)

fuerzas_internas = np.array([0, 0, 0, 0, 0, 0, 0])
for ángulo, carga in zip(rango_de_ángulos, rango_de_cargas):
    Ax, Ay, Dx = cálculo_reacciones(carga, ángulo, distancia_AC)
    CE, BC = cálculo_fuerzas_nodo_C(carga, ángulo)
    BA, BE = cálculo_fuerzas_nodo_B(BC)
    EA, ED = cálculo_fuerzas_nodo_E(CE, ángulo_barra_inclinada_AE)
    DE, DA = cálculo_fuerzas_nodo_D(Dx)
    AE, AB = cálculo_fuerzas_nodo_A(Ax, Ay, ángulo_barra_inclinada_AE)
    fuerzas_internas = np.array([("AB", BA), ("AE", EA), ("AD", DA), ("BC", BC), ("BE", BE), ("CE", CE), ("ED", ED)])
    elemento_falla = np.where(np.absolute(np.array([x[1] for x in fuerzas_internas], dtype=np.float64)) > fuerza_maxima)
    if len(elemento_falla) > 0:
        print(f"la primera falla se presenta en la carga de: ", carga, "KN en la dirección de ", ángulo, "grados")
        for i in elemento_falla[0]:
            fuerza = fuerzas_internas[i][1].astype(float)
            causa_falla = "compresión" if fuerza < 0 else "tracción"
            print(f"El elemento {fuerzas_internas[i][0]} falla bajo la carga {np.absolute(fuerza)} KN a {causa_falla}")
        print()

#### TERCER PUNTO
## Código 2174157, valor de L tomado es el 7,porque el quinto valor es 1
L = 7
aumento_distancia = float("0.0" + str(L))
A = 5
B = 7
AB = 60
AB = round(AB / 10) * 10
AB, aumento_distancia
print("Ángulo AB:",AB,",","Aumento de distancia de:" ,aumento_distancia)

for ángulo, carga in zip(rango_de_ángulos, rango_de_cargas):
    distancia_AC += aumento_distancia
    distancia_AB += aumento_distancia
    ángulo_barra_inclinada_AE = np.arctan(distancia_AD / distancia_AB)
    Ax, Ay, Dx = cálculo_reacciones(carga, ángulo, distancia_AC)
    CE, BC = cálculo_fuerzas_nodo_C(carga, ángulo)
    BA, BE = cálculo_fuerzas_nodo_B(BC)
    EA, ED = cálculo_fuerzas_nodo_E(CE, ángulo_barra_inclinada_AE)
    DE, DA = cálculo_fuerzas_nodo_D(Dx)
    AE, AB = cálculo_fuerzas_nodo_A(Ax, Ay, ángulo_barra_inclinada_AE)
    fuerzas_internas = [("AB", BA, carga), ("AE", EA, carga), ("AD", DA, carga), ("BC", BC, carga), ("BE", BE, carga), ("CE", CE, carga), ("ED", ED, carga)]
    print(f"Aumento de distancia de:", distancia_AB - 2, "Las fuerzas calculadas son:")
    print(fuerzas_internas)
    print()
import python as np
import sympy as sp

# 1.0 UNIDADES 
# 1.1 METROS
longitud_cables = 1.5
# 1.2 PULGADAS
diametro_cables = 1 / 4
# 1.3 MILIMETROS CUADRADOS
area_cables = np.pi * (diametro_cables * 25.4 / 1000) ** 2 / 4
# 1.4 PULGADAS
diametro_pasadores = 1 / 2
# 1.5 MILIMETROS CUADRADOS
area_pasadores = np.pi * (diametro_pasadores * 25.4 / 1000) ** 2 / 4

Codigo = "2174157"
x = sum([int(i) for i in Codigo]) 
# BARRA EN METROS
L = 0.2 * x

FS_fluencia = float(f"1.{x}")
FS_resistencia_Ult = float(f"1{x + 3}")

angulo_BC = 40 * np.pi / 180
angulo_BE = 50 * np.pi / 180

# PROPIEDADES DEL MATERIAL
# UNIDADES EN MPa
esfuerzo_ultimo = 400
esfuerzo_fluencia_tracion = 250
esfuerzo_fluencia_compresion = 145
# GPa
E = 200


# SUMATORIA DE FUERZAS EN Y AND X
Ax , Ay , BC , BD , BE , W = sp.symbols('Ax Ay BC BD BE W')
# SUMATORIA EN X
Fx = Ax - BC * sp.sin(angulo_BC) + BE * sp.sin(angulo_BE)
equilibrio_X = sp.Eq(Fx, 0)
equilibrio_X
print(equilibrio_X)

# SUMATORIA EN Y
Fy = Ay - W * L + BC * sp.cos(angulo_BC) + BE * sp.cos(angulo_BE) + BD
equilibrio_Y = sp.Eq(Fy, 0)
equilibrio_Y


# SUMATORIA DE MOMENTO
M_a = W * L **2 / 2 + BC * sp.cos(angulo_BC) * L + BE * sp.cos(angulo_BE) * L + BD * L
equilibrio_momento = sp.Eq(M_a, 0)
equilibrio_momento


# ECUACION DE COMPATIBILIDAD POR DESPLAZAMIENTO
igualando_BC = sp.Eq(BC - BD * sp.cos(angulo_BC), 0)
ecuacion_igualacion_BC = BD * sp.cos(angulo_BC)
igualando_BC

igualando_BE = sp.Eq(BE - BD * sp.cos(angulo_BE), 0)
ecuacion_igualacion_BE = BD * sp.cos(angulo_BE)
igualando_BE

# DESPEJAR BD DE LA SUMATORIA DE MOMENTO
ecuacion = sp.Eq(M_a.subs([(BC, ecuacion_igualacion_BC), (BE, ecuacion_igualacion_BE)]), 0)
ecuacion


BD_vs_W = sp.solve(ecuacion)[0][BD]
BD_vs_W


F_BC = ecuacion_igualacion_BC.subs(BD, BD_vs_W)
F_BE = ecuacion_igualacion_BE.subs(BD, BD_vs_W)

fuerza_cables = [F_BC, BD_vs_W, F_BE]
fuerza_cables

w_max = []
for fuerza in fuerza_cables:
    # ESFUERZO AXIAL DE LOS CABLES
    ecuacion_esfuerzo_axial = sp.Eq(fuerza / area_cables, esfuerzo_ultimo * 10 ** 6 / FS_resistencia_Ult)
    w_cables = sp.solve(ecuacion_esfuerzo_axial)[0] / 1000
    # PARA LOS PASADORES SUPERIORES EN CORTANTE SIMPLE
    ecuacion_esfuerzo_cortante = sp.Eq(fuerza / area_pasadores, esfuerzo_fluencia_compresion * 10 ** 6 / FS_fluencia)
    w_pasadores = sp.solve(ecuacion_esfuerzo_cortante)[0] / 1000
    w_max.append(w_cables)
    w_max.append(w_pasadores)

print("carga distribuida maxima en los cables [kN / m]", w_cables)
min(w_max)

# PARA LOS APOYOS EN CORTANTE DOBLE
Ay_eq = Fy.subs([(BC, ecuacion_igualacion_BC), (BE, ecuacion_igualacion_BE), (BD, BD_vs_W)])
Ay_vs_w = sp.solve(sp.Eq(Ay_eq, 0))
Ay_vs_w = Ay_vs_w[0][Ay]
Ay_vs_w
print(Ay_vs_w)

Ax_eq = Fx.subs([(BC, ecuacion_igualacion_BC), (BE, ecuacion_igualacion_BE), (BD, BD_vs_W)])
Ax_vs_w = sp.solve(sp.Eq(Ax_eq, 0))
Ax_vs_w = Ax_vs_w[0][Ax]
Ax_vs_w
print(Ax_vs_w)


# LAS FUERZAS EN EL APOYO A
F_apoyo_A = sp.sqrt(Ay_vs_w ** 2 + Ax_vs_w ** 2)
F_apoyo_A
print(F_apoyo_A)

# FUERZA DE CONEXION DE LOS CABLES EN EL PUNTO 
F_x_cortante_doble = ecuacion_igualacion_BE.subs(BD, BD_vs_W) * sp.sin(angulo_BE) - ecuacion_igualacion_BC.subs(BD, BD_vs_W) * sp.sin(angulo_BC)
F_y_cortante_doble = ecuacion_igualacion_BE.subs(BD, BD_vs_W) * sp.sin(angulo_BE) + ecuacion_igualacion_BC.subs(BD, BD_vs_W)* sp.sin(angulo_BC)
F_apoyo_B = sp.sqrt(F_x_cortante_doble * 2 + F_y_cortante_doble * 2)
F_apoyo_B
print(F_apoyo_B)


fuerzas_apoyo = [F_apoyo_A, F_apoyo_B]
w_max_apoyo = []
for fuerza in fuerzas_apoyo:
    ecuacion_esfuerzo_cortante = sp.Eq(fuerza / area_pasadores, esfuerzo_fluencia_compresion * 10 ** 6 / FS_fluencia)
    w_pasadores = sp.solve(ecuacion_esfuerzo_cortante)[0] / 1000
    w_max_apoyo.append(w_pasadores)

print("carga distribuida maxima por pasadores (kN / m",)
carga_maxima = min(w_max)
carga_maxima
print(carga_maxima)


# CALCULO DE FUERZAS Y DEFORMACIONES

valor_BE = ecuacion_igualacion_BE.subs([(BD, BD_vs_W), (W, carga_maxima)])
deformacion_BE = (valor_BE * 1000 * longitud_cables) / (area_cables * E * 10 ** 9)
# FUERZA EN kN y DEFORMACION mm
valor_BE, deformacion_BE * 1000
print(deformacion_BE)
print(valor_BE)

valor_BC = ecuacion_igualacion_BC.subs([(BD, BD_vs_W),(W, carga_maxima)])
deformacion_BC = (valor_BC * 1000 * longitud_cables) / (area_cables * E * 10 ** 9)
valor_BC, deformacion_BC * 1000
print(valor_BC)
print(deformacion_BC)

valor_BD = BD_vs_W.subs(W, carga_maxima)
deformacion_BD = (valor_BD * 1000 * longitud_cables) / (area_cables * E * 10 ** 9)
valor_BD, deformacion_BD * 1000
print(valor_BD)
print(deformacion_BD)

valor_Ax = Fx.subs([(BE, valor_BE), (BC, valor_BC)])
valor_Ax = sp.solve(sp.Eq(valor_Ax, 0))[0]
valor_Ax
print(valor_Ax)

valor_Ay = Fy.subs([(BE, valor_BE), (BC, valor_BC), (BD, valor_BD), (W, carga_maxima)])
valor_Ay = sp.solve(sp.Eq(valor_Ay, 0))[0]
valor_Ay
print(valor_Ay)










