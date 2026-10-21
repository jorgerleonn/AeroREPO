# Resumen Problema 2 — Pórtico 3D con apoyo móvil sobre plano inclinado

Documento de traspaso para el siguiente agente que toque `Problema2_MEF.m`.
Última actualización: primera resolución del Problema 2.

---

## 1. Qué es

Variante del **Problema 1** (`Graders/Grader 1.1/Grader1_1.m`). Es **la misma
estructura** (mismo material, sección, malla y cargas), con **una sola
modificación**:

> El empotramiento del **nodo 1** se transforma en un **empotramiento móvil
> sobre un plano inclinado perpendicular al eje de la viga del elemento 1**.

Todo lo demás es invariante respecto del Problema 1.

- 17 nodos, 16 elementos, 102 GDL (6 por nodo: `[u, v, w, θx, θy, θz]`).
- Método: barras 3D de 2 nodos, rigidez local `Ke_viga3D` y transformación
  `Te` de `Geom_viga3D` (idénticas a `Grader1_1.m`).

### Ficheros relacionados
- `Graders/Grader 2/Problema2_MEF.m` — programa principal (esta solución).
- `Graders/Grader 1.1/Grader1_1.m` — problema base del que se parte.
- `Graders/Grader 1/ResumenGrader1.md` — notas del Grader 1 (2D/otra estructura),
  útil para entender convenios generales del método.

---

## 2. Salidas que pide el enunciado

| Variable | Descripción | Formato |
|---|---|---|
| `U_quest` | Vector de desplazamientos nodales (mm) | **columna `[102×1]`** |
| `FR1_quest` | Vector de reacciones en el apoyo 1 (nodo 1) | **columna `[6×1]`** |

Se usan exactamente esos nombres. En este problema **no** hay desplazamiento
impuesto (no existe `uload` como en el Grader 1). El enunciado pide no usar
`clear all` (el script usa sólo `close all;`).

---

## 3. Material, sección y cargas (invariantes del Problema 1)

### Sección cajón rectangular (función `inercia_viga_cajon_rectangular`)
- `b = L/40 = 25 mm` (ancho), `a = 2b = L/20 = 50 mm` (altura), `t = L/400 = 2.5 mm`.
- Devuelve `I = [Iy, Iz]` **y** el área `A`:
  - `Iy = (b·a³ − (b−2t)(a−2t)³)/12` — **eje fuerte** (a paralelo al eje local Y).
  - `Iz = (a·b³ − (a−2t)(b−2t)³)/12` — eje débil.
  - `A = a·b − (a−2t)(b−2t) = 350 mm²`.
- En el script: `Iz = I(1)` (fuerte) e `Iy = I(2)` (débil). **Ojo con los nombres:
  están intercambiados respecto a la notación geométrica.** No "arreglarlos":
  el bloque `x-y` (`v,θz`) usa `Iz = I(1)` (fuerte) y el bloque `x-z` (`w,θy`)
  usa `Iy = I(2)` (débil); así se replicó el problema 1.

### Material
- `E = 100e3 N/mm²`, `nu = 1/3`, `G = E/(2(1+nu))`.
- `J = 2·t·(a−t)²·(b−t)²/(a+b−2t)` (Bredt, sección cerrada de pared delgada).

### Cargas
- `q0x = 5 N/mm` distribuida en **+X global**, perfil triangular `q(z) = q0x·z/L`.
  - columna n1–n4 (elementos 1–3): 0 en n1 → `q0x` en n4.
  - elementos a `z = L` (7,8,9,10): constante `q0x`.
  - ramas verticales (11–13 y 14–16): `q0x` arriba → 0 abajo.
  - **la viga superior n4–n7 (elementos 4,5,6) NO lleva carga distribuida**
    (se excluye con `ismember(e,[4 5 6])`).
  - Total de carga distribuida: `12795.085 N` en +X.
- `P1 = 10000 N` en **−Z** en el nodo 8 (`P1` al GDL `6·7+3 = 45`).
- `P2 = 5000 N` en **−Z** en el nodo 13 (`P2` al GDL `6·12+3 = 75`).
- `F_total = Fn + Fneq` (puntuales + distribuidas). `Fneq_quest = Fneq` (sólo distribuidas).

---

## 4. El apoyo del nodo 1 (núcleo del problema 2)

El eje de la viga del elemento 1 (n1→n2, `(0,0,0)→(L/6,0,L/3)`) es:

```
e1 = (1, 0, 2)/√5     → normal al plano inclinado
```

Ése es el eje perpendicular al plano. Base local del apoyo:
- `e1 = (1,0,2)/√5`  → **traslación restringida** (normal al plano).
- `e2 = (0,1,0)`     → eje Y global.
- `e3 = e1 × e2 = (-2,0,1)/√5` → 2ª dirección libre.

Interpretación **literal** del enunciado (3 restricciones en el nodo 1):

| GDL nodo 1 | Estado | Justificación |
|---|---|---|
| traslación según `e1` (`e1·(u,v,w)=0`) | **restringido** | "restringe todos los demás movimientos de traslación" |
| traslación `v` (Y) y según `e3` | **libres** | "permite el movimiento traslacional paralelo al plano" |
| giro `θy` (eje ⊥ al plano de simetría X–Z) | **libre** | "articulado… alrededor de su eje perpendicular" |
| giros `θx`, `θz` | **restringidos** | "impide el giro en el resto de planos-ejes" |

**Nota:** bloquear los giros `θx` y `θz` globales equivale a bloquear los dos
giros contenidos en el plano X–Z (da igual plantearlo en ejes del apoyo o
globales, porque el giro libre es justo alrededor de Y). Por eso el bloque de
giros de `T_apoyo` se deja en **ejes globales** y sólo las **traslaciones** se
rotan a la base del plano.

### Implementación (transparencia de coordenadas)
- `R_tras = [e1, e2, e3]` (columnas = ejes locales en globales).
- `T_apoyo = blkdiag(R_tras, eye(3))` y `T_bc = eye(102)` con `T_bc(1:6,1:6)=T_apoyo`.
- `K_t = T_bc'·K_quest·T_bc`, `F_t = T_bc'·F_total`.
- `gdl_apoyo = [1, 4, 6, 6·(12-1)+(1:6), 6·(17-1)+(1:6)]`
  (los GDL `1,4,6` son locales del nodo 1; `67:72` y `97:102` son los
  empotramientos de n12 y n17, que **no cambian**).
- Resolver en la base transformada y recuperar `U_quest = T_bc·U_t`.
- Reacciones: `FR1_quest = (K_quest·U_quest)(1:6) − F_total(1:6)` (globales).

---

## 5. Geometría, malla y `v_ref`

Coordenadas (`L = 1000`):

| Nodo | Coordenada |
|---|---|
| n1 | (0, 0, 0) |
| n2 | (L/6, 0, L/3) |
| n3 | (L/3, 0, 2L/3) |
| n4 | (L/2, 0, L) |
| n5 | (5L/6, 0, L) |
| n6 | (7L/6, 0, L) |
| n7 | (3L/2, 0, L) |
| n8 | (3L/2, L/4, L) |
| n9 | (3L/2, L/2, L) |
| n10 | (3L/2, L/2, 2L/3) |
| n11 | (3L/2, L/2, L/3) |
| n12 | (3L/2, L/2, 0) |
| n13 | (3L/2, −L/4, L) |
| n14 | (3L/2, −L/2, L) |
| n15 | (3L/2, −L/2, 2L/3) |
| n16 | (3L/2, −L/2, L/3) |
| n17 | (3L/2, −L/2, 0) |

Conectividad: `[1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8; 8 9; 7 13; 13 14; 9 10; 10 11; 11 12; 14 15; 15 16; 16 17]`.

`v_ref` (dirección del eje local y / "arriba" de la sección):
- 1:3 `[-1;0;0]` (columna inclinada), 4:6 `[0;0;1]` (viga superior),
- 7:8 y 9:10 `[0;0;1]` (ramas horizontales ±Y),
- 11:13 y 14:16 `[0;-1;0]` (ramas verticales).

---

## 6. Convenios clave (no romper)

1. GDL del elemento: `[u1 v1 w1 θx1 θy1 θz1 | u2 v2 w2 θx2 θy2 θz2]`.
   Nodo `i` → GDL `6(i−1)+1 .. 6(i−1)+6`.
2. `Geom_viga3D` devuelve `Te` (global→local, 12×12): `Te = blkdiag(Re,Re,Re,Re)`.
   - `K_global_elem = Te'·Ke·Te`, `f_global = Te'·f_loc`.
3. `Ke_viga3D` usa la "Forma A" (idéntica al Grader 1, que pasó el grader).
4. Las cargas no afectan a `K_quest`; las condiciones de contorno tampoco.
5. En MATLAB un salto de línea dentro de `[...]` = **nueva fila** (usar `...`).

---

## 7. Estado del script

- §1 Parámetros, §2 geometría/malla, §2.1 `Le/Te/conectividad_g`: **OK** (igual que Problema 1).
- §3 `K_quest`: **OK** (ensamblada `Te'·Ke·Te`; `rank(K_quest) = 96` = 102 − 6 modos rígidos).
- §4 `Fneq` / `Fneq_quest` / `F_total`: **OK** (total distribuida +X = 12795.085 N; P1 y P2 en −Z).
- §5 condiciones de contorno: **OK** — apoyo móvil del nodo 1 vía `T_bc` + empotramientos n12 y n17.
- §6 resolución: **OK** — `U_quest` (columna) y `FR1_quest` (columna). `rank(K_LL) = 87/87` (estable).
- §7 gráfica 3D original vs deformada: **OK**.

---

## 8. Resultados obtenidos (primera resolución)

`FR1_quest` (global, N y N·mm):

```
[ -1033.26648 ;      % Fx
        0        ;    % Fy  (GDL libre)
  -2066.53296 ;      % Fz
   5171.91051 ;      % Mx
        0        ;    % My  (GDL libre)
   -131.69357 ]      % Mz
```

- La fuerza resultante = `2310.45 N`, **exactamente según `e1 = (1,0,2)/√5`**.
  En la base del apoyo `FR1 = [-2310.454, 0, 0]` (sólo componente axial).
- `Fy = My = 0` (coherente con los GDL libres: traslación Y y giro θy).
- Verificación de restricciones en n1: `e1·u = 0`, `θx = 0`, `θz = 0` (residuo 0).
- Equilibrio global: reacciones + cargas externas = 0 (X: 12795.08, Z: 15000).
- `max|U_quest| = 472.887 mm` (nodo 1, GDL 1). El nodo 1 desliza mucho porque
  el apoyo es muy poco restrictivo (sólo coarta axial + giros X–Z).

Nodos destacados de `U_quest` (mm; columnas u, v, w, θx, θy, θz):

| Nodo | u | v | w | θx | θy | θz |
|---|---|---|---|---|---|---|
| 1 | 472.887 | −2.304 | −236.444 | 0 | −0.217 | 0 |
| 2 | 400.428 | −2.212 | −200.187 | ~0 | −0.217 | ~0 |
| 4 | 260.245 | −1.470 | −130.055 | ~0 | −0.193 | ~0 |
| 9 | 238.081 | −1.526 | −0.281 | 0.0286 | 0.2326 | 0.0632 |
| 12 | 0 | 0 | 0 | 0 | 0 | 0 |
| 17 | 0 | 0 | 0 | 0 | 0 | 0 |

---

## 9. Pendientes / a confirmar

1. **Interpretación del apoyo (lo único discutible).** Se ha tomado el
   enunciado literal: 2 traslaciones libres (todo el plano) + 1 giro libre
   (`θy`), 3 GDL restringidos. Si la asignatura interpreta el "empotramiento
   móvil" de forma más restrictiva, basta cambiar `gdl_apoyo`:
   - bloquear además la traslación `v` (Y): añadir el GDL local `2`.
   - bloquear todos los giros: añadir el GDL local `5`.
   Con cargas simétricas en Y, bloquear o no la traslación Y apenas cambia el
   resultado (ese desplazamiento es ≈0), pero bloquear `θy` sí cambia bastante.
2. **Convenio de salida de `FR1_quest`.** Se entrega en **ejes globales**
   (igual que el Problema 1). Si el grader lo quiere en la **base del apoyo**,
   es `[e1;e2;e3]'·FR1_quest(1:3)` para la fuerza (da `[-2310.454;0;0]`) más
   los momentos `θ` en la base correspondiente.
3. **Formato fila vs columna.** Aquí `U_quest` y `FR1_quest` se entregan en
   **columna**. El Grader 1 pedía `Fneq_quest`/`Fe*_quest` en fila; conviene
   confirmar el formato que espera el grader del Problema 2.

---

## 10. Cómo verificar rápido

```matlab
% cd a Practicas de MEF/Graders/Grader 2
run('Problema2_MEF.m')
rank(K_quest)                                   % 96
rank(K_t(gdl_libres,gdl_libres))                % 87  (= nº de GDL libres)
e1_apoyo'*U_quest(1:3)                          % ~0  (traslación axial nula)
[U_quest(4); U_quest(6)]                        % ~0  (θx y θz nulos)
% FR1 en la base del apoyo:
[e1_apoyo'; e2_apoyo'; e3_apoyo']*FR1_quest(1:3)  % [-2310.454; 0; 0]
% Equilibrio: reacciones (n1+n12+n17) + cargas externas = 0
```
