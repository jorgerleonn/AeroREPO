# Resumen Grader 1 — Pórtico 3D de MEF

Documento de traspaso para el siguiente agente que toque `Grader1_MEF.m`.
Última actualización: sesión de trabajo sobre el Grader 1.

---

## 1. Qué es

Script de MATLAB (`Grader1_MEF.m`) que resuelve un **pórtico 3D de 2 nodos por elemento**
mediante el **Método de los Elementos Finitos** (barras 3D con 6 GDL por nodo:
`[u, v, w, θx, θy, θz]`).

- 7 nodos, 6 elementos, 42 GDL.
- El objetivo es reproducir las variables de salida que pide el **grader** de la asignatura.
- Comentarios y variables en español. Las líneas del tipo `consultar página NN de las diapositivas`
  apuntan a apuntes externos (no están en el repo).

### Ficheros relacionados
- `Practicas de MEF/Graders/Grader1_MEF.m` — programa principal.
- `Practicas de MEF/Graders/Funciones_Graders/` — `GaussSeidel.m`, `Geom_viga2D.m`, `Ke_viga2D.m`
  (versiones 2D de referencia; el script define sus propias funciones locales 3D al final).
- `Practicas de MEF/Clases/Clase 4/Clase4profesor.m` — referencia 2D del método (cargas distribuidas,
  condensación estática, reacciones, autovalores).

---

## 2. Salidas que pide el grader

| Variable | Descripción | Estado |
|---|---|---|
| `K_quest` | Matriz de rigidez global (ejes globales, 42×42) | **OK** (grader la da por correcta) |
| `Fneq_quest` | Vector de fuerzas nodales equivalentes (filas, `[1,42]`) | **OK** (distribuidas) |
| `U_quest` | Vector de desplazamientos nodales | **Pendiente** |
| `FR1_quest` | Reacciones en el apoyo 1 | **Pendiente** |
| `FR7_quest` | Reacciones en el apoyo 7 | **Pendiente** |
| `F_uload` | Fuerza equivalente del desplazamiento impuesto `uload` | **Pendiente** |
| `Fe1_quest` | Fuerzas internas **locales** del elemento 1 (vector fila) | **Pendiente** |
| `Fe4_quest` | Fuerzas internas **locales** del elemento 4 (vector fila) | **Pendiente** |
| `Fe6_quest` | Fuerzas internas **locales** del elemento 6 (vector fila) | **Pendiente** |

**Formato importante:** el grader quiere `Fneq_quest` como **fila `[1,42]`** (y pide
`Fe*_quest` "en vector fila"). Internamente conviene calcular en **columna `[42,1]`**
porque `K\u`, `K*U`, `Te'*f_loc` exigen columna; se transpone **solo al final** para la salida
(`Fneq_quest = Fneq';`). En MATLAB fila y columna **no** son intercambiables en `*` y `\`.

---

## 3. Material, sección y cargas

### Sección doble T (función `inercia_viga_doble_T`)
- `b = L/40 = 0.025 m` (ancho de alas), `a = 2b = 0.05 m` (altura), `t = L/200 = 0.005 m` (espesor).
- `[I, A] = inercia_viga_doble_T(a,b,t)` devuelve `I = [Ix, Iy]` **y** el área `A`.
  - `Ix = (b·a³ − (b−t)(a−2t)³)/12`  (eje **fuerte**)
  - `Iy = (2·t·b³ + (a−2t)·t³)/12`  (eje **débil**)
  - `A  = 2·b·t + (a−2t)·t = 4.5e-4 m²`
- Valores: `Ix = 1.5375e-7`, `Iy = 1.3438e-8`.

### Material
- `E = 100e9 Pa`, `nu = 0.3`, `G = E/(2(1+nu)) = 3.8462e10 Pa`.
- `J = (2·b·t³ + (a−2t)·t³)/3 = 3.75e-9` (torsión, sección abierta de pared delgada).
  No hay torsión por las cargas, pero `J` es necesario para armar `Ke`.
- Mapeo de inercias a ejes locales: **`Iz = Ix`** (fuerte→local z), **`Iy = I(2)`** (débil→local y).

### Cargas
- **Distribuida triangular** sobre la columna (elementos 1, 2 y 3, de n1 a n4):
  valor `0` en n1 y máximo `q0y = 500 N/m` en n4, en dirección **`+Y` global**.
- **Puntual `P = 1000 N` en el nodo 5** (dirección/GDL exactos **pendientes de confirmar**).
- **Desplazamiento impuesto `uload = L/100 = 0.01`** (nodo/GDL y signo **pendientes de confirmar**).

### Apoyos (según indicaciones del usuario)
- **Nodo 1 → apoyo simple:** se restringen **desplazamientos** (`u,v,w` = GDL `1,2,3`), giros libres.
- **Nodo 7 → empotramiento:** se restringen los **6 GDL** (`37..42`).
- **OJO:** en el script actual `gdlR = [1,2,3,19,20,21]` es **incorrecto** (apunta al nodo 4).
  Debe ser `gdlR = [1,2,3, 37:42]`.

---

## 4. Geometría y decisiones tomadas

Coordenadas con `L = 1`, `L1=1`, `L2=2`, `L3=4`, `L4=3`:

| Nodo | Coordenada |
|---|---|
| n1 | (1, 0, 0) |
| n2 | (1, 0, 4/3) |
| n3 | (1, 0, 8/3) |
| n4 | (1, 0, 4) |
| n5 | (1, 1, 3.5) |
| n6 | (1, 2, 3) |
| n7 | (0, 2, 3) |

Conectividad: `[1 2; 2 3; 3 4; 4 5; 5 6; 6 7]`.

- **n5 fue modificado a `(L3-L4)/2 + L4 = 3.5`** (original era `(L3-L4)/2 = 0.5`).
  El usuario confirma que **está bien así**: n5 es el punto medio de n4–n6, por lo que los
  elementos **4 y 5 son colineales** (una barra recta partida en dos). Esto es intencionado.
- Longitudes: `Le(1:3) = 4/3 ≈ 1.3333`, `Le(4) = Le(5) = √1.25 ≈ 1.1180`, `Le(6) = 1`.

### `v_ref` por elemento (orienta la sección; es la dirección del **eje local y / alma**)
- Elementos 1–3 (columna vertical): `[0; -1; 0]`.
- Elementos 4 y 5 (colineales): `[0; sin(α); cos(α)]` con `α = atan((L3-L4)/L2) = 26.565°`.
- Elemento 6 (viga horizontal en −X): `[0; 0; 1]`.
- Con estos `v_ref`, **todos salen perpendiculares a su barra (90°)** y `det(Re) = +1`.

---

## 5. Convenios clave (no romper)

1. **GDL del elemento (orden):** `[u1 v1 w1 θx1 θy1 θz1 | u2 v2 w2 θx2 θy2 θz2]`.
   En la estructura, nodo `i` → GDL `6(i−1)+1 .. 6(i−1)+6`.
2. **Transformación:** `Geom_viga3D` devuelve `Te` (12×12) que pasa **global → local**
   (`u_local = Te·u_global`). Por tanto:
   - `K_global_elem = Te' · Ke_local · Te`
   - `f_global = Te' · f_loc`
   - `u_local = Te · U(gdl_e)`
   `Te = blkdiag(Re,Re,Re,Re)` con `Re = [n_x'; n_y'; n_z']`.
3. **`Ke_viga3D` usa la "Forma A"** para flexión en `x-z` (`w`–`θy`), físicamente consistente
   con `θy = −dw/dx`. Se probó la "Forma B" (patrón idéntico al plano `x-y`) y el grader **no**
   la aceptó; **quedarse con la Forma A**.
   - Bloque `x-y` (`v`,`θz`) y bloque `x-z` (`w`,`θy`) con la matriz Hermite estándar.
   - Verificado: simétrica, rango 6, y su parte en el plano coincide exactamente con `Ke_viga2D`.
4. **Las cargas NO afectan a `K`.** Las condiciones de contorno **tampoco** modifican `K`
   (solo se usan en la condensación para resolver).
5. **Cargas distribuidas (vector consistente):** el vector local `f_loc` va en coordenadas
   **locales**; los cortantes en `v` (GDL 2,8) y los momentos en `θz` (GDL 6,12).
   - La carga es `+Y` global, y para la columna `+Y global = −y local`, por eso el vector
     queda con signos negativos (fuerzas) y signos opuestos entre `θz1` y `θz2`.
   - Se usaron las fórmulas estándar **descompuestas** (uniforme `q1` + triangular `dq = q2−q1`):
     ```
     V1 = q1·L/2 + dq·3L/20
     M1 = q1·L²/12 + dq·L²/30
     V2 = q1·L/2 + dq·7L/20
     M2 = q1·L²/12 + dq·L²/20      (el signo global se aplica en f_loc)
     ```
     (Equivalen a `V1=(7q1+3q2)L/20`, `M1=(3q1+2q2)L²/60`, etc.; el `/60` es solo denominador común.)
   - En `Fneq_quest` global, los momentos aparecen en **`θx` global** porque el eje local `z` de
     la columna es `+X` global.
6. **Comillas de salida:** ver §2 (fila `[1,N]`).

---

## 6. Estado por secciones del script

- **§1 Parámetros** — OK (`L,E,b,a,t,I,A,Ix,Iy,nu,G,J,Iz,P,q0y,uload`).
- **§2 Geometría/malla** — OK (nodos, conectividad, `gdlxnod=6`, `n_grados_libertad=42`).
- **§2.1 `Le`, `Te`, `conectividad_g`** — OK. Ojo: dentro de `[...]` en MATLAB un salto de línea
  es **separador de fila**, por eso `gdl_i` usa `...` para concatenar en horizontal.
- **§3 `K_quest`** — **OK**. Ensamblada como `Te_i'·Ke_i·Te_i`. Verificado:
  simétrica (`1.6e-13`) y **rango 36** = 42 − 6 modos de sólido rígido. El grader la acepta.
- **§4 `Fneq` / `Fneq_quest`** — **OK** (solo distribuidas). Verificado: la suma de fuerzas
  globales en Y = **1000 N** (área del triángulo = ½·L3·q0y). `Fneq_quest = Fneq';` (fila).
  Entradas no nulas: GDL `2,4,8,10,14,16,20,22`.
- **§5 (mal rotulada "condiciones de contorno")** — **ROTA**. Usa `K` y `F_total` que **no existen**
  (K se renombró a `K_quest`; `F_total` nunca se definió). Falla en la línea 130.
  Además `gdlR` es incorrecto (§3) y `Ux/Uy/Uz = U(1:3:end)...` es indexado **2D**, incorrecto
  para 3 DOF/nodo; para 3D sería `U(1:6:end), U(2:6:end), U(3:6:end)`.
- **§6 resolución** — vacía (debería contener lo que ahora está mal puesto en §5).
- **§7 autovalores** — vacía.
- **§8 gráficas** — §8.1 solo pinta la estructura sin deformar; falta la deformada (3D).
  `c_nod_def` se declaró como `nnod×2` pero la estructura es 3D.
- **§9 funciones** — `Ke_viga3D`, `Geom_viga3D`, `GaussSeidel`, `inercia_viga_doble_T` (OK).

---

## 7. Pendientes / siguiente trabajo sugerido

1. **Arreglar la resolución (antigua §5 / §6):**
   - Usar `K_quest` (no `K`).
   - Construir el vector de acciones totales:
     `F = Fneq + Fn(P en nodo 5) + F_uload` (o el efecto del desplazamiento en la condensación).
   - `gdlR = [1 2 3, 37:42]`; `uR` con el valor de `uload` en su GDL (0 en el resto).
   - Condensación: `KLL = K_quest(gdlL,gdlL)`, etc.; `uL = KLL\(fL − KLR·uR)`.
2. **`U_quest`** (columna internamente; fila si el grader lo pide).
3. **Reacciones** `FR1_quest` (GDL 1..3) y `FR7_quest` (GDL 37..42):
   `FR = (K_quest·U)(gdlR) − acciones externas en gdlR`.
4. **`F_uload`:** fuerza equivalente del desplazamiento impuesto. Típicamente el vector
   `K_quest(:,gdl_uload)·uload` (o `K_quest·U_imposed`, con `U_imposed` nulo salvo en el GDL de `uload`).
5. **`Fe1_quest`, `Fe4_quest`, `Fe6_quest`:** fuerzas internas **locales** por elemento, en fila:
   `f_loc_e = Ke_viga3D(...) · Te(:,:,e) · U(conectividad_g(e,:))` (fila `[1,12]`).
6. Confirmar con el usuario:
   - GDL y dirección de `P` (nodo 5).
   - GDL y signo de `uload`.
   - Si `Fneq_quest` debe incluir también `P`/`uload` o solo la distribuida (de momento **solo distribuida**).
7. Autovalores/energía (§7) y gráfica de la deformada 3D (§8).

---

## 8. Cómo verificar rápido

```matlab
run('Grader1_MEF.m')
rank(K_quest)                       % 36
norm(K_quest - K_quest')            % ~1e-13
sum(Fneq_quest(2:6:end))            % 1000  (fuerza total en Y)
```

(Con `Fneq_quest` en fila, indexar los GDL de Y con `Fneq_quest(2:6:end)`.)

---

## 9. Resumen de decisiones de la sesión (por si reaparecen dudas)

- `Ke_viga3D`: **Forma A** confirmada (la B no pasó el grader).
- `Ke_viga3D` no era el problema de la `K`; lo que faltaba era el ensamblaje/renombrado y los apoyos.
- Un salto de línea dentro de `[...]` en MATLAB = **nueva fila** (usar `...` para concatenar en horizontal).
- Cargas en **ejes locales** con `Te'` a global; los apoyos no tocan `K`.
- Los `*_quest` que pide el grader van en **fila**; el cálculo interno en **columna**.
