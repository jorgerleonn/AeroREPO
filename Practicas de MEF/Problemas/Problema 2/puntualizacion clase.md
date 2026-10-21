El Jacobiano no es constante en el interior del elemento, 

El Jacobiano en cada pnto del elemento tiene un valor, 
pues hay un parámetro en MEF que es la relación entre el Jacobiano minimo y el máximo del elemento
esto lo que me dice es que el Jacobiano no tiene el mismo valor dentro del elemento.

El Jacobiano (o más rigurosamente, la matriz jacobiana $J$) **no representa una homotecia general**, sino la **mejor aproximación lineal local** de una transformación diferenciable $\mathbf{f}: \mathbb{R}^n \to \mathbb{R}^m$.

En términos del elemento diferencial:

* **Elemento de línea (vectorial):** El vector diferencial $d\mathbf{x}$ se transforma mediante la acción directa de la matriz jacobiana:

$$d\mathbf{y} = J \, d\mathbf{x}$$


donde $J_{ij} = \frac{\partial f_i}{\partial x_j}$.
* **Dilatación / Escalamiento local:**
* **Elemento de longitud escalar ($ds$):** No escala uniformemente salvo que la transformación sea conforme. En general, depende de la dirección:

$$ds^2 = d\mathbf{y}^T d\mathbf{y} = d\mathbf{x}^T (J^T J) d\mathbf{x}$$


donde el tensor métrico inducido es $g = J^T J$.
* **Elemento de volumen ($dV$):** Aquí sí actúa como factor de escala escalar global local mediante el **determinante jacobiano**:

$$dV_y = \vert{}\det(J)\vert{} \, dV_x$$


Para que la transformación del elemento de línea sea estrictamente una **homotecia local** (mismo factor de escala en todas las direcciones y conservación de ángulos), $J$ debe satisfacer la condición de conformalidad:


$$J^T J = \lambda^2 I \quad (\lambda > 0)$$


donde $\lambda = \vert{}\det(J)\vert{}^{1/n}$ es la razón de la homotecia local.