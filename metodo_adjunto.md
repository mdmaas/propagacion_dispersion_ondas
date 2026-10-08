# Plan de implementación — Unidad 5: El método adjunto

Documento de trabajo que describe el contenido y la implementación de `unidad-5.ipynb`.
El notebook es una clase completa (MyST/Jupyter) con teoría, derivaciones y scripts de
Python ejecutables y validados.

## 1. Objetivo

Partiendo de lo visto en la unidad anterior (ajuste de un modelo de velocidades
parametrizado por un polinomio de pocos coeficientes, resuelto con diferencias finitas
y Gauss–Newton), mostrar por qué ese enfoque **no escala** cuando el número de parámetros
$p$ es enorme (FWI, optimización de formas, diseño inverso) y presentar el **método
adjunto** como la herramienta que calcula el gradiente completo con **una** resolución
del problema directo y **una** del adjunto, sin importar $p$.

## 2. Estructura de la clase (celdas del notebook)

Hilo pedagógico: **todos los casos con EDP se trabajan primero en dimensión finita**
(discretizar → aplicar el adjunto → interpretar la ecuación adjunta discreta como ecuación
diferencial adjunta) y **recién al final** se discute el punto de vista en dimensión infinita.

### Bloque 0 — Motivación
- Markdown: recordatorio de la unidad previa (pocos parámetros, Gauss–Newton) y el muro
  computacional cuando $p \sim 10^6$–$10^9$.
- Idea: calcular $\nabla_m J$ sin formar la matriz de sensibilidad $\partial x/\partial m$.

### Bloque 1 — Fundamento: multiplicadores de Lagrange
Problema restringido
$$\min_{x,m} J(x,m)\quad \text{s.a. } g(x,m)=0,\qquad x\in\mathbb R^n,\ m\in\mathbb R^p .$$
- Reformulación como problema sin restricciones $\min_m \hat J(m)=J(x(m),m)$ (teorema de la
  función implícita) y la regla de la cadena
  $\dfrac{d\hat J}{dm}=\dfrac{\partial J}{\partial m}+\dfrac{\partial J}{\partial x}\dfrac{dx}{dm}$.
- Lagrangiano $\mathcal L(x,m,\lambda)=J(x,m)+\lambda^{T}g(x,m)$.
- Condiciones de estacionariedad:
  1. $\nabla_\lambda\mathcal L=g(x,m)=0$ (problema **directo**).
  2. $\nabla_x\mathcal L=\nabla_x J+(Dg)^T\lambda=0$ (problema **adjunto**: sistema lineal).
  3. Derivada total $d\hat J/dm$ y **cancelación** del término de sensibilidad
     $\partial x/\partial m$ gracias a la ecuación adjunta.
- Resultado final:
  $$\boxed{\ \frac{d\hat J}{dm}=\frac{\partial J}{\partial m}+\lambda^{T}\frac{\partial g}{\partial m}\ }$$
  Costo: resolver directo + adjunto + evaluación de derivadas parciales. **Independiente de $p$.**

### Bloque 2 — Ejemplo de control discreto (completo)
- Sistema dinámico discreto $x_{k+1}=a x_k+b m_k$, $x_0$ fijo, $k=0,\dots,N-1$.
- Objetivo $\displaystyle E=\sum_{k=1}^{N}(x_k-x_k^d)^2+\frac{\alpha}{2}\sum_{k=0}^{N-1}m_k^2$.
- Lagrangiano y derivación de la **recurrencia adjunta hacia atrás en el tiempo**:
  $$\lambda_N=-2(x_N-x_N^d),\qquad \lambda_k=a\,\lambda_{k+1}-2(x_k-x_k^d),$$
  y el gradiente $\nabla_{m_k}E=\alpha m_k-b\,\lambda_{k+1}$.
- **Script 1** (`discrete_control`): resuelve directo/adjunto, **verifica el gradiente
  contra diferencias finitas** (error relativo $\sim10^{-9}$) y corre descenso por gradiente.
  Resultado validado: $E$ baja de $19.25$ a $\sim3.1\times10^{-3}$, error de seguimiento
  $\max_k|x_k-x_k^d|\approx3.8\times10^{-3}$.

### Bloque 3 (§3) — Poisson: discretizar primero, adjuntar, interpretar
- Discretización de $-u''=m$, $u(0)=u(1)=0$: sistema tridiagonal $Au=m$ (stencil de 3
  puntos) y funcional $J_h$ por cuadratura.
- Lagrangiano discreto con el **producto interno con matriz de masa**
  $\langle a,b\rangle_h=h\,a^Tb$.
- Ecuación adjunta $A\lambda=-2(u-u_d)$ y gradiente $\nabla_m\hat J=h(\alpha m-\lambda)$.
- **Interpretación:** $A\lambda=-2(u-u_d)$ es la discretización (mismo stencil) de la EDO
  adjunta $-\lambda''=-2(u-u_d)$ con $\lambda(0)=\lambda(1)=0$; se identifica la ecuación.
- **Script 2** (`poisson`): verificación contra diferencias finitas (error relativo
  $\sim10^{-7}$) y descenso con búsqueda de línea. Validado: $N=120$, $\alpha=10^{-3}$,
  $\max|u-u_d|\approx4.6\times10^{-2}$.

### Bloque 4 (§4) — Helmholtz: discretizar primero, adjuntar, interpretar
Ecuación de Helmholtz 2D con **incidencia de onda plana** $u^{\rm inc}=e^{ik_0 x}$ ($+x$,
convención $e^{-i\omega t}$) y capa absorbente; el campo dispersado resuelve
$(\Delta+k_0^2\varepsilon)u^s=-k_0^2(\varepsilon-1)u^{\rm inc}$, $u=u^{\rm inc}+u^s$.
- **Caso A — sísmica / scattering inverso (FWI):** discrepancia con datos
  $J=\frac12\sum_s\sum_r|u_s(r)-d_s(r)|^2$. Gradiente adjunto (convención
  $e^{-i\omega t}$, $A=-\Delta_h-k_0^2\,\mathrm{diag}(\varepsilon)-i\,\mathrm{diag}(\sigma)$)
  $\nabla_\varepsilon J=+k_0^2\sum_s\mathrm{Re}(\,\overline{\lambda_s}\odot u_s\,)$,
  con $A^H\lambda_s=r_s$ (residuo). Es el **kernel de migración** de FWI.
- **Caso B — focalización de una onda plana:** se maximiza la **intensidad en un punto**
  $\mathbf x_f=(14,7)$,
  $$J(\varepsilon)=|u(\mathbf x_f)|^2 .$$
- Derivación adjunta: para la onda plana, $\delta u=k_0^2A^{-1}\mathrm{diag}(u)\,
  \delta\varepsilon$ (con $u$ el campo total), $\nabla_\varepsilon J=k_0^2\mathrm{Re}[\overline{\lambda}\odot u]$
  y $A^H\lambda=g_u$ con $g_u=2\,u(\mathbf x_f)\,\mathbf e_f$ (una fuente puntual **en el
  foco**, retropropagada).
- **Interpretación:** $A^H\lambda=g_u$ es la discretización de la **ecuación de Helmholtz
  adjunta** (la misma ecuación con una fuente en el foco y absorción conjugada): el campo
  adjunto es el campo **retropropagado** (fase-conjugado) desde el foco.
- **Regularización:** maximizar $J-\frac{\mu}{2}\int|\nabla\varepsilon|^2$; en la práctica,
  suavizado gaussiano del gradiente (ancho $s$) y cotas $\varepsilon\in[0.25,4]$.
- **Script 3** (`focalizacion`): malla $161\times161$, $L=16\lambda$, onda plana, lente (losa
  $x\in[5,8],\,y\in[3,11]$) y foco en $(14,7)$. Verificación del gradiente contra
  diferencias finitas (error relativo $\sim10^{-7}$) y ascenso con búsqueda de línea y
  suavizado ($s=2$). Validado: $J$ de $1$ a $\approx17$ ($|u|$ en el foco $\approx4.2$) en
  $50$ iteraciones ($\sim3$–$4$ min), con lente suave.
- **Película (§4.6):** se guardan los fotogramas de cada iteración y se genera
  `focalizacion_iter.mp4` (3 paneles: $\varepsilon$, $|u|$ y $J$), embebida en el notebook
  mediante `Video(..., embed=True)`, para ver la formación de la lente y del foco.
- El óptimo **no es único** (problema no convexo y subdeterminado); la regularización
  **selecciona** una solución lisa y fabricable. Con $s=2$ el diseño quedó suave; con menos
  suavizado aparecen estructuras de alto contraste.

### Bloque 5 (§5) — El punto de vista en dimensión infinita
- Teorema de Lagrange en espacios de Hilbert y el **operador adjunto** $(D_zG)^{*}$.
- **Derivada de Gateaux** y ejemplos ($\nabla_u J_1=2(u-u_d)$, $\nabla_m J_2=\alpha m$).
- Ecuación adjunta continua por integración por partes para Poisson: se recupera (3.7).
- Discusión **DTO vs. OTD** (discretizar-y-optimizar vs. optimizar-y-discretizar): en general
  los gradientes no coinciden exactamente (diferencias de cuadratura, condiciones de borde y
  del operador); en Poisson coinciden exactamente gracias a la matriz de masa. En la práctica
  se usa DTO.

### Bloque 6 — Cierre
- Comparación de costos: $p$ resoluciones (sensibilidad) vs. directo + adjunto.
- Ejercicios propuestos (incluye DTO vs. OTD).

## 3. Decisiones técnicas clave

| Tema | Decisión | Motivo |
|---|---|---|
| Formulación del control discreto | recurrencia temporal explícita | Muestra el adjunto como “marcha atrás en el tiempo”. |
| EDPs | **discretizar primero** y usar el producto interno con matriz de masa ($h\,a^Tb$) | Permite interpretar el adjunto discreto como la discretización del adjunto continuo. |
| Orden de la presentación | dimensión finita primero; dimensión infinita (Gateaux/Hilbert/OTD) al final | El gradiente discreto es el exacto del funcional; la interpretación continua cierra el círculo. |
| Verificación | diferencias finitas centradas en nodos aleatorios | Prueba rigurosa y barata del gradiente adjunto. |
| Helmholtz | diferencias finitas + capa absorbente (imaginaria) y `scipy.sparse.linalg.splu` | Simple, rápido y suficiente para la demostración. |
| Objetivo de diseño | **intensidad en el foco** $J=|u(\mathbf x_f)|^2$ | Muy visual: la onda plana converge a un punto brillante. |
| Incidencia | **onda plana** desde la izquierda ($+x$), formulación campo total/dispersado | Frente colimado limpio que la lente focaliza en $(14,7)$. |
| Regularización | suavizado gaussiano del gradiente (ancho $s$) y cotas de $\varepsilon$ | El diseño inverso es subdeterminado; permite el compromiso focalización/suavidad. |
| Optimizador | descenso/ascente por gradiente con backtracking + suavizado gaussiano | Didáctico (el objetivo de la unidad es el gradiente) y estable. |
| Reutilización | matriz complejo-simétrica $\Rightarrow A^H=$ conjugada; se factorizan ambas una vez por iteración | Muestra el ahorro del adjunto. |

## 4. Archivos

- `metodo_adjunto.md` (este documento).
- `unidad-5.ipynb`: notebook final, construido con `nbformat`, con las celdas del bloque 0–6.
- Scripts de prototipo en `/tmp/opencode/proto_*.py` (validación, no se versionan).

## 5. Criterios de aceptación

1. El notebook es JSON válido y abre como notebook de Jupyter/MyST.
2. Las tres verificaciones de gradiente contra diferencias finitas dan error relativo
   $<10^{-5}$.
3. Los scripts corren de punta a punta sin errores.
4. La narrativa sigue el hilo **dimensión finita primero**: motivación → Lagrange en
   $\mathbb R^n$ → control discreto → Poisson (discretizar, adjuntar, interpretar) →
   Helmholtz (idem) → dimensión infinita (Hilbert, Gateaux, DTO vs. OTD).
5. La funcional de la EDP de focalización queda derivada e implementada, con regularización,
   visualización del foco y película de las iteraciones.
6. Cada EDP incluye la **interpretación** del adjunto discreto como discretización de la
   ecuación diferencial adjunta correspondiente.
