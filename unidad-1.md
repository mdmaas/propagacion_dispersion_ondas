# Unidad 1: La ecuación de ondas y el principio de Huygens–Fresnel

La teoría y algunos ejercicios de esta unidad se desarrollan en la [unidad_1](pdf/chap_1.pdf).

El principal resultado de esta unidad es la Fórmula de representación de Kirchhoff, que es una formalización matemática, deducible a partir de la ecuación de ondas, del principio de Huygens–Fresnel:

```{math}
:label: eq-kirchhoff
u(\mathbf{r}) =
\int_\Gamma \left[
    \Phi_k(\mathbf{r},\mathbf{r}')\,\frac{\partial u}{\partial n}(\mathbf{r}')
    - u(\mathbf{r}')\,\frac{\partial \Phi_k(\mathbf{r},\mathbf{r}')}{\partial n(\mathbf{r}')}
\right] dS(\mathbf{r}')
```

