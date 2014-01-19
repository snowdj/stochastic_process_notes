---
layout: post
title: "Lecture 5: More on inference"
category: 'Lecture'
---

Instructor: Alexandre Bouchard-C&ocirc;t&eacute;

Editor: TBA


#### Moments of a DP

In this section, we derive the  first and second moments of $G(A)$, for $G \sim \dirp(\alpha\_0, G\_0)$ and  $A \subset \Omega$.  To do that, we use the Kolmogorov definition and consider the partition $(A, A^c)$ of $\Omega$.  We get:
\\begin{eqnarray}
(G(A), G(A^c)) \sim \mbox{Dir}(\alpha\_0G\_0(A), \alpha\_0G\_0(A^c)).
\\end{eqnarray}
This implies that:
\\begin{eqnarray}
G(A) \sim \mbox{Beta}(F, G),
\\end{eqnarray}
where $x$ denotes $\alpha\_0G\_0(A)$, and $y$ denotes $\alpha\_0G\_0(A^c)$.

The first moment of  $G(A)$ is therefore
\\begin{eqnarray}
\E[G(A)]=\frac{x}{x+y}=\frac{\alpha\_0G\_0(A)}{\alpha\_0G\_0(A)+\alpha\_0G\_0(A^c)}=G\_0(A), 
\\end{eqnarray}
and the second moment of  $G(A)$ is 
\\begin{eqnarray}
\var[G(A)] &=& \frac{xy}{(x+y)^2(1+x+y)} \\\\
&=&\frac{\alpha\_0^2G\_0(A)(1-G\_0(A))}{\alpha\_0^2(\alpha\_0+1)} \\\\
&=&\frac{G\_0(A)(1-G\_0(A))}{\alpha\_0+1}.
\\end{eqnarray}

This gives an interpretation of $\alpha\_0$ as a precision parameter for the Dirichlet process. 


### Supplementary references and notes

**Under construction**
