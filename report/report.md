---
title: "PHSCS 513R: Final Project"
author: Eduardo Alvarez, Tyler Critchfield, Ryan Anderson
date: 9 Apr 2021
geometry: margin=1in
output: pdf_document
---

# Section Heading 

Here's a nifty equation:

$$
V(r) = A_1 \cos{r} + A_2 \cos{2r} - Fr
$$

where $A_1 = 5$, $A_2 = 1$, and $F = 1.5$ (initially)

## Subsection Heading

Now let's include a figure. Note that the path must be a relative path.

![Plot of the washboard potential over $[-10,10]$.\label{fig-potential-global}](figures/potential_global.pdf){ width=80% }

Now we can refer to Fig. \ref{fig-potential-global}.

Now, let's include some code here:

```julia
minimum value: -9.426199462493551
maximum value: 6.125688042278578
barrier:       15.551887504772129
```
