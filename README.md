## 📐 Locally-Exact Homogenization Theory (LEHT)

### 
The Locally-Exact Homogenization Theory (LEHT) is an analytical approach based on the Trefftz concept, in which the local fields are represented by series expansions that satisfy the governing differential equations. The solution is obtained by imposing continuity conditions at the fiber–matrix interface and periodicity conditions on the unit cell. This methodology allows the effective thermal conductivity of materials with inclusions to be determined. The LEHT formulation implemented in this repository is based on the concepts presented in **DOI:** [https://doi.org/10.1016/j.ijheatmasstransfer.2020.119477].

### Syntax

The `LEHT_Thermal.m` function computes the effective thermal conductivity matrix, generates the 2D temperature field, and extracts 1D temperature profiles for a composite material considering a circular inclusion within a square matrix.
* LEHT_Thermal(k_m, k_i, frac, field, x_cut, y_cut)


**Table 1:** Inputs parameters' declaration - LEHT
---
| Parameter | Description | Accepted Values |
| :---: | :--- | :---: |
| **`k_m`** | Thermal conductivity of the matrix phase. | `> 0` |
| **`k_i`** | Thermal conductivity of the inclusion phase. | `> 0` |
| **`frac`** | Volume fraction of the inclusion in the Representative Unit Cell (RUC). | `[0.05, 0.75]` |
| **`field`** | Enables or disables the plotting of the total 2D temperature field. | `0` (disable) or `1` (enable) |
| **`x_cut`** | Coordinate to extract the vertical temperature profile. | `0` (disable) or `0 < x_cut <= 1` |
| **`y_cut`** | Coordinate to extract the horizontal temperature profile. | `0` (disable) or `0 < y_cut <= 1` |


### Usage Example

To run the analysis with a matrix conductivity of $0.5 \ W/(m \cdot °C)$, inclusion conductivity of $4.5 \ W/(m \cdot °C)$, and a volume fraction of 60 %, while also generating the 2D temperature field and extracting profiles at $x_1 = 0.25$ and $x_2 = 0.55$, execute the following command:
* LEHT(0.5, 4.5, 0.6, 1, 0.25, 0.55)

***Command Window Output:***
```text
EFFECTIVE THERMAL CONDUCTIVITY MATRIX (K*)
    1.4722    0.0000
   -0.0000    1.4722
```

***Graphical Results:***

The command above will also generate the following plots:

| <img src="imagens/fieldLEHT.png" width="500"> | <img src="imagens/tempV.png" width="500"> | <img src="imagens/tempH.png" width="500"> |
| :---: | :---:  |:---: |

---

