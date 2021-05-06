# SQ_Weight Appendix

## Description

SQ\_WEIGHT determines if the evaluation of HA, HB, VA, VB, H4 and V4
shapes could be skipped based on the cost of SQ, H and V shapes.
Specifically:

-   skip HA, HB and H4 if (valid SQ and H) and (H\_COST \>
    (SQ\_WEIGHT \* SQ\_COST) / 100)

-   skip VA, VB and V4 if (valid SQ and V) and (V\_COST \>
    (SQ\_WEIGHT \* SQ\_COST) / 100)

-   `The lower the SQ_WEIGHT, the higher the chance to skip NSQ`

## SQ\_WEIGHT Derivation

SQ\_WEIGHT = **Base** + **Offset**

**Base = f(Preset)**

  **Preset(s)**  | **BASE**
  ---------------| ---------
  MR             | ∞
  M0             | 105
  M1             | 105
  M2             | 100
  M3 & beyond    | 95

**Offset = f (Target Shape, Block/PIC Type, Coeff Info, QP)**

![sq_weight_fig1](./img/sq_weight_fig1.png)

PS.
![sq_weight_fig2](./img/sq_weight_fig2.png)

## Notes

The feature settings that are described in this document were compiled at v0.8.3 of the code and may not reflect the current status of the code. The description in this document represents an example showing  how features would interact with the SVT architecture. For the most up-to-date settings, it's recommended to review the section of the code implementing this feature.
