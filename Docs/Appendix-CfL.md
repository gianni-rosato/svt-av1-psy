[Top level](../README.md)

# Chroma from Luma Prediction

## 1. Description of the algorithm

The general idea behind the chroma from luma (CfL) prediction feature is to
exploit the correlation between luma and chroma to express the Intra prediction
of chroma sample values as an affine function of the corresponding
reconstructed luma sample values, where the reconstructed luma samples are
sub-sampled to match the chroma sub-sampling. The chroma prediction is given by

$`Chroma_{pred}=\alpha*Luma_{recon}+\beta`$

where $`Chroma_{pred}`$ and
$`Luma_{recon}`$ are predicted chroma
and reconstructed luma samples, respectively. The parameters
$`\alpha`$ and
$`\beta`$ can be determined (at least
theoretically) using least squares regression. The feature provides gains in
screen sharing applications.

In practice, the CfL prediction is performed as illustrated in Figure 1 below.

![cfl_appendix_fig1](./img/cfl_appendix_fig1.png)

##### Figure 1. Block diagram of the chroma from luma prediction process.

The steps illustrated in the diagram above can be summarized as follows:

  - Consider the reconstructed luma sample values.

  - Reconstructed luma samples are sub-sampled to match the chroma
    sub-sampling.

  - Calculate the $`DC_{Luma\_recon}`$ (i.e. average) of the
    reconstructed luma sample values.

  - Subtract the $`DC_{Luma\_recon}`$ from the reconstructed luma sample values to generate the AC reconstructed luma sample values, $`Luma_{recon\_AC}`$, which has a zero average.

  - Compute $`\alpha_{AC}`$ using the AC reconstructed luma sample values.

  - Compute the intra DC mode chroma prediction, $`DC_{chroma}`$. The final chroma
    from luma prediction is then given by:

$`Chroma_{pred} = \alpha_{AC} * Luma_{recon,AC} + DC_{Chroma}`$

## 2. Implementation of the algorithm

**Inputs**: luma inverse quantized residuals

**Outputs**: Best $`\alpha`$ and chroma residuals

**Control macros/flags**:

##### Table 1. Control flags related to the CfL prediction.
| **Flag**          | **Level**     | **Description**                                                                      |
| ----------------- | ------------- | ------------------------------------------------------------------------------------ |
| cfl_level         | Picture       | Describes the CfL level of the encoder.                                              |

**Details of the implementation**

![CfL_fig1](./img/CfL_fig1.png)

##### Figure 2. High-level encoder pipeline dataflow with CfL feature.


![CfL_fig2](./img/CfL_fig2.png)

##### Figure 3. The main function calls leading to CfL prediction. The functions highlighted in blue are where CfL prediction takes place.


![CfL_fig3](./img/CfL_fig3.png)

##### Figure 4. Continuation of Figure 2 showing the details of CfL processing in the function CfLPrediction.

The high level dataflow of CfL in SVT-AV1 is shown in Figure 2. CfL prediction takes place in MD through the function ```CflPrediction```
and in the encode pass through the function ```Av1EncodeLoop/Av1EncodeLoop16bit```. The details of the CfL prediction in the function ```CflPrediction``` are presented in Figure 4.
Similar flow is also followed in the function ```Av1EncodeLoop/Av1EncodeLoop16bit```, except for the fact that $`\alpha`$
is calculated only in MD and the encode pass would use the same $`\alpha`$
to perform the final CfL prediction. In the following, the details of the CfL processing in the function ```CflPrediction``` are presented.

For an intra coded block, the function ```CflPrediction``` is called when the ```intra_chroma_mode``` is set to ```UV_CFL_PRED```. There are four steps in the function:

**Step 1**: Reconstruct the Luma samples (```AV1PerformInverseTransformReconLuma```)

The first step is to reconstruct the luma samples, since the latter would be
used to generate the chroma prediction. At this stage in the encoder pipeline,
the luma residuals are transformed, quantized and inverse quantized. In this
step, the inverse transform is applied, and the reconstructed luma residuals
are added to the prediction to build the reconstructed samples.

**Step 2**: Compute the AC component of the luma intra prediction

In this step, the luma reconstructed samples are down sampled to match
the size of chroma samples using the ``` cfl_luma_subsampling_420 ```
function. Then the AC luma values are calculated by subtracting the DC luma
value using the ```svt_subtract_average``` function. The resulting AC values are stored
in the ```pred_buf_q3 buffer```.

**Step 3**: Find the best $`\alpha`$

The best $`\alpha`$ values for the chroma components are calculated by
minimizing the overall full cost. The algorithm performs a search over the 16 possible
values of $`\alpha`$ and finds the best value that minimizes the joint prediction cost.
The search is performed in the context of a joint sign between the two chroma components.
After the best value for $`\alpha`$ is calculated, the joint cost is compared with the cost of DC prediction and the winner is selected.


**Step 4**: Generate the chroma prediction

After the best $`\alpha`$ is selected, the prediction using the
CfL mode is performed using the ```svt_cfl_predict``` function. The chroma
residuals are then calculated using the function ```residual_kernel```.

## 3. Optimization of the algorithm

Finding the best $`\alpha`$ requires searching different
values in the set of allowed $`\alpha`$ values and calculating the cost
associated with each value. Performing this $`\alpha`$ search
process in MD for every luma mode and block size
at MD would be very costly. In order to find the best quality-speed
trade offs for the feature, CfL and UV (i.e. chroma) control signals are defined with multiple levels.
Table 2 shows the CfL control signals and their descriptions.
The CfL control signals are set in the function ```set_cfl_ctrls``` based on the ```cfl_level``` value.

##### Table 2. CfL control signals description.

| **Signal**        | **Description**                                                                                                                                                                                           |
| ----------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| enabled           | 0/1: Disable/Enable CfL candidate injection                                                                                                                                                               |
| itr_th            | Threshold to indicate the minimum number of α values to try. However if a large enough number of α values are evaluated without improvements in the overall rate-distortion cost, the search would stop.  |

Table 3 shows the CfL-related UV control signal and its description. The signal is set in the function ```set_chroma_controls``` based on the chroma level ```uv_level```.

##### Table 3. CfL-related UV control signals description.

| **Signal**        | **Description**                                                                                                       |
| ----------------- | ------------------------------------------------------------------------------------------------------------          |
| uv_cfl_th         | Threshold to skip CfL if the ratio of the best intra cost to the best inter cost is greater than uv_cfl_th.           |

The CfL and UV levels are set according to the encoder preset, PD_PASS, temporal layer index, slice type and screen content class.

## 4. Signaling

CfL is an Intra chroma mode that is allowed only for blocks with height and width of 32 or smaller.
The entropy encoder signals the chroma mode per block and if the mode is CfL,
extra parameters are included in the bit stream:

  - ```cfl_alpha_signs``` contains the sign of the alpha values for U and
    V packed together into a single syntax element with 8 possible
    values. (The combination of two zero signs is prohibited as it is
    redundant with DC Intra prediction.)

## Notes

The feature settings that are described in this document were compiled at
v1.7.0 of the code and may not reflect the current status of the code. The
description in this document represents an example showing how features would
interact with the SVT architecture. For the most up-to-date settings, it's
recommended to review the section of the code implementing this feature.

## References

[1] Luc N. Trudeau, Nathan E. Egge and David Barr,
“Predicting Chroma from Luma in AV1”, Data Compression Conference, 2017.
