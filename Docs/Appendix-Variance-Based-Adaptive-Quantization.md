# Variance Based Adaptive Quantization

It is well known that the visibility of compression artifacts in a patch of  an encoded image or video largely depends on the spatial/spatiotemporal contrast of the patch region.  This process is formally known as perceptual contrast masking, and the extent of distortion masking of a spatial/spatiotemporal patch is reported by measuring the contrast masking threshold of the patch through subjective experiments. The purpose of adaptive quantization in video encoding is to exploit the property of perceptual contrast masking to adjust the quantization parameter (QP) value of a region according to its contrast masking threshold. Typically, an adaptive quantization algorithm seeks to use a lower QP value for encoding regions of low contrast where the threshold is low, so that these regions can be rendered at a better perceptual quality, while using higher QP values for high contrast (higher spatial frequency) regions which can mask distortions to a larger extent.

## 1. Description of the algorithm

In SVT-AV1, the implementation of adaptive quantization relies on block variances, which are used as gross estimates of the perceptual contrast of the spatial blocks. In particular, the QP can be adjusted at the coding unit (CU) level. The syntax for the segmentation feature in AV1 allows the frames to be divided into a maximum of 8 segments, where each CU can be assigned a segment id. For each segment, a quantization parameter, loop filter strength, prediction reference frame and skip mode can be specified.  This segmentation is at the core of the adaptive quantization algorithm in SVT-AV1 and is described next.

### Step 1: Finding segmentation parameters using variance distribution

The frame level segmentation step in SVT-AV1 utilizes variances of 8✕8, 16✕16, 32✕32 and 64✕64 blocks that are computed during the statistics/pre-processing stage, and are thus available  to the segmentation process. For each frame to be encoded, the segmentation parameters for the frame are calculated first, before any block mode decisions are made.  The variance <img src="https://render.githubusercontent.com/render/math?math=\sigma_i"> of each 64✕64 superblock  is computed as the average of the variances of all local 8✕8 blocks that are enclosed by the superblock, i.e.

<img src="https://render.githubusercontent.com/render/math?math=\sigma_i = \frac{1}{64}\sum_{j=0}^{63}\sigma_{ij}  \ \text{for} \ i=1 \cdots K">

where <img src="https://render.githubusercontent.com/render/math?math=\sigma_{ij}"> are the variances of the 8✕8 blocks within superblock *i* and *K* is the total number of superblocks in the video frame. The superblock variances thus computed are pooled across the frame, and the logarithm of their maximum, minimum and average are computed.  Let these quantities be <img src="https://render.githubusercontent.com/render/math?math=\sigma_{max}">, <img src="https://render.githubusercontent.com/render/math?math=\sigma_{min}">, and <img src="https://render.githubusercontent.com/render/math?math=\sigma_{avg}"> respectively.

Subsequently the histogram bin width is computed as <img src="https://render.githubusercontent.com/render/math?math=w=(\sigma_{max}-\sigma_{min})/N">, where *w* is the bin width and *N* is the total number of segments (*N*=8). The bin edges of the histogram are given by <img src="https://render.githubusercontent.com/render/math?math=e_i=(\sigma_{min} %2B i \cdot w)"> for  <img src="https://render.githubusercontent.com/render/math?math=i=0 \ \ \cdots \ \ N">. Each histogram bin corresponds to a segment, and for each segment, the deviation of the segment's q-index from the base q-index for the frame, denoted by Δq is assigned as follows:

<img src="https://render.githubusercontent.com/render/math?math=$\Delta q_{i %2B 1}=\frac{(e_i %2B e_{i %2B 1})} {2} - \sigma_{avg} \ \ \  \text{for} \ \ i=0 \ \cdots \ N-1">

The values <img src="https://render.githubusercontent.com/render/math?math=e_i">as well as <img src="https://render.githubusercontent.com/render/math?math=\Delta q_i"> for <img src="https://render.githubusercontent.com/render/math?math=i = 1 \cdots N"> are saved as segmentation parameters of the frame at the end of Step 1.

### Step 2: Segmentation and QP assignment to CUs

In the encode pass, once the superblock partition decisions are made, and CU sizes are available for each superblock, the quantization parameter to be used to encode each CU is individually determined using the segmentation parameters computed Step 1. This involves a table lookup to get the variance of the current CU, using the geometry of the current CU and the variance pointer for the current superblock.  Since only variances of square blocks are calculated during the pre-processing stage, the variance of a rectangular CU is computed as the average of the variances of the two adjacent square blocks of the next smaller size that make up the rectangular CU. Further, as the smallest block size for which block variances are available is 8✕8, all smaller block variances are approximated by the variance of the 8x8 block that encloses such smaller blocks.

Once the variance of the CU is computed as described above, the histogram bin corresponding to the current CU is determined.  A CU is with variance <img src="https://render.githubusercontent.com/render/math?math=\sigma_c">is assigned to segment *i*  if <img src="https://render.githubusercontent.com/render/math?math=\sigma_c \leq  e_i"> and <img src="https://render.githubusercontent.com/render/math?math=\sigma_c > e_j \ \ \forall j > i">. The Δq value for the segment to which the CU is assigned is looked up, and the q-index for the CU is obtained by adding the Δq to the base q-index of the frame. Finally, the q-index is mapped to the QP value using the appropriate look-up table, and this QP value is subsequently used to encode the CU.

## 2. Implementation of the algorithm

**Inputs**: superblock pointer, CU pointer, CU block geometry and segmentation parameters.

**Outputs**: CU segment id and QP value.

**Control macros/flags**: `adaptive-quantization` is the high level flag used to enable/disable adaptive quantization (default: disabled).

### Implementation details

The segmentation parameters are defined by the structure `SegmentationParams` in EbSegmentation.h. An array member of  `SegmentationParams`, defined by `feature_data[MAX_SEGMENTS][SEG_LVL_MAX]` contains the feature data for each segment; here, the macro `MAX_SEGMENTS` defines the number of allowed segments for each frame and is set to 8 in SVT-AV1, while `SEG_LVL_MAX` is the enum for the maximum number of features per segment. In the current implementation in SVT-AV1, the segmentation map is always updated for every frame, and thus the function *setup_segmentation()*  in EbSegmentation.c sets the corresponding segmentation parameters as follows:

```c
segmentation_params->segmentation_update_data = 1
segmentation_params->segmentation_update_map = 1;
segmentation_params->segmentation_temporal_update = EB_FALSE;
```

After computing the bin edges and bin centers of the variance histogram, the bin edge as well as the *Δq* of each segment of the frame are updated as follows in the function *find_segment_qps()* in EbSegmentation.c  :

```C
for (int i = 0; i < MAX_SEGMENTS; i++) {
        segmentation_params->variance_bin_edge[i] = POW2(bin_edge);
        segmentation_params->feature_data[i][SEG_LVL_ALT_Q] =
            ROUND((uint16_t)strength * (MAX(1, bin_center) - avg_var));
        bin_edge += step_size;
        bin_center += step_size;
    }
```

Here, `SEG_LVL_ALT_Q` is the enum for the position corresponding to *Δq* in the `feature_data` data structure. The logarithm operation previously used to compute the variance bin edges is inverted by an exponential operation when the bin edges are saved as segmentation parameters.

In the encoding phase, the variance of each CU is found through a table lookup by the function *get_variance_for_cu()* and subsequently adaptive quantization takes place in the function *apply_segmentation_based_quantization()* by finding the histogram bin corresponding to the CU variance and assigning the corresponding segment id and QP value to the CU as follows:

```C
uint16_t variance = get_variance_for_cu(blk_geom, variance_ptr);
    for (int i = 0; i < MAX_SEGMENTS; i++) {
        if (Log2f(variance) <= segmentation_params->variance_bin_edge[i]) {
            cu_ptr->segment_id = i;
            break;
        }
    }
    int32_t q_index = picture_control_set_ptr->parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx +
                      picture_control_set_ptr->parent_pcs_ptr->frm_hdr.segmentation_params.feature_data[cu_ptr->segment_id][SEG_LVL_ALT_Q];
    cu_ptr->qp = q_index_to_quantizer[q_index];
```

### Relevant files and functions in the codebase

The main source files that implement adaptive quantization are located in Source/Lib/Encoder/Codec and are:

- EbSegmentation.c
- EbSegmentation.h

Additionally, the following files which define the data structures for segmentation parameters and are located in Source/Lib/Common/Codec:

- EbSegmentationParams.c
- EbSegmentationParams.h

The relevant functions used in EbSegmentation.c to perform the  adaptive quantization are described in listed below:

| Function                                  | Description                                                  |
| ----------------------------------------- | ------------------------------------------------------------ |
| *find_segment_qps()*                      | computes the distribution of superblock variances over each frame and initializes the segmentation parameters corresponding to variance bin edges and segment *Δq* values. |
| *setup_segmentation()*                    | initializes the remaining segmentation parameters frame by frame. |
| *get_variance_for_cu()*                   | performs a table lookup to find the variance of the current CU while encoding. |
| *apply_segmentation_based_quantization()* | finds the histogram bin for the current CU and assigns segment id as well as the QP value for the CU. |



## 3. Optimization of the algorithm

The adaptive quantization algorithm is computationally efficient since it uses block variances as contrast descriptors, which are available from the pre-processing statistics and needs very little extra computation to compute the variance distribution and histogram bin edges. Thereby, no optimizations are performed for higher speed presets.

## 4. Signaling

The relevant flags of the frame header OBU that are associated with the segmentation required for adaptive quantization are described in the following table:

| Flag                                                         | Description                                                  |
| :----------------------------------------------------------- | ------------------------------------------------------------ |
| `segmentation_enabled`                                       | 1: indicates segmentation is enabled for the frame (always 1 for adaptive quantization).  <br/>0: Indicates that the frame does not use segmentation. |
| `segmentation_update_map`                                    | 1: Indicates that the segmentation map are updated during the decoding of this frame.          <br/>0: Indicates that the segmentation map from the previous frame is used. |
| `segmentation_temporal_update`                               | 1: Indicates that the updates to the segmentation map are coded< relative to the existing segmentation map.<br/>0: Indicates that the new segmentation map is coded without reference to the existing segmentation map. |
| `segmentation_update_data`                                   | 1: Indicates that new parameters are about to be specified for each segment. <br/>0: Indicates that the segmentation parameters should keep their existing values. |
| `feature_enabled[i][j]` <br/> 0 ≤ i < `MAX_SEGMENTS` <br/> 0 ≤ j < `SEG_LVL_MAX` | 0: Indicates the feature j disabled for segment i. <br/>1: Indicates the feature j disabled for segment i. |
| `feature_data[i][j]` <br/> 0 ≤ i < `MAX_SEGMENTS` <br/> 0 ≤ j < `SEG_LVL_MAX` | Specifies the j<sup>th</sup> feature data for segment i.     |
