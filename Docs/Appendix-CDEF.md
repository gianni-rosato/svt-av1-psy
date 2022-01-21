# Constrained Directional Enhancement Filter (CDEF) Appendix

## 1.  Description of the algorithm

The constrained directional enhancement filter (CDEF) is applied after
the deblocking filter and aims at improving the reconstructed picture by
addressing ringing artifacts. CDEF is a combination of the directional
de-ringing filter from the Daala codec (Mozilla) and the Constrained Low
Pass Filter (CLPF) from the Thor codec (Cisco).

Filtering is applied on an 8x8 block level, which is a large enough
block size to reliably detect edges, but small enough block size to
accurately detect edge directions. Filtering is applied to both luma and
chroma samples. For a given block, the algorithm consists of the two
main steps outlined below:

1. Identify the direction **d** of the block (i.e. direction of edges). Eight directions {0,…,7} could be identified.

2. Filtering
    * Apply a nonlinear filter along the edge in the identified direction. Filter taps are aligned in the direction of the block. The main goal is to address ringing artifacts.

    * Filter mildly along a 45 degree direction from the edge.

The two steps are outlined in more detail in the following.

***Step 1 – Identification of edge direction***. Eight edge directions
could be considered. The directions are indexed with d=0,…,7 as
indicated in Figure 1 below.

![image24](./img/image24.png)

##### Figure 1. Block directions.

For a given input block, the identification of edge direction is
performed as follows:

  - Input pixels along each of the lines k=0,1,… are averaged. The
    original sample values along each on the lines k=0,1,… are replaced
    by the average value of the samples on the line. The resulting block
    is referred to as the averaged block.

  - The variance of the error between the source block and the averaged
    block is computed.

  - The operations above are repeated for each of the eight directions.

  - The direction with the lowest variance is selected as the filtering
    direction.

The example in Figure 2 below illustrates this step for an 8x8 input block. In this
example, direction 0 results in the smallest error variance and hence
was selected as the block direction.

![image25](./img/image25.png)

##### Figure 2. Example of block direction identification.

***Step 2 – Filtering***. The filtering operation consists of two main
steps, namely a primary filtering operation and a secondary filtering
operation. The primary filter acts along the identified block direction.
The secondary filter acts at ![math](http://latex.codecogs.com/gif.latex?45^o) from the identified
direction. In the example shown in Figure 3 below, the block direction
is d=0 ![math](http://latex.codecogs.com/gif.latex?45^o). The sample to be filtered is highlighted in
red. The samples to be considered when filtering the red sample in
primary filtering are highlighted in green (a total of four samples).
Those considered in the secondary filtering of the red sample are
located at ![math](http://latex.codecogs.com/gif.latex?45^o) angle from the block direction and are
highlighted in blue (a total of eight samples).

![image26](./img/image26.png)

##### Figure 3. Example of primary and secondary filtering directions.

The nonlinear low-pass filters used in the filtering process are given
by [\[1\]](#ref-1):

![math](http://latex.codecogs.com/gif.latex?P_{filtered}(i,j)=p(i,j)+\sum_{m,n}w_{m,n}f(p(m,n)-p(i,j),S,D))

Where ![image](http://latex.codecogs.com/gif.latex?p(i,j)) is the sample to be filtered,
![image](http://latex.codecogs.com/gif.latex?P_{filtered}(i,j)) is the filtered value of
sample ![image](http://latex.codecogs.com/gif.latex?p(i,j)), ![image](http://latex.codecogs.com/gif.latex?w) are the filter weights,
![image](http://latex.codecogs.com/gif.latex?f) is a nonlinear constraint function, ![image](http://latex.codecogs.com/gif.latex?S) is the
filter strength and ![image](http://latex.codecogs.com/gif.latex?D) is the filter damping.

When ![image](http://latex.codecogs.com/gif.latex?|p(m,n)-p(i,j)|) is small, ![image](http://latex.codecogs.com/gif.latex?f(p(m,n)-p(i,j),S,D)=p(m,n)-p(i,j)), implying that the filter
behaves as an FIR filter. When ![image](http://latex.codecogs.com/gif.latex?|p(m,n)-p(i,j)|) is large,
![image](http://latex.codecogs.com/gif.latex?f(p(m,n)-p(i,j),S,D)=0) and no filtering is applied to the sample. The function ![image](http://latex.codecogs.com/gif.latex?f)
de-emphasizes neighboring pixels ![image](http://latex.codecogs.com/gif.latex?p(m,n)) with large
contrast (i.e. large magnitude of ![image](http://latex.codecogs.com/gif.latex?p(m,n)-p(i,j))). The strength
![image](http://latex.codecogs.com/gif.latex?S) and damping ![image](http://latex.codecogs.com/gif.latex?D) parameters control the attenuation of the large magnitude differences.


Filtering along the identified block direction is referred to as primary
filtering, and makes use of primary filtering strength and primary
damping values. The filter weights for primary filtering are shown in
Figure 4 below. The sample to be filtered is shown in blue. For even
strengths, a = 2 and b = 4, whereas for odd strengths a = 3 and b = 3.

![image27](./img/image27.png)

##### Figure 4. Filter weights for primary filtering.

Filtering of samples that are at 45 degrees from the identified edge
direction for the block is referred to as secondary filtering. Secondary
filtering makes use of secondary strength and secondary damping values.
The filter weights for secondary filtering are indicated in Figure 5
below as a function of the block direction. The sample to be filtered is
shown in blue.

![image28](./img/image28.png)

##### Figure 5. Filter weights for secondary filtering.


## 2.  Implementation

**Inputs to cdef\_kernel**: Output frame from the deblocking filter.

**Outputs of cdef\_kernel**: CDEF filtered frame, filter parameters.

**Controlling flags**:

Control flags associated with CDEF are listed in Table 1 below.

##### Table 1. Control flags for CDEF.

| **Flag**                        | **Level**      | **Description**                                                                                                            |
| ------------------------------- | -------------- | -------------------------------------------------------------------------------------------------------------------------- |
| --enable-cdef                     | Configuration  | CDEF filter control (0:OFF , 1: ON (Default))                           |
| cdef\_level                     | Sequence       | Indicates whether to use CDEF for the whole sequence.                                                                      |
| cdef\_level                     | Picture        | Indicates the level of complexity of the CDEF strength search as a function of the encoder mode (enc\_mode).               |

**Implementation details**

Important function calls associated with CDEF are highlighted in Figure 6 below. The function calls are
organized according to the depth of the function call.

![image29](./img/image29.png)

##### Figure 6. The main function calls associated with CDEF.

The main steps involved in the implementation of the algorithm are
outlined below, followed by more details on some of the important
functions.

Step 1 - Splitting the frame into segments

- The frame to be filtered is divided into segments to allow for parallel filtering operations on
  different parts of the frame. The segments are set according to the following
  (see ```load_default_buffer_configuration_settings``` in ```EbEncHandle.c```)
  The number of segment rows is set to 1 if ```(luma height/64)<6```, else it set to 6.
- The number of segment columns is set to 1 if ```(luma width/64)<10```, else it set to 6.

The segments are processed in ```cdef_kernel```. Each segment is split into
64x64 filter blocks.

Step 2: Perform CDEF search for each segment \[each running on a
separate thread\]. Each segment goes through a filter search operation
through the function (```cdef_seg_search```). For a given 64x64 filter block
in a segment, the main purpose of the search is to identify the
directions of the 8x8 blocks in the filter block, and the best filter
(Primary strength, Secondary strength) pair to use in filtering the
filter block. The primary filter strength takes value in {0,…,15},
whereas the secondary filter strength takes value in {0, 1, 2, 4}. The
(primary strength, secondary strength) pairs are then indexed and
ordered as indicated in Table 2 below:

##### Table 2. (primary strength, secondary strength) pairs.

| **Filter Strength Index** | **(Primary Strength, Secondary Strength) Pair** |
| ------------------------- | ----------------------------------------------- |
| 0                         | (0,0)                                           |
| 1                         | (0,1)                                           |
| 2                         | (0,2)                                           |
| 3                         | (0,4)                                           |
| 4                         | (1,0)                                           |
| 5                         | (1,1)                                           |
| …                         | (...,…)                                         |
| 63                        | (15,4)                                          |

The search for the best (Primary strength, Secondary strength) pair to
use is equivalent to the search for the index for such pair.

The primary luma damping (```pri_damping```) and secondary luma damping (```sec_damping```)
values are set as a function of the base qindex for the picture and are given by:

```c
pri_damping = 3 + (base_qindex/64);
sec_damping = 3 + (picture_control_set_ptr->parent_pcs_ptr->base_qindex/64);
```

Chroma damping values are always one less the luma damping value.

The CDEF search (```cdef_seg_search```) proceeds along the following steps.

  - Loop over all 64x64 filter blocks in the segment.

  - Loop over the picture planes

  - Loop over the specified filter strengths (see “Reducing Number of Filter Strengths Tested” in
    the Optimization section for more info on how to set which strengths to test.)

    - Perform the following for each 8x8 non-skip block (```cdef_filter_fb```):

      - Perform the following for each 8x8 non-skip block (```cdef_filter_fb```):

        - Find the direction for each 8x8 block (```cdef_find_dir```).
        - Filter the 8x8 block according to the identified direction using the set filter strength
          (```cdef_filter_block```, C only version ```cdef_filter_block_c```.
          More details on ```cdef_filter_block_c``` are provided below.).


    - Compute the filtering mse for the filter block corresponding to the filter strength being considered
      (```compute_cdef_dist```).

Step 3: Select a subset of filter strengths to use in the final filtering of the 64x64
filter blocks in the frame based on the filtering results from step 2
(```finish_cdef_search```. More details on ```finish_cdef_search``` are provided below.).
This step is frame-based and is performed by only one thread.

Step 4: Complete the filtering of the frame based on the selected set of filtering strengths from Step 3.
(```av1_cdef_frame```. More details on ```av1_cdef_frame``` are provided below.)

**More details about cdef\_filter\_block\_c**

For a given 8x8 block, filtering is applied to all samples in the 8x8
block. Filtering is to be applied according to the identified direction
for the 8x8 block. For a given sample to be filtered in the block, the
position of the neighboring samples to be considered in the filtering
operation are given by the array `cdef_directions` according to the
identified direction.

The primary and secondary filter coefficients are given in Section 1.

**More details on finish\_cdef\_search in step 3**

For each 64x64 filter block, the output from Step 2 is an array of distortion values
corresponding to different filter strength pairs (Primary strength, Secondary strength).
To reduce the overhead associated with the signaling of the individual filter strength index for each
64x64 filter block, only a subset of the identified filter strength pairs is selected. Final filtering of the
64x64 filter blocks in the frame is then redone using the best among the selected subset of filter strengths.
The encoder needs to signal to the decoder only the selected subset of filter strengths for the decoder to use
in the filtering operation. The encoder could signal a set that consists of only 1, 2, 4, or 8 different
(Primary strength, Secondary strength) pairs to be used for the frame. The specific pair to use for each 64x64
filter block is signaled separately. The search performed in ```finish_cdef_search``` is to find the best RDO option
(i.e. 1, 2, 4, or 8 filter strength pairs for the frame) to work with.

  - Loop over the cardinality of the set of the strength pair options (1
    then 2 then 4 then 8)

<!-- end list -->

  - Call the function ```joint_strength_search_dual``` to determine the best such set for each of the options
    based on filtering distortion (```joint_strength_search_dual``` function makes use of a greedy search
    algorithm). Compute the RDO cost of each of the options and keep track of the best option (i.e. the best
    number of bits and the corresponding set of best (Primary strength, Secondary strength) pairs. The latter are
    stored in the ```cdef_strengths```  and ```cdef_uv_strengths```  arrays.

  - Loop over the filter blocks in the frame and select for each filter block the best (Primary strength,
    Secondary strength) pair. The selected pair is signaled in ```mbmi.cdef_strength``` whereas damping values
    are stored in ```cdef_pri_damping``` and ```cdef_sec_damping```.

  - The most used filter strength pair in the filter blocks for the frame is then identified as the frame strength
    and its corresponding index is stored in ```pPcs->cdef_frame_strength```.

**More details on av1\_cdef\_frame**

```
Loop over the 64x64 filter blocks
   Loop over the three picture planes
    Call cdef_filter_fb to filter the samples in the filter block using the selected filter strength pairs
```
(selected in ```finish_cdef_search```).


## 3.  Optimization of the algorithm

The search for the best filter strength pair for each 64x64 block can be algorithmically optimized using the
features described below.  The aggressiveness of the CDEF algorithm depends on the CDEF filter mode
(```picture_control_set_ptr->cdef_level```), which is specified based on the encoder preset
(```picture_control_set_ptr->enc_mode```).

### Reducing Number of Filter Strengths Tested

The search in ```cdef_seg_search``` for the filter strength is performed by considering a subset of the allowable
filter strength indices [0,63].  For each ```cdef_level```, a set of primary and secondary filter strengths are
specified to be tested.  The search is performed in two stages:

1st stage: Test the specified primary filter strengths.
The number of primary filter strengths to test is specified by ```cdef_ctrls->first_pass_fs_num```, and
the values of the primary strengths are set in the array ```cdef_ctrls->default_first_pass_fs```.

2nd stage: Test the specified secondary filter strengths.
The number of secondary filter strengths to test is specified by ```cdef_ctrls->default_second_pass_fs_num```, and
the values of the primary strengths are set in the array ```cdef_ctrls->default_second_pass_fs```.

### Reducing Number of Rows Used in CDEF search

The CDEF search can be performed on subsampled blocks to reduce the number of required computations.  A subsampling factor is specified using
```cdef_ctrls->subsampling_factor```, according to the allowable values in Table 3.


##### Table 3. Allowable subsampling factors in CDEF search.

|**Subsampling_Factor** | **Action** |
| --------------------- | ---------- |
| 1                     | No subsampling|
| 2                     | Subsample each block by 2 (i.e. perform CDEF filtering on every 2nd row)|
| 3                     | Subsample each block by 4 (i.e. perform CDEF filtering on every 4th row)|

### Using Reference Frame Info to Reduce CDEF Search

Information from the nearest reference frames can be used to reduce the number of filter
strengths tested for each frame.  The selected CDEF filter strengths for each frame are saved in the
```EbReferenceObject``` to be used by subsequent frames.

When ```cdef_ctrls->search_best_ref_fs``` is enabled, only the best filter strengths from the reference frames are tested,
as follows:

```
If (list0_best_filter == list1_best_filter)
    Skip CDEF search; use the filter selected by the reference frames
Else if (list0_best_filter == OFF && list1_best_filter == OFF)
    Skip CDEF search; assume CDEF is off for this frame
Else {
    add list0_best_filter to be tested
    if (list1_best_filter != list0_best_filter)
        add list1_best_filter to be tested
    if (list0_best_filter_uv == OFF && list1_best_filter_uv == OFF)
        Disallow CDEF search for the chroma planes only
}
```

When ```cdef_ctrls->use_reference_cdef_fs``` is enabled, CDEF search is skipped and the filter strengths are
set to the average of the lowest and highest filter strengths of the reference frames, as follows:

```
Lowest_fs = MIN(lowest_selected_fs_from_list0, lowest_selected_fs_from_list1)
Highest_fs = MAX(highest_selected_fs_from_list0, highest_selected_fs_from_list1)
Luma_cdef_fs = MIN( 63, (lowest_fs + highest_fs) / 2)
Chroma_cdef_fs = 0
```

When ```cdef_ctrls->use_skip_detector``` is enabled, CDEF will be disabled if the skip area percentage of the
nearest reference frames (i.e. the percentage of zero coefficients in the nearest ref frames) is above 75%.

### Cost Biasing to Reduce CDEF Application

After the CDEF search, the best selected filters must be applied to each SB; however,
if the best selected filter for an SB is (0,0), then no filtering is required.
By enabling ```cdef_ctrls->zero_fs_cost_bias```, the cost of the (0,0) filter can be scaled down to make it more
favourable, resulting in fewer SBs requiring filtering in ```svt_av1_cdef_frame```.

When ```cdef_ctrls->zero_fs_cost_bias``` is non-zero, the cost of the (0,0) filter for each SB will be
scaled by (```cdef_ctrls->zero_fs_cost_bias /64```).

4.  **Signaling**

At the frame level, the algorithm signals the luma damping value and up to 8 different filter strength presets to
choose from. Each preset includes luma primary preset, chroma primary preset, luma secondary preset, a chroma
secondary preset and the number of bits used to signal the 64x64 level preset.
Table 4 summarizes the parameters signaled at the frame level.

At the 64x64 filter block level, the algorithm signals the index for the specific preset to work with for the
64x64 filter block from among the set of presets specified at the frame level.
Table 5 summarizes the parameters signaled at the filter block level.

##### Table 4. CDEF parameters signaled at the frame level.

| **Frame level Parameters**                                        | **Values (for 8-bit content)** |
| ----------------------------------------------------------------- | ------------------------------ |
| Luma Damping D                                                    | {3, 4, 5, 6}                   |
| Number of bits used for filter block signaling                    | {0,..,3}                       |
| List of 1, 2, 4 or 8 presets. Each preset contains the following: |                                |
| Luma primary strength                                             | {0,…,15}                       |
| Chroma primary strength                                           | {0,…,15}                       |
| Luma secondary strength                                           | {0,1,2,3}                      |
| Chroma secondary strength                                         | {0,1,2,3}                      |

##### Table 5. CDEF parameters signaled at the filter block level.
| **Filter-Block-level Parameters**                                 | **Values**                     |
| ----------------------------------------------------------------- | ------------------------------ |
| Index for the preset to use                                       | Up to 7                        |


## Notes

The feature settings that are described in this document were compiled at v0.9.0 of the code and may not reflect the current status of the code. The description in this document represents an example showing  how features would interact with the SVT architecture. For the most up-to-date settings, it's recommended to review the section of the code implementing this feature.


## References

<a name = "ref-1"> </a>
\[1\] Steiner Midtskogen and Jean-Marc Valin, The AV1 Constrained
Directional Enhancement Filter (CDEF), 2017.
