/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbAvailability.h"

// Bottom left availability function
EbBool  is_bottom_left_available(uint32_t depth,
    uint32_t part_index)
{
    EbBool  availability;

    uint32_t   number_of_partitions_in_width = 1 << depth;
    uint32_t   horizontal_partition_index = part_index & (number_of_partitions_in_width - 1);
    uint32_t   vertical_partition_index = (part_index - horizontal_partition_index) >> depth;

    if (depth < 4) {
        availability = (EbBool)((horizontal_partition_index    &   (number_of_partitions_in_width - 1)) == 0);                                                                    // left column border check
        availability = (EbBool)(availability | ((((horizontal_partition_index) & 0x01) == 0) && (((vertical_partition_index) & 0x01) == 0)));                           // x is even and y is even
        availability = (EbBool)(availability | ((((horizontal_partition_index) & 0x03) == 0) && (((vertical_partition_index) & 0x03) == 0x01)));                        // x is a multiple of 4  and y mod 4 is equal to 1
        availability = (EbBool)(availability &  (((vertical_partition_index != (number_of_partitions_in_width - 1)))));                                                           // bottom row border check
    }
    else {
        availability = (EbBool)((horizontal_partition_index    &   (number_of_partitions_in_width - 1)) == 0);                                                                    // left column border check
        availability = (EbBool)(availability | ((((horizontal_partition_index) & 0x01) == 0) && (((vertical_partition_index) & 0x01) == 0)));                               // x is even and y is even
        availability = (EbBool)(availability | ((((horizontal_partition_index) & 0x03) == 0) && (((vertical_partition_index + 1) & 0x03) != 0x00)));                        // x is a multiple of 4  and (y + 1) % 4 is different from 1
        availability = (EbBool)(availability | ((((horizontal_partition_index) & 0x07) == 0) && (vertical_partition_index != ((number_of_partitions_in_width >> 1) - 1))));     // x is a multiple of 8  and y is different from ((number_of_partitions_in_width) / 2 - 1)
        availability = (EbBool)(availability &  (((vertical_partition_index != (number_of_partitions_in_width - 1)))));                                                           // bottom row border check
    }

    return availability;
}

// Upper right availability function
EbBool is_upper_right_available(
    uint32_t depth,
    uint32_t part_index)
{
    EbBool  availability;

    uint32_t   number_of_partitions_in_width = 1 << depth;
    uint32_t   horizontal_partition_index = part_index & (number_of_partitions_in_width - 1);
    uint32_t   vertical_partition_index = (part_index - horizontal_partition_index) >> depth;

    if (depth < 4) {
        availability = (EbBool)(horizontal_partition_index == (number_of_partitions_in_width - 1));                                                                               // right column border check
        availability = (EbBool)(availability | ((((horizontal_partition_index) & 0x01) == 0x01) && (((vertical_partition_index) & 0x01) == 0x01)));                    // x is odd and y is odd
        availability = (EbBool)(availability | ((((horizontal_partition_index) & 0x03) == 0x03) && (((vertical_partition_index) & 0x03) == 0x02)));                    // x mod 4 is equal to 3 and y mod 4 is equal to 2
        availability = (EbBool)(availability &  (vertical_partition_index != 0));                                                                                             // top row border check
    }
    else {
        availability = (EbBool)(horizontal_partition_index == (number_of_partitions_in_width - 1));                                                                                // right column border check
        availability = (EbBool)(availability | ((((horizontal_partition_index) & 0x01) == 0x01) && (((vertical_partition_index) & 0x01) == 0x01)));                     // x is odd and y is odd
        availability = (EbBool)(availability | ((((horizontal_partition_index) & 0x03) == 0x03) && (vertical_partition_index % (number_of_partitions_in_width >> 2) != 0)));    // x mod 4 is equal to 3 and y is a multiple of (number_of_partitions_in_width / 4)
        availability = (EbBool)(availability | ((((horizontal_partition_index) & 0x07) == 0x07) && (vertical_partition_index != (number_of_partitions_in_width >> 1))));        // x mod 8 is equal to 7 and y is different from (number_of_partitions_in_width / 2)
        availability = (EbBool)(availability &  (vertical_partition_index != 0));                                                                                             // top row border check
    }

    return (EbBool)(!availability);
}
