# svt-av1 ffmpeg plugin installation

1. Build and install SVT-AV1 
- git clone https://github.com/OpenVisualCloud/SVT-AV1
- cd SVT-AV1
- cd Build && cmake .. && make -j `nproc` && sudo make install

2. Apply SVT-AV1 plugin and enable libsvtav1 to FFmpeg
- git clone https://github.com/FFmpeg/FFmpeg ffmpeg
- cd ffmpeg
- git checkout release/4.1
- export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/lib
- export PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:/usr/local/lib/pkgconfig
- SVT-AV1 alone:
   - git apply SVT-AV1/ffmpeg_plugin/0001-Add-ability-for-ffmpeg-to-run-svt-av1.patch
   - ./configure --enable-libsvtav1
- Based on SVT-HEVC:
   - git apply SVT-HEVC/ffmpeg_plugin/0001-lavc-svt_hevc-add-libsvt-hevc-encoder-wrapper.patch
   - git apply SVT-AV1/ffmpeg_plugin/0001-Add-ability-for-ffmpeg-to-run-svt-av1-with-svt-hevc.patch
   - ./configure --enable-libsvthevc --enable-libsvtav1
- make -j `nproc`

3. Verify
- ./ffmpeg  -i input.mp4 -c:v libsvt_av1 -g 30 -vframes 1000 -y test.ivf
- ./ffmpeg  -video_size 720x480 -pixel_format yuv420p -f rawvideo -i input.yuv -c:v libsvt_av1 -y test.mp4

