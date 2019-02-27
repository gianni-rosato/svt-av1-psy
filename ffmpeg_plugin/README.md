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
If no SVT-HEVC plugin in FFmpeg
- git apply ../SVT-AV1/ffmpeg_plugin/0001-Add-ability-for-ffmpeg-to-run-svt-av1.patch
- ./configure --enable-libsvtav1
Or if based on the SVT-HEVC plugin patch in FFmpeg
- git apply ../SVT-AV1/ffmpeg_plugin/0001-Add-ability-for-ffmpeg-to-run-svt-av1-with-svt-hevc.patch
- ./configure --enable-libsvthevc (optional) --enable-libsvtav1
- make -j `nproc`

3. Verify
>> ffmpeg is now built with svt-av1, sample command line: 
./ffmpeg  -i input.mp4 -c:v libsvt_av1 -rc 1 -b:v 10M -preset 1  -y test.265
./ffmpeg  -i input.mp4 -vframes 1000 -c:v libsvt_av1 -y test.mp4

