# svt-av1 ffmpeg plugin installation
## Linux

1\. Build and install SVT-AV1

``` bash
git clone --depth=1 https://github.com/OpenVisualCloud/SVT-AV1
cd SVT-AV1
cd Build
cmake .. -G"Unix Makefiles" -DCMAKE_BUILD_TYPE=Release
make -j $(nproc)
sudo make install
```

2\. Apply SVT-AV1 plugin and enable libsvtav1 to FFmpeg

``` bash
git clone -b release/4.2 --depth=1 https://github.com/FFmpeg/FFmpeg ffmpeg
cd ffmpeg
export LD_LIBRARY_PATH+=":/usr/local/lib"
export PKG_CONFIG_PATH+=":/usr/local/lib/pkgconfig"
```

- With SVT-AV1 alone:

``` bash
git apply SVT-AV1/ffmpeg_plugin/0001-Add-ability-for-ffmpeg-to-run-svt-av1.patch
./configure --enable-libsvtav1
```

- With SVT-HEVC and SVT-AV1:

``` bash
git apply SVT-HEVC/ffmpeg_plugin/0001-lavc-svt_hevc-add-libsvt-hevc-encoder-wrapper.patch
git apply SVT-AV1/ffmpeg_plugin/0001-Add-ability-for-ffmpeg-to-run-svt-av1-with-svt-hevc.patch
./configure --enable-libsvthevc --enable-libsvtav1
make -j $(nproc)
```

3\. Verify

``` bash
./ffmpeg -i input.mp4 -c:v libsvt_av1 -g 30 -vframes 1000 -y test.ivf
./ffmpeg -video_size 720x480 -pixel_format yuv420p -f rawvideo -i input.yuv -c:v libsvt_av1 -y test.mp4
```
