#!/usr/bin/env bash

if [[ $NOTROOT != 1 ]]; then
  SUDO=sudo
fi

[ -d "$HOME/.ccache" ] && $SUDO chown -R $USER: $HOME/.ccache
export parent_dir=$PWD
mkdir -p $parent_dir/Build/linux/${build_type:=Release}
cmake -S $parent_dir -B $parent_dir/Build/linux/$build_type -G"${generator:-Unix Makefiles}" -DCMAKE_BUILD_TYPE=$build_type -DBUILD_SHARED_LIBS=${shared_libs:-ON} -DBUILD_TESTING=${testing:-OFF} ${CMAKE_EFLAGS}
cmake -j$(nproc) --build $parent_dir/Build/linux/$build_type
$SUDO cmake --build $parent_dir/Build/linux/$build_type --target install
$parent_dir/Bin/$build_type/SvtAv1EncApp -help
