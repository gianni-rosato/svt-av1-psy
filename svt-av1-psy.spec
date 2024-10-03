
Name:           SVT-AV1-PSY
Version:        2.2.1-A
Release:        1%{?dist}
Summary:        The Scalable Video Technology for AV1 (SVT-AV1 Encoder and Decoder) with perceptual enhancements for psychovisually optimal AV1 encoding \
License:        (BSD-3-Clause-Clear) AND (BSD-2-Clause)
URL:            https://github.com/gianni-rosato/svt-av1-psy
Source0:     https://github.com/gianni-rosato/svt-av1-psy/archive/refs/tags/v%{version}.tar.gz

BuildRequires:  cmake
BuildRequires:  yasm-devel
BuildRequires:  clang
BuildRequires:  clang++
BuildRequires: libdovi-devel

%description
SVT-AV1-PSY is the Scalable Video Technology for AV1 (SVT-AV1 Encoder and Decoder) with perceptual enhancements for psychovisually optimal AV1 encoding.


%build
%cmake -DLIBDOVI_FOUND=1 
%cmake_build

%install
%cmake_install

%check
%ctest

%files
%license LICENSE
%{_bindir}/%{name}

%changelog
* Thu Sep 26 2024 Gianni Rosato <giannirosato@proton.me> - 2.2.1-A
PSY Updates
Features

    Encoding with odd (non-mod2) dimensions is now possible
    Encoding at resolutions lower than 64x64 is now possible, down to as small as 4x4

Quality & Performance

    A new variance boost curve for Tune 4 has been introduced, optimized for still image encoding performance
    Improved color reproduction & overall picture quality in Tune 4 & Tune 2 through chroma qindex scaling (similar functionality is already present in Tune 3)
    Higher quality presets are available when encoding 16K
    Default --chroma-qm-min has been updated to 8 (from 0)

Documentation

    Updated Handbrake link to new, cross-platform Handbrake PSY builds

Bug Fixes

    The version number should now be correct when building from an archive instead of cloning with git
    AVIFs with dimensions exceeding 4K resolution should decode correctly in applications that previously refused to decode these images when they were produced by SVT-AV1(-PSY)
