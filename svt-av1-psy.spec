
Name:           SVT-AV1-PSY
Version:        2.2.1.A
Release:        1%{?dist}
Summary:        SVT-AV1  with perceptual enhancements  optimal AV1 encoding
License:        BSD-3-Clause-Clear AND BSD-2-Clause
URL:            https://github.com/gianni-rosato/svt-av1-psy
Source0:     https://github.com/gianni-rosato/svt-av1-psy/archive/refs/tags/v2.2.1-A.tar.gz

BuildRequires:  cmake
BuildRequires:  yasm-devel
BuildRequires:  clang
BuildRequires:  clang++
BuildRequires: libdovi-devel

%description
SVT-AV1-PSY is the Scalable Video Technology for AV1
(SVT-AV1 Encoder and Decoder) with perceptual enhancements for psychovisually
optimal AV1 encoding.


%build
%cmake -DLIBDOVI_FOUND=1 
%cmake_build

%install
%cmake_install

%check
%ctest

%files
%license LICENSE PATENTS.md
%{_bindir}/%{name}
%doc CHANGELOG.md README.md Docs/*.md


%changelog
%autochangelog
