%global toolchain clang

Name:           svt-av1-psy
Version:        2.2.1.A
Release:        1%{?dist}
Summary:        SVT-AV1  with perceptual enhancements  optimal AV1 encoding
License:        BSD-3-Clause-Clear AND BSD-2-Clause
URL:            https://github.com/gianni-rosato/svt-av1-psy
Source0:     https://github.com/gianni-rosato/svt-av1-psy/archive/refs/tags/v2.2.1-A.tar.gz

BuildRequires:  cmake
BuildRequires:  yasm-devel
BuildRequires:  clang
BuildRequires: libdovi-devel
BuildRequires: ninja-build

%description
NOTE: not compiled with hdr10plus
SVT-AV1-PSY is the Scalable Video Technology for AV1
(SVT-AV1 Encoder and Decoder) with perceptual enhancements for psychovisually
optimal AV1 encoding.

applications that use %{name}.

%package devel
Summary:        Development files for %{name}
Requires:       %{name}%{?_isa} = %{?epoch:%{epoch}:}%{version}-%{release}

%description devel
The %{name}-devel package contains  header files for developing


%prep
%autosetup -n %{name}-2.2.1-A

%build

%cmake -DLIBDOVI_FOUND=1 \
                 -G Ninja
%cmake_build

%install
%cmake_install

%check
%ctest

%files
%license LICENSE.md PATENTS.md
%doc CHANGELOG.md README.md Docs/
%{_libdir}/libSvtAv1Enc.so.2{,.*}
%{_bindir}/SvtAv1EncApp


%files devel
%{_includedir}/svt-av1/
%{_libdir}/pkgconfig/*.pc
%{_libdir}/libSvtAv1Enc.so

%changelog
%autochangelog
