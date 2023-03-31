FROM buildpack-deps:kinetic as BUILDER

WORKDIR /app

RUN apt-get update && apt-get install -y \
    cmake \
    yasm \
    && rm -rf /var/lib/apt/lists/*

ENV LDFLAGS "-static -static-libgcc"

RUN git clone --depth 1 https://gitlab.com/AOMediaCodec/SVT-AV1.git && \
    cd ./SVT-AV1/Build/linux && \
    ./build.sh release static



FROM alpine as RELEASE

COPY --from=BUILDER /app/SVT-AV1/Bin/Release/SvtAv1EncApp /app/

ENTRYPOINT ["/app/SvtAv1EncApp"]
CMD [ "--help" ]
