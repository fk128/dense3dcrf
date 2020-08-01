FROM ubuntu:18.04

RUN apt-get update && \
    apt-get install -y build-essential git cmake autoconf libtool pkg-config wget \ 
    libboost-all-dev libtbb-dev \
    && rm -rf /var/lib/apt/lists/*


ADD https://github.com/InsightSoftwareConsortium/ITK/releases/download/v5.1.0/InsightToolkit-5.1.0.tar.gz . 
RUN tar xvf InsightToolkit-5.1.0.tar.gz && rm InsightToolkit-5.1.0.tar.gz && cd InsightToolkit-5.1.0 && mkdir build && \
    cd build && cmake -DCMAKE_INSTALL_PREFIX=../../itk  \
    -DITK_BUILD_DEFAULT_MODULES:BOOL=OFF \
    -DBUILD_EXAMPLES:BOOL=OFF \
    -DModule_ITKIONIFTI:BOOL=ON \
    -DBUILD_TESTING:BOOL=OFF --target clean  .. \
    && make -j4 && make install 
    


RUN mkdir -p /3rdparty 
ENV THIRD_PARTY_DIR=/3rdparty
WORKDIR /3rdparty
ADD https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz . 
RUN mkdir eigen && tar xf eigen-3.3.7.tar.gz -C eigen --strip-components=1 && rm eigen-3.3.7.tar.gz


COPY . /densecrf
WORKDIR /densecrf

RUN mkdir build && cd build && cmake -DCMAKE_FIND_DEBUG_MODE=ON .. && make -j4

RUN find -type f -executable -exec mv {} /usr/local/bin/ \;