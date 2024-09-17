FROM fedora:39 AS base

RUN dnf -y update

# install git
RUN dnf -y install git

# install python 3.11
RUN dnf -y install python3 python3-devel gcc graphviz-devel

# install pip for python 3.11
RUN python3 -m ensurepip
RUN pip3 install --upgrade pip setuptools wheel

# Activate Jupyter extensions
RUN dnf -y --refresh install npm

# install juyter lab and other modules
RUN pip3 install \
    pytest \
    jupyter \
    jupyterlab \
    jupyter_contrib_nbextensions \
    nbconvert \
    nbformat \
    openpyxl \
    scipy \
    matplotlib\
    pandas\
    poetry\
    paho-mqtt

EXPOSE 8888

# install cimpy develop mode
RUN mkdir git
RUN cd /git && git clone https://github.com/sogno-platform/cimpy.git
RUN cd /git/cimpy/ && \
    git checkout 23-typeerror-exception-when-importing-dynamic-data-of-synchronous-machines && \
    python3 -m pip install -e .

# install pyvolt
RUN cd /git && git clone https://github.com/martinmoraga/pyvolt.git
RUN cd /git/pyvolt/ && python3 -m pip install -e .

# add DSE folder
COPY . /DSE/