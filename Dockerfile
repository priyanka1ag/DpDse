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
    #matplotlib\
    pandas\
    poetry\
    paho-mqtt\
    scikit-learn\
    django-environ

EXPOSE 8888

# install cimpy develop mode
RUN mkdir git
RUN cd /git && git clone https://github.com/sogno-platform/cimpy.git
RUN cd /git/cimpy/ && python3 -m pip install -e .

# install pyvolt
RUN cd /git && git clone https://github.com/martinmoraga/pyvolt.git   
RUN cd /git/pyvolt/ && python3 -m pip install -e . 

# install platform
RUN cd /git && git clone https://github.com/SEGuRo-Projekt/Platform.git
RUN cd /git/Platform/ && python3 -m pip install -e .

# environment variable
env TLS_CACERT /DSE/keys/ca.crt
env TLS_CERT /DSE/keys/cert-acs-dse.pem   
env TLS_KEY /DSE/keys/key-acs-dse.pem     

# add DSE folder
COPY . /DSE/