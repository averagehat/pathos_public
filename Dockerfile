FROM centos:centos7
MAINTAINER Michael Panciera

ARG GIT_USER
ARG GIT_TOKEN

RUN yum -y update \
    && yum -y install curl bzip2 git wget make gcc-c++ \
    && curl -sSL https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /usr/local/ \
    && rm -rf /tmp/miniconda.sh \
    && conda install -y python=2 \
    && conda update conda \
    && conda clean --all --yes \
    && rpm -e --nodeps curl bzip2 \
    && yum clean all

ENV PATH "/usr/local/bin/:$PATH"

RUN mkdir /HERE && cd /HERE && git clone https://${GIT_USER}:${GIT_TOKEN}@github.com/averagehat/pathos.git \
    && cd pathos/install \
    && sh assume-conda-install.sh /usr/local/bin/ \ 
    && cd .. \
    && sh mk_yaml.sh > test.yaml \
    && cd databases \
    && make krona 
