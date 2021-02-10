FROM gitpod/workspace-full

USER gitpod

RUN sudo apt-get update \
    && sudo apt-get install -y gfortran gnuplot cmake \
    && sudo rm -rf /var/lib/apt/lists/*

# Give control back to Gitpod Layer
USER root
