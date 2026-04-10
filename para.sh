#!/bin/bash

# mode: wave or seis
export MODE=seis

export XMAX=500.0
export ZMAX=300.0
export NX=201
export NZ=201

export XSRC=250.0
export ZSRC=5.0

export XR=330.0
export ZR=150.0

export CFL=0.5
export TMAX=1.50

export F0=40.0
export T0=0.1

export ACCURACY=2

# boundary control
export USE_ABSORB=1
export W=60
export A=0.0053

export MSKIP=5
export FPS=20
export INTERVAL=50

export OUT_MOVIE=cpp_wave_animation.mp4
export OUT_SEIS=seismogram.txt

# external model files
export VS_FILE=vs_model.txt
export RHO_FILE=rho_model.txt