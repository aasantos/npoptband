#!/bin/bash
nvcc -arch=sm_86 -I"../include" maingpu.cu -o vix
