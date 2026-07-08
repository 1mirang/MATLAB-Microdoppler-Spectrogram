# MATLAB Micro-Doppler Spectrogram Analysis

This repository contains MATLAB scripts for modeling, simulating, and analyzing Micro-Doppler effects in wireless channels. 

## Project Structure
* **`patzold.m`**: Simulates the Channel Transfer Function (CTF) for a swinging pendulum and generates the theoretical Micro-Doppler spectrogram using STFT.
* **`testdata.m`**: Script for processing and visualizing real measured CSI data.

## Features
* 3D Geometrical Model Simulation.
* Moving Target Indicator (MTI) filter for static clutter removal.
* Dynamic colormap and background noise suppression.

## Reference
The algorithms and channel models implemented in this project are based on the work by:
* **A. Abdelgawwad, A. Borhani, and M. Pätzold**, *"Modelling, Analysis, and Simulation of the Micro-Doppler Effect in Wideband Indoor Channels with Confirmation Through Pendulum Experiments,"* Sensors, vol. 20, no. 4, p. 1049, Feb. 2020.
* **Link to Paper:** [https://doi.org/10.3390/s20041049](https://doi.org/10.3390/s20041049)
