#!/usr/bin/env bash

mcrun -c --mpi=2 BifrostFullInstrument_v1_inelastic_linewidths.instr -d v11_open_t_offset_0 -n 1000000  WaveMin=3.9  WaveMax=5.72  L_0=3.924  chopPulseOpening=0.004  DivSlit0_width=0.1  DivSlit1_width=0.1  DivSlit2_width=0.1  DivSlit3_width=0.1  Npulse=1  print=0  makeVirtualSource=0  printMValues=0  power=2  BWopen=161  sample_radius=0.005  A3=0  A4=45  Etrans=0  t_Focus_widt=1.2

