#!/bin/bash

BIN_DIR=/home/ameliajb/workarea/CxAODframework/stack_hists_signal
SAMPLE_DIR=/home/ameliajb/workarea/outputs/forPlots_v2 

mkdir stacked_plots
rm stacked_plots/*
cd stacked_plots
#ttbarall nonall, singletop, Wlv(ev/muv/taumu),Zll, Zvv ,gamma, dijet,monoW,monoH


$BIN_DIR/stack_hists all_samples Sample \
$SAMPLE_DIR/hist-monoWjjIsrDesD1m50.root monoWjjIsrDesD1m50 1 \
$SAMPLE_DIR/hist-monoWjjIsrDesD1m1300.root monoWjjIsrDesD1m1300 0 \
$SAMPLE_DIR/hist-monoWjjIsrDesD5m50.root monoWjjIsrDesD5m50 0 \
$SAMPLE_DIR/hist-monoWjjIsrDesD5m1300.root monoWjjIsrDesD5m1300 0 \
$SAMPLE_DIR/hist-monoWjjIsrDesD9m50.root monoWjjIsrDesD9m50 0 \
$SAMPLE_DIR/hist-monoWjjIsrDesD9m1300.root monoWjjIsrDesD9m1300 0 \
$SAMPLE_DIR/hist-monoZjjIsrD1m50.root monoZjjIsrD1m50 0 \
$SAMPLE_DIR/hist-monoZjjIsrD1m1300.root monoZjjIsrD1m1300 0 \
$SAMPLE_DIR/hist-monoZjjIsrD5m50.root monoZjjIsrD5m50 0 \
$SAMPLE_DIR/hist-monoZjjIsrD5m1300.root monoZjjIsrD5m1300 0 \
$SAMPLE_DIR/hist-monoZjjIsrD9m50.root monoZjjIsrD9m50 0 \
$SAMPLE_DIR/hist-monoZjjIsrD9m1300.root monoZjjIsrD9m1300 0 \
$SAMPLE_DIR/hist-monoHbb_mx1_xdxhDh.root monoHbb_mx1_xdxhDh 0 \
$SAMPLE_DIR/hist-monoHbb_mx65_xdxhDh.root monoHbb_mx65_xdxhDh 0 \
$SAMPLE_DIR/hist-monoHbb_mx1000_xdxhDh.root monoHbb_mx1000_xdxhDh 0 \
$SAMPLE_DIR/hist-monoHbb_mx1_xgxFhDh.root monoHbb_mx1_xgxFhDh 0 \
$SAMPLE_DIR/hist-monoHbb_mx65_xgxFhDh.root monoHbb_mx65_xgxFhDh 0 \
$SAMPLE_DIR/hist-monoHbb_mx1000_xgxFhDh.root monoHbb_mx1000_xgxFhDh 0 \
$SAMPLE_DIR/hist-monoHbb_mx1_zpzp100.root monoHbb_mx1_zpzp100 0 \
$SAMPLE_DIR/hist-monoHbb_mx65_zpzp100.root monoHbb_mx65_zpzp100 0 \
$SAMPLE_DIR/hist-monoHbb_mx1000_zpzp100.root monoHbb_mx1000_zpzp100 0 \

#$SAMPLE_DIR/hist-Zvv.root Zvv 1 \


cd ..
