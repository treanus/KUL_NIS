#!/bin/bash
version="v0.1 - dd 29/11/2021"

kul_main_dir=$(dirname "$0")
script=$(basename "$0")

source /usr/local/KUL_apps/KUL_NIS/KUL_main_functions.sh

rm -rf KUL_LOG
mkdir -p KUL_LOG/$script

participant="John"

task_in[0]="echo hello0; echo hello1"
task_participant[0]="Stefan"
task_in[1]="sleep 10"
task_participant[1]="Silvia"
task_in[2]="sleep 4"
task_participant[2]="Radwan"
KUL_task_exec 1 "Three tasks" "test1"

task_in="sleep 15"
unset task_participant
KUL_task_exec 1 "One tasks" "test2"
