#!/bin/bash

# put this in project launch configuration properties of eclipse (or other IDE) instead of
# project executable to clean img files prior to launching the application

echo starting...


sleep 1
#rm ./images/*.png
#sleep 1

EXPERIMENT="4a"
echo
echo
echo  --- Experiment$EXPERIMENT ---
echo
cp sm.experiment$EXPERIMENT.ini sm.ini
./Debug/cli2 | tee -a $EXPERIMENT_results

EXPERIMENT="5a"
echo
echo
echo  --- Experiment$EXPERIMENT ---
echo
cp sm.experiment$EXPERIMENT.ini sm.ini
./Debug/cli2 | tee -a $EXPERIMENT_results

EXPERIMENT="6a"
echo
echo
echo  --- Experiment$EXPERIMENT ---
echo
cp sm.experiment$EXPERIMENT.ini sm.ini
./Debug/cli2 | tee -a $EXPERIMENT_results

EXPERIMENT="7a"
echo
echo
echo  --- Experiment$EXPERIMENT ---
echo
cp sm.experiment$EXPERIMENT.ini sm.ini
./Debug/cli2 | tee -a $EXPERIMENT_results

EXPERIMENT="8a"
echo
echo
echo  --- Experiment$EXPERIMENT ---
echo
cp sm.experiment$EXPERIMENT.ini sm.ini
./Debug/cli2 | tee -a $EXPERIMENT_results



#rm ./images/*.tmp
echo finished
