#!/bin/bash
java -jar -Xmx64m  /home/pharting/workspace/DawidColombo/Replication/BatchParallelModel.jar 1000 industryScenario 1 -0.1 strategyParameter 3.0 marketEntryExit false strategyExperiment true marketEntryNu 3.0 fractionEffectivityInnovators 0.05 fractionEffectivityImitators 0.05 rateLocationDecision 0.02 locationCosts 0.01 /home/pharting/workspace/DawidColombo/Replication >>/dev/null
