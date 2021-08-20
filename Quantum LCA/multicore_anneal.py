import multiprocessing
import neal
from dwave.system import LeapHybridSampler
def anneal(bqm, dwave=False, num_reads=1000, num_sweeps=1000, best=False, cores=8):
    #If best is True, then it selects any samplesets with the best energy
    processes = []
    manager = multiprocessing.Manager()
    returnDict = manager.dict()
    if(dwave):
        cores = 1
    def singleAnneal(bqm, returnDict, label):
        if(dwave):
            print("[%s] Sending to DWave quantum computer"%(label))
            sampler = LeapHybridSampler()
            sampleset = sampler.sample(bqm)
        else:
            print("[%s] Simulating annealing..."%(label))
            sampler = neal.SimulatedAnnealingSampler()
            sampleset = sampler.sample(bqm, num_reads=num_reads, num_sweeps=num_sweeps)
        returnDict[label] = sampleset
    for i in range(cores):
        p = multiprocessing.Process(target=singleAnneal, args=(bqm, returnDict, "Core %s"%(i)))
        processes.append(p)
        p.start()
    for process in processes:
        process.join()
    if(not best):
        return returnDict
    else:
        energy = list(returnDict.values())[0].first.energy
        for sampleset in returnDict.values():
            if(sampleset.first.energy < energy):
                energy = sampleset.first.energy
        bestReturnDict = {}
        for label in returnDict:
            if(returnDict[label].first.energy == energy):
                bestReturnDict[label] = returnDict[label]
        return bestReturnDict
