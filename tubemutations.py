
######################################
#             imports                #
######################################
import random
#import matplotlib.pyplot as plt

######################################
#            parameters              #
######################################
random.seed()
standardDeathChance=0.1
standardReproductionChance=0.1
mutatedDeathChance=0.1
mutatedReproductionChance=0.11
killAge=20
stepnum=30
tubenormnum=10
tubemutnum=10
simnum = 10

######################################
#           tube object              #
######################################
class Tube:
    def __init__(self, mutation):
        self.alive = True
        self.age = 0
        self.mut = mutation
    def __repr__(self):
        return "Tube {"+" Alive: "+str(self.alive)+" Age: "+str(self.age)+" Mutated: "+str(self.mut)+" }"
    def isAlive(self):
        return self.alive
    def getAge(self):
        return self.age
    def plusEpoch(self):
        self.age += 1
        return True
    def getMutation(self):
        return self.mut
    def kill(self):
        self.alive = False
        return True
    def reproduce(self):
        return Tube(self.mut)
     

def runSim():
    tubes = []
    for i in range(tubenormnum):
        tubes.append(Tube(False))
    for i in range(tubemutnum):
        tubes.append(Tube(True))
        
        
    for step in range(0,stepnum):
        for tube in tubes:
            if tube.isAlive() == True:
                tube.plusEpoch()
                if tube.getMutation() == True:
                    dc = mutatedDeathChance
                    rc = mutatedReproductionChance
                    if tube.getAge() == killAge:
                        tube.kill()
                else:
                    dc = standardDeathChance
                    rc = standardReproductionChance
                stat = random.choices(['broken', 'alive'], [dc, 1-dc], k=1)
                if stat[0] == 'broken':
                    tube.kill()
                division = random.choices(['reproduced', 'stayedcalm'], [rc, 1-rc], k=1)
                if division[0] == 'reproduced':
                    tubes.append(tube.reproduce())
    
    
    counter = 0    
    for tube in tubes:
        if tube.isAlive():
            #display(tube)
            counter+=1
    #display(counter) 
    
    return tubes

def stat(rt):
    mutantDomination = 0
    countmut = 0
    countstd = 0
    for tube in rt:
        if tube.isAlive() == True:
            if tube.getMutation() == True: 
                countmut += 1
            else:
                countstd += 1
    if countmut > countstd:
        mutantDomination = 1
    return mutantDomination

######################################
#                main                #
######################################
mdarr = []
countd = 0
for k in range(simnum):
    rt = runSim()
    mdarr.append(stat(rt))
    countd += stat(rt)
print(mdarr)
print(str(sum(mdarr)/len(mdarr)*100) + '% mutant domination')
