from tulip import *
from graphclass import *
from graphfuncs import *
from measures import *
import numpy as np
import scipy.stats as sc
import csv
import openpyxl as ox
import os.path

################Ranking functions
#Function to rank the publications by their coeff
def RankPubsByCoeff(g, coeff, limit = float('inf')):
    print coeff
    H = {n : coeff[n] for n in g.graph.getNodes() if g.graph["node_type"][n] == "publication"} 
    l = sorted(H, key=H.__getitem__, reverse=True)

    #now remove all the elements of equal, bottom value
    botval = coeff[l[len(l)-1]]
    while coeff[l[len(l)-1]] == botval:
        l = l[0:len(l)-1]
    
    if limit != float('inf'):
        return l[:limit]
    else:
        return l


#Reduction functions from publications to authors
    
def PubToAuthorSum(g, coeff, autcoeff):
    for n in g.graph.getNodes():
        laut = g.graph["authors"][n]
        for aut in laut:
            autcoeff[root.id2aut[aut]] += coeff[n]

def PubToAuthorPonderedSum(g, coeff, autcoeff):
    for n in g.graph.getNodes():
        laut = g.graph["authors"][n]
        for aut in laut:
            autcoeff[root.id2aut[aut]] += coeff[n]/len(laut)

def PubToAuthorAverage(g, coeff, autcoeff):
    npubs = {aut : 0 for aut in autcoeff} #the number of publications of each author
    for n in g.graph.getNodes():
        laut = g.graph["authors"][n]
        for aut in laut:
            autcoeff[root.id2aut[aut]] += coeff[n]
            npubs[root.id2aut[aut]] +=1
                    
            
    for aut in autcoeff:
        autcoeff[aut] = autcoeff[aut]/npubs[aut]

def PubToAuthorPonderedAverage(g, coeff, autcoeff):
    npubs = {aut : 0 for aut in autcoeff} #the number of publications of each author
    for n in g.graph.getNodes():
        laut = g.graph["authors"][n]
        for aut in laut:
            autcoeff[root.id2aut[aut]] += coeff[n]/len(laut)
            npubs[root.id2aut[aut]] +=1


            
    for aut in autcoeff:
        autcoeff[aut] = autcoeff[aut]/npubs[aut]
            
            
        
    
def ReductionPubToAuthor(g, coeff, ReductionFunc, limit = float('inf')): #Deduce the top ranked authors based on the coeffs of the publications. Need to have calculated the coeffs with a measure before using this function
    #transfer the flow from the publications to the authors
    autlistgraph = []
    for n in g.graph.getNodes():
        autlistgraph+=g.graph["authors"][n]
    autlistgraph = list(set(autlistgraph))

    autlist = [root.id2aut[n] for n in autlistgraph]#get list of NODES of authors
    
    #table which associates author to coeff
    autcoeff = {}
    for aut in autlist: #initialize it
        autcoeff[aut] = 0.0

    #now cycle through the nodes and distribute coeff
    ReductionFunc(g,coeff,autcoeff)
    
    #now rank the table 
    l = sorted(autcoeff, key=autcoeff.__getitem__, reverse=True)

    #now remove all the elements of equal, bottom value
    botval = autcoeff[l[len(l)-1]]
    while autcoeff[l[len(l)-1]] == botval:
        
        l = l[0:len(l)-1]
    
    if limit != float('inf'):
        return l[:limit]
    else:
        return l




#rank authors by coeff
def RankAuthorsByCoeff(g, coeff, limit = float('inf')):
    H = {n : coeff[n] for n in g.graph.getNodes() if g.graph["node_type"][n] == "author"} 
    l = sorted(H, key=H.__getitem__, reverse=True)

    #now remove all the elements of equal, bottom value
    botval = coeff[l[len(l)-1]]
    while coeff[l[len(l)-1]] == botval:
        l = l[0:len(l)-1]
    
    if limit != float('inf'):
        return l[:limit]
    else:
        return l

################STATISTICS ON THE COEFFS
#Todo : average, ecart type, max, min, etc.


################STATISTICS ON THE GRAPH
#todo : average deg, min max deg, etc.
def AverageNodeCitDeg(g):
    k = 0.0
    #make sure the orientation is correct
    reverseFlag = False
    if g.edgeDirection == "reversed":
        reverseFlag = True

    if reverseFlag == False:
        for n in g.graph.getNodes():
            k+=g.graph.indeg(n)
        k=k/g.nNodes() #make sure this is a double and not an integer
    else:
        for n in g.graph.getNodes():
            k+=g.graph.outdeg(n)
        k=k/g.nNodes() #make sure this is a double and not an integer
        
    return k

def MaxNodeCitDeg(g):
    k = 0
    node = None
    #make sure the orientation is correct
    reverseFlag = False
    if g.edgeDirection == "reversed":
        reverseFlag = True

    if reverseFlag == False:
        for n in g.graph.getNodes():
            if g.graph.indeg(n) > k:
                node = n
                k = g.graph.indeg(n)
    else:
        for n in g.graph.getNodes():
            if g.graph.outdeg(n) > k:
                node = n
                k = g.graph.outdeg(n)
        
    return node

def MinNodeCitDeg(g):
    k = float('inf')
    node = None
    #make sure the orientation is correct
    reverseFlag = False
    if g.edgeDirection == "reversed":
        reverseFlag = True

    if reverseFlag == False:
        for n in g.graph.getNodes():
            if g.graph.indeg(n) < k:
                node = n
                k = g.graph.indeg(n)
    else:
        for n in g.graph.getNodes():
            if g.graph.outdeg(n) < k:
                node = n
                k = g.graph.outdeg(n)

    return node

################COEFFICIENT CALCULATIONS

def KendallTauDistance(l1, l2):

    def round_to_one(f): #rounds to 1 the pathological cases due to floating point arithmetic errors
        if abs(f - 1.0) < 0.0000000001:
            return 1.0
        else:
            return f
    
    H1 = {}
    H2 = {}

    k=1
    for l in l1:
        H1[l] = k
        k+=1

    k=1
    for l in l2:
        H2[l] = k
        k+=1


    keylist = list(set(H1) & set(H2))
    #now create the two ranking lists
    ll1 = []
    ll2 = []
    
    for h in keylist:
        ll1+= [H1[h]]
        ll2+= [H2[h]]
        
    return round_to_one(sc.kendalltau(ll1,ll2)[0])


#PearsonCoeffs and SpearmanRankCoeffs give the same results
def PearsonCoeff(l1,l2): 
    #we need to first convert the list of nodes to a list of integers
    H1 = {}
    H2 = {}

    k=1
    for l in l1:
        H1[l] = k
        k+=1

    k=1
    for l in l2:
        H2[l] = k
        k+=1

    keylist = list(set(H1) & set(H2))
    #now create the two ranking lists
    ll1 = []
    ll2 = []
    
    for h in keylist:
        ll1+= [H1[h]]
        ll2+= [H2[h]]
        
    return np.corrcoef(ll1,ll2)[0,1]

def SpearmanRankCoeff(l1, l2):
    H1 = {}
    H2 = {}
    
    k=1
    for l in l1:
        H1[l] = k
        k+=1

    k=1
    for l in l2:
        H2[l] = k
        k+=1


    keylist = list(set(H1) & set(H2))

    #now create the two ranking lists
    ll1 = []
    ll2 = []
    
    for h in keylist:
        ll1+= [H1[h]]
        ll2+= [H2[h]]
        
    return sc.spearmanr(ll1,ll2)[0]


"""
def KendallTauDistance(l1, l2, depth = 0, samples = 0): #use depth to only compare N random elements, use samples to do it a certain number of times and average it out
    
    H1 = {}
    H2 = {}
    
    for i in range(len(l1)):
        H1[l1[i]] = i
    for i in range(len(l2)):
        H2[l2[i]] = i

    K = 0.0

    keylist = list(set(H1) & set(H2))
    if depth != 0: #can do better. Randomize the choice of keys?
        keylist = keylist[:depth]
    
    k=0
    for i in keylist:
        #print k
        for j in keylist:
            if (H1[i] < H1[j] and H2[i] > H2[j]) or (H1[i] > H1[j] and H2[i] < H2[j]):
                K+=1.0
        k+=1

    K = K/(len(keylist)*(len(keylist)-1))

    #TODO : implement depth and samples
    
    return K
    
def SpearmanRankDistance(l1, l2): #no need for depth or samples since it's so quick

    H1 = {}
    H2 = {}

    for i in range(len(l1)):
        H1[l1[i]] = i
    for i in range(len(l2)):
        H2[l2[i]] = i

    K = 0.0

    keylist = list(set(H1) & set(H2))

    
    k=0
    for i in keylist:
        #print k
        K+=abs(H1[i] - H2[i])
        k+=1

    K = K/(len(keylist)*(len(keylist)-1))

    return 1-K

def SpearmanCoeffDistance(l1, l2, coeff1, coeff2): #no need for depth or samples since it's so quick
    #TO TEST
    
    K = 0.0

    keylist = list(set(l1) & set(l2))

    
    k=0
    for i in keylist:
        #print k
        K+=abs(coeff1[i] - coeff2[i])
        k+=1

    K = K/(len(keylist)*(len(keylist)-1))

    return K"""



#TODO : give the top 10 list of the pub or authors who are most distant in different rankings

##########Functions to compare the measures

#CompareMeasures : 
def CompareMeasures(g, names, CorrelationFunc):
    l = [RankPubsByCoeff(g,g.graph[i]) for i in names] #get all the lists

    for i in range(len(l)):
        for j in range(i):
            print names[i]+" x "+names[j],CorrelationFunc(l[i],l[j])

def PrintTopAuthors(g, names, ReducFunc):
    l = [ReductionPubToAuthor(g,g.graph[i], ReducFunc) for i in names] 

    for i in range(len(l)):
        print names[i]
        for j in range(20):
            print j," : ", root.graph["node_property"][l[i][j]]

def CompareAuthorRankings(g, coeffname, ReducFuncs, CorrelationFunc):
    lr = [ReductionPubToAuthor(g, g.graph[coeffname], r) for r in ReducFuncs]
    for i in range(len(lr)):
        for j in range(i):
            print str(ReducFuncs[i])+" x "+str(ReducFuncs[j]),CorrelationFunc(lr[i],lr[j])

def CompareAuthorMeasures(g, names, CorrelationFunc, ReducFunc):
    l = [ReductionPubToAuthor(g,g.graph[i], ReducFunc) for i in names]

    for i in range(len(l)):
        for j in range(i):
            print names[i]+" x "+names[j],CorrelationFunc(l[i],l[j])

#####Function to export the coeffs calculated into an excel file
#This function checks first if the excel file exists. If it does it doesn't recalculate whatever is already in it. Orders the coeff by alphabetical order
def ExportCoeffs(g, names, RankingFunc, CorrelationFunc, filename):
    #we assume that all the coeffs included in the loaded excel file have already been calculated 
    def get_names(ws):
        i=2
        l=[]
        while not (ws.cell(row=i, column=1).value is None):
            l += [str(ws.cell(row=i, column=1).value)]
            i+=1
        return (i,l)
    
    #wb = ox.Workbook()
    alreadythere = []
    size = 2
    if os.path.isfile(filename+".xlsx") :
        wb = ox.load_workbook(filename+'.xlsx')
        ws = wb.active
        size,alreadythere = get_names(ws)
    else:
        wb = ox.Workbook()
        ws = wb.active

    ws.title = filename
    #get the new names to add
    toadd = list(set(names) - set(alreadythere))
    toadd = sorted(toadd, key=lambda x: (x[0].isdigit(), x[1].isdigit(), x))
    totalnames = list(set(alreadythere+names))
    print toadd

    #get all the rankings
    ranks = {n : RankingFunc(g,g.graph[n]) for n in totalnames}
    
    #fill in the names of the measures to add and calculate the coeffs with the already added names
    for n in toadd:
        print "adding "+n
        size,alreadythere = get_names(ws)
        #add the name
        ws.cell(row = 1, column = size).value = ''.join([c for c in n if (c.isupper() or c.isdigit())])
        ws.cell(row = size, column = 1).value = n

        #add the name to alreadythere
        alreadythere+=[n]

        #calculate the coeff
        for i in range(len(alreadythere)):
            ws.cell(row = size, column = 2+i).value = CorrelationFunc(ranks[n],ranks[alreadythere[i]])

    wb.save(filename+".xlsx")


#This function doesn't check if an excel file exists, it recalculates everything and overwrites. Orders the coeffs in the order given by names.   
def ExportCoeffs_ordered(g, names, RankingFunc, CorrelationFunc, filename):
    #get the workbook
    wb = ox.Workbook()
    ws = wb.active
    ws.title = filename
    
    #get all the rankings
    ranks = {n : RankingFunc(g,g.graph[n]) for n in names}
    
    #fill in the names of the measures to add and calculate the coeffs with the already added names
    k = 2
    for n in names:
        print "adding "+n
        #add the name
        ws.cell(row = 1, column = k).value = ''.join([c for c in n if (c.isupper() or c.isdigit())])
        ws.cell(row = k, column = 1).value = n

        #calculate the coeff
        for i in range(k-1):
            ws.cell(row = k, column = 2+i).value = CorrelationFunc(ranks[n],ranks[names[i]])
        k+=1

    wb.save(filename+".xlsx")
