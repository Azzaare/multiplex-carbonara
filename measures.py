from math import *
from tulip import *
from graphclass import *
from graphfuncs import *
import networkx as nx
import scipy as sc


#Enfin peut etre une derniere categorie de flots, qui en fait seraient des flots au sein des sous graphes multiplexes, mais qui prendraient en compte le fait que certaine aretes sont dans d'autres sous graphes grace a des ponderations
#ESTCEQUIL FAUT PONDERER??

#mettre le rapport sur overleaf
#faire des mesures pour les auteurs p.ex. eigenvector networkx rendu undirected
#mesures basees sur le temps : profondeur temporelle. mesure de popularite : profondeur/tps passe pour arrive a cette profondeur


#####Classic measures of citation analysis
#H-INDEX
def HIndex_static(g, coeffname = "HIndexCoeff"): #nombre de publications ayant plus de n citations
    if g.__class__.__name__ != "PubAuthorGraph": #to change to make it compatible with other graphs
        print "must be used with a PubAuthorGraph"
        return

    autlist = [root.id2aut[n] for n in root.autlist] #get list of author nodes

    coeff = g.graph.getIntegerProperty(coeffname)
    
    #for each author we need to check the number of publications above a certain citation threshold
    w = 0
    for a in autlist: #go through the authors
        pubs = {}
        for n in g.graph.getInNodes(a):
            pubs[n] = g.graph.indeg(n) #fill in the data concerning the number of citations of the publications associated to an author
        #order it into a list from the most cited to the least cited
        lpubs = sorted(pubs, key=pubs.__getitem__, reverse=True)
        
        hindex = 0
        
        for i in lpubs:
            if pubs[i] > hindex:
                hindex+=1
            else:
                break
            
        coeff[a] = hindex
        w+=1
        print w


#A tester    
def HIndexPub_static(g, coeffname = "HIndexPubCoeff"): #TODO
    coeff = g.graph.getIntegerProperty(coeffname)
    
    #for each pubs we need to check the number of its citations above a certain citation threshold
    w = 0
    for a in g.graph.getNodes(): #go through the authors
        pubs = {}
        for n in g.graph.getInNodes(a):
            pubs[n] = g.graph.indeg(n) #fill in the data concerning the number of citations of the publications associated to an author
        #order it into a list from the most cited to the least cited
        lpubs = sorted(pubs, key=pubs.__getitem__, reverse=True)
        
        hindex = 0
        for i in lpubs:
            if pubs[i] > hindex:
                hindex+=1
            else:
                break
            
        coeff[a] = hindex
        w+=1
        print w

#1-degree
def Degree_static(g, coeffname = "DegreeCoeff"):
    k = 0
    coeff = g.graph.getDoubleProperty(coeffname)
    for n in g.graph.getNodes():
        #node_list = GetSons(g,n, direction = "reversed") #must have reversed direction
        coeff[n] = g.graph.indeg(n)#len(node_list)-1
        print k," : ",coeff[n]
        k+=1


#Degree pondered by the num of cited pubs. The same thing as 1-monoplexflow
def DegreePondered_static(g, coeffname = "DegreePonderedCoeff"):
    k = 0
    coeff = g.graph.getDoubleProperty(coeffname)
    for n in g.graph.getNodes():
        coeff[n] = 0.0
        for m in g.graph.getInNodes(n):
            coeff[n] += 1.0/g.graph.outdeg(m)
        print k, " : ", coeff[n]
        k+=1
    
        

#####Node influence metrics from the literature

#Epidemic measure
def ExForce_static(g, coeffname = "ExForceCoeff"):
    if g.edgeDirection != "direct": #eventually implement it for reversed graphs
        print "must be used with direct graph"
        return
    coeff = g.graph.getDoubleProperty(coeffname)
    k = 0
    for n in g.graph.getNodes():
        exf = 0.0
        lclusters = InfectClusters(g,n)
        dj = []
        dtot = 0.0
        for l in lclusters:
            outedges = []
            for m in l:
                for e in g.graph.getOutEdges(m): #get the list of all the outedges
                    outedges+=[e]
            outedges = list(set(outedges))#remove duplicates. Is there a more efficient way to do this? 
            dj+=[len(outedges)]
            dtot+=1.0*len(outedges)

        #now we can calculate the coeff value
        val = 0.0
        for i in dj:
            val+=-(i/dtot)*log((i/dtot))

        coeff[n] = val #and assign it
        k+=1
        print k
        
#AccessibilityMeasure
def AccessibilityMeasure_static(g):
    pass

#Probability of passage measure
def PassageProbability_static(g):
    pass


#####Centrality measures

#determine the measure of the number of nodes in the DAG associated to a node.
def KDegree_static(g, depth = float('inf'), coeffname = "KDegreeCoeff"): 
    k = 0
    coeff = g.graph.getDoubleProperty(coeffname)
    for n in g.graph.getNodes():
        node_list = SelectDag(g,n,depth, direction = "reversed") #must have reversed direction
        coeff[n] = len(node_list)-1
        print k," : ",coeff[n]
        k+=1     

#Recalculate kdegree after adding a node        
def KDegree_addNode(g, node, depth = float('inf'), coeffname = "KDegreeCoeff"):
    #this is the function which recalculates the values of the coeffs according to which node has been added
    #we assume "node" has just been added

    coeff = g.graph.getDoubleProperty(coeffname)
    node_list = SelectDag(g, node, depth, direction = "direct") #select dag with direct direction
    #for each node in node_list, increase its coeff by 1
    for n in node_list:
        coeff[n] +=1
    coeff[node] -=1 #we shouldnt add +1 to the node we have just added

#dynamic version of kdegree. Decomposes the graph with getconstructionorder and calculates the measure dynamically.     
def KDegree_dynamic(g, depth = float('inf'), coeffname= "KDegreeCoeff"):
    coeff = g.graph.getDoubleProperty(coeffname)
    coeff.SetAllNodeValue(0.0)
    gr = g.addSubGraph(EmptyGraph)
    construction_order,edges = GetConstructionOrder(g)
    #this is the function which creates an empty subgraph and slowly adds the elements, while updating the measure each time dynamically

    """#try to check if some part of the graph already exists
    starting_point = 0
    for n in construction_order:
        if g.graph["DagSizeFairCoeff"][n] !=0.0: #if the coeff has already been calculated for the node add it directly to the graph, no need to recalculate
            gr.graph.addNode(n) 
            for t in edges[n]:
                gr.graph.addEdge(n,t)
        else:
            print "starting at "+str(starting_point)
            break
        starting_point+=1
        
    construction_order = construction_order[starting_point:] #get rid of the nodes already added"""
    
    k=0
    for n in construction_order: #for all remaining nodes
        #add the node and its out edges
        gr.graph.addNode(n)
        for t in edges[n]:
            e = g.graph.existEdge(n,t)
            gr.graph.addEdge(e) 
        k+=1
        #update the flow of the graph
        KDegree_addNode(gr, n, depth)
        if k%100 == 0:
            print k
        
    #now the graph's coeffs should be uptodate. Delete the subgraph, we don't need it anymore
    g.delSubGraph(gr)


#Measure which calculates the depth of the induced dag for each node   
def DagDepth_static(g, coeffname = "DagDepthCoeff"):
    def get_depth(node):
        node_list = SelectDagDistance(g, node, depth=float('inf'), direction ="reversed")
        return max(node_list.values())

    coeff = g.graph.getDoubleProperty(coeffname)
    k=0
    for n in g.graph.getNodes():
        d = get_depth(n)
        coeff[n] = d
        print k," : ",coeff[n]
        k+=1
    
###NetworkX centrality measures

#Closeness Centrality
def Closeness_static(g, coeffname = "ClosenessCoeff"):
    n = g.addSubGraph(NetworkxGraph)
    H = nx.closeness_centrality(n.graph)
    n.transferCoeff(H, coeffname)
    g.delSubGraph(n)
    
#Betweenness Centrality
def Betweenness_static(g, coeffname = "BetweennessCoeff"):
    n = g.addSubGraph(NetworkxGraph)
    H = nx.betweenness_centrality(n.graph)
    n.transferCoeff(H, coeffname)
    g.delSubGraph(n)
    
#Katz Centrality
def Katz_static(g, coeffname = "KatzCoeff"):
    n = g.addSubGraph(NetworkxGraph)
    H = nx.katz_centrality(n.graph)
    n.transferCoeff(H, coeffname)
    g.delSubGraph(n)

#Harmonic Centrality
def Harmonic_static(g, coeffname = "HarmonicCoeff"):
    n = g.addSubGraph(NetworkxGraph)
    H = nx.harmonic_centrality(n.graph)
    n.transferCoeff(H, coeffname)
    g.delSubGraph(n)

#Indegree Centrality
def Indegree_static(g, coeffname = "IndegreeCoeff"):
    n = g.addSubGraph(NetworkxGraph)
    H = nx.in_degree_centrality(n.graph)
    n.transferCoeff(H, coeffname)
    g.delSubGraph(n)

###Information Centrality
#DOES NOT WORK
def InformationCloseness_static(g, coeffname = "InformationClosenessCoeff"):
    #equivalent a currentflow closeness centrality cf brandes and fleischer 2005
    #bon se renseigner quand meme sur la validite du truc
    n = g.addSubGraph(NetworkxGraph)
    H = information_centrality(n.graph)#nx.current_flow_closeness_centrality(n.graph)
    n.transferCoeff(H, coeffname)
    g.delSubGraph(n)

###Pondered Betweenness Centrality
#DOES NOT WORK
def InformationBetweenness_static(g, coeffname = "InformationBetweennessCoeff"):
    #se renseigner sur la validite du truc parce que bon voila quoi
    n = g.addSubGraph(NetworkxGraph)
    H = nx.current_flow_betweenness_centrality(n.graph)
    n.transferCoeff(H, coeffname)
    g.delSubGraph(n)


###Pondered Degree Centrality
#DOES NOT WORK
def InformationDegree_static(g, coeffname = "InformationDegreeCoeff"):
    #Do the same than with informatino centrality
    pass


        
#####MONOPLEX FLOW RELATED FUNCTIONS

#measure which calculates the monoplex flow        
def MonoplexFlow_static(g, coeffname="FlowCoeff"):
    def get_source_nodes(g):
	l = []	
	for n in g.graph.getNodes():
		if g.graph.indeg(n) == 0 :
			l+=[n]
	return l
    def calculate_flow(g):
        metric = g.graph[coeffname]
	if g.nNodes() == 0:	
		return 
	source_list = get_source_nodes(g)
	for n in source_list:
		for e in g.graph.getOutEdges(n):
                    m = g.graph.target(e) #get the target of the edge
                    metric[m] = metric[m]+metric[n]/g.graph.outdeg(n)
                    metric[e] = metric[n]/g.graph.outdeg(n)
		g.graph.delNode(n)	
	return calculate_flow(g)

    if g.edgeDirection != "direct" and g.edgeDirection != "reversed":
        print "Must be an oriented publication graph"
        return
    
    if tlp.AcyclicTest.isAcyclic(g.graph) == False:
        print "Graph not acyclic!"
        return
    
    gr = g.addSubGraph(CloneGraph)
            
    metric = g.graph.getDoubleProperty(coeffname)
    metric.setAllNodeValue(1.0)
    calculate_flow(gr)
    g.delSubGraph(gr)


#measure which calculates the k-monoplex flow    
def KMonoplexFlow_static(g, depth=float('inf'), coeffname="KFlowCoeff"): #resultat different probablement du aux erreurs de calcul flottants
    if depth==float('inf'):
        MonoplexFlow_static(g)
        return

    if depth==0:
        print "bug"
        return
    
    def get_source_nodes(g):
	l = []	
	for n in g.graph.getNodes():
		if g.graph.indeg(n) == 0 :
			l+=[n]
	return l

    flow_ancestry = {m : {i : 0.0 for i in range(depth)} for m in g.graph.getNodes()} #associates to each node the dict [age of flow => quantity of flow]. By convention the age of the flow coming from a neighbor is 0 and it increases by 1 at each transfer.
    
    def calculate_flow(g):
        metric = g.graph[coeffname]
	if g.nNodes() == 0:	
            return 
	source_list = get_source_nodes(g)
        print g.nNodes()

        if g.nNodes()==1:
            for k in range(depth):
                print flow_ancestry[source_list[0]][k]
        
	for n in source_list:
            for e in g.graph.getOutEdges(n):
                m = g.graph.target(e)
                flowint = 1.0/(g.graph.outdeg(n))
                flow_ancestry[m][0] += flowint #the new flow generated by the node n
                metric[m] += flowint
                metric[e] += flowint
                
                for k in range(depth-1): #now transfer the flow received in n, updating its age
                    kflow = flow_ancestry[n][k]/g.graph.outdeg(n)
                    flow_ancestry[m][k+1] = kflow
                    metric[m] += kflow
                    metric[e] += kflow
                    
                        
            g.graph.delNode(n)
            del flow_ancestry[n]
	return calculate_flow(g)

    if g.edgeDirection != "direct" and g.edgeDirection != "reversed":
        print "Must be an oriented publication graph"
        return
    
    if tlp.AcyclicTest.isAcyclic(g.graph) == False:
        print "Graph not acyclic!"
        return
    
    gr = g.addSubGraph(CloneGraph)
            
    metric = g.graph.getDoubleProperty(coeffname)
    metric.setAllNodeValue(1.0)
    metric.setAllEdgeValue(0.0)
    calculate_flow(gr)
    g.delSubGraph(gr)



        
#Recalculate the monoplex flow after adding a node
def MonoplexFlow_addNode(g, node, depth = float('inf'), coeffname="FlowCoeff"):
    #WE NEED TO SEPARATE THE CASES WHERE WE HAVE A PUBGRAPH AND A REVERSEPUBGRAPH
    #FOR THE WHILE ONLY DO PUBGRAPH, DO REVERSEPUBGRAPH LATER
    #IS THERE A MORE EFFICIENT WAY TO DO THIS??
    
    def get_source_nodes(dag):
	l = []	
	for n in dag.graph.getNodes():
		if dag.graph.indeg(n) == 0 :
			l+=[n]
	return l
    
    def calculate_flow(dag):
        metric = dag.graph["TempCoeff"]
        coeff = g.graph[coeffname]
	if dag.nNodes() == 0:	
		return
	source_list = get_source_nodes(dag)
	for n in source_list:
		for e in dag.graph.getOutEdges(n):
                    m = g.graph.target(e)
                    metric[m] = metric[m]+metric[n]/dag.graph.outdeg(n)
                    metric[e] = metric[n]/dag.graph.outdeg(n)

                coeff[n]+=metric[n] #update the coeff in the graph
		dag.graph.delNode(n)	
	return calculate_flow(dag)

    
    def transfer_flow(depth):
        dag = g.addSubGraph(InducedDagGraph,node=node, depth = depth, direction = "direct") #get the induced dag
        coeff = dag.graph.getDoubleProperty("TempCoeff")
        coeff_old = g.graph[coeffname]

        coeff.setAllNodeValue(0.0) #put all values of the dag to 0. This is so that we only calculate the flow of the added node
        coeff[node] = 1.0 #set the value of the coeff to 1 in the dag
        coeff_old[node] = 0.0 #set it to 0 in the graph
        calculate_flow(dag) #calculate the flow on the dag

        g.delSubGraph(dag)

    def transfer_flow_reversed():
        pass
    
    #we assume the node has just been added
    if g.edgeDirection == "direct":
        transfer_flow(depth)
    elif g.edgeDirection == "reversed":
        print "not yet implemented"
        return
    else:
        print "Must be an oriented publication graph"
        return

    

#Calculate the monoplex flow dynamically    
#Can add depth!
def MonoplexFlow_dynamic(g, depth = float('inf'), coeffname="FlowCoeff"):
    #Do something to check whether of the coeffs have already been calculated or not
    if g.edgeDirection != "direct" and g.edgeDirection != "reversed":
        print "Must be an oriented publication graph"
        return
    
    coeff = g.graph.getDoubleProperty(coeffname)
    gr = g.addSubGraph(EmptyGraph)
    construction_order,edges = GetConstructionOrder(g)

    k=0
    for n in construction_order: #for all remaining nodes
        #add the node and its out edges
        gr.graph.addNode(n)
        for t in edges[n]:
            e = g.graph.existEdge(n,t) #have to do this otherwise it creates duplicate edges...
            gr.graph.addEdge(e) 
        k+=1
        #update the flow of the graph
        MonoplexFlowAddNode(gr, n, depth = depth, coeffname = coeffname)
        if k%100 == 0:
            print k
        
    #now the graph's coeffs should be uptodate. Delete the subgraph, we don't need it anymore
    g.delSubGraph(gr)


####Same functions as the monoplex flow but we don't divide the flow by the outdegree when we transfer it
def MonoplexFlowUnbounded_static(g, coeffname="FlowUnboundedCoeff"):
    def get_source_nodes(g):
	l = []	
	for n in g.graph.getNodes():
		if g.graph.indeg(n) == 0 :
			l+=[n]
	return l
    def calculate_flow_unbounded(g):
        metric = g.graph[coeffname]
        if g.nNodes() == 0:	
		return 
	source_list = get_source_nodes(g)
	for n in source_list:
            for e in g.graph.getOutEdges(n):
                m = g.graph.target(e)
                metric[m] = metric[m]+metric[n]
                metric[e] = metric[n]
            g.graph.delNode(n)	
	return calculate_flow_unbounded(g)
    
    if tlp.AcyclicTest.isAcyclic(g.graph) == False:
        print "Graph not acyclic!"
        return

    if g.edgeDirection != "direct" and g.edgeDirection != "reversed":
        print "Must be an oriented publication graph"
        return

    metric = g.graph.getDoubleProperty(coeffname)
    metric.setAllNodeValue(1.0)
    
    gr = g.addSubGraph(CloneGraph)
            

    calculate_flow_unbounded(gr)
    g.delSubGraph(gr)


#TODO add the depth
def MonoplexFlowUnbounded_addNode(g, node, depth = float('inf'), coeffname="FlowUnboundedCoeff"):
    def get_source_nodes(dag):
	l = []	
	for n in dag.graph.getNodes():
		if dag.graph.indeg(n) == 0 :
			l+=[n]
	return l
    
    def calculate_flow_unbounded(dag):
        metric = dag.graph["TempCoeff"]
        if dag.nNodes() == 0:	
		return 
	source_list = get_source_nodes(dag)
	for n in source_list:
		for e in dag.graph.getOutEdges(n):
                    m = dag.graph.target(e)
                    metric[m] = metric[m]+metric[n]
                    metric[e] = metric[n]

                g.graph[coeffname][n]+=metric[n]
                dag.graph.delNode(n)	
	return calculate_flow_unbounded(dag)
    
    def transfer_flow(depth):
        dag = g.addSubGraph(InducedDagGraph,node=node, depth = depth, direction = "direct") #get the induced dag
        coeff_old = g.graph[coeffname]
        coeff = dag.graph.getDoubleProperty("TempCoeff")

        coeff.setAllNodeValue(0.0)
        coeff[node] = 1.0 #add the node of the coeff
        coeff_old[node] = 0.0
        
        calculate_flow_unbounded(dag) #calculate the flow on the dag

        g.delSubGraph(dag)

    def transfer_flow_reversed():
        pass
    
    #we assume the node has just been added
    #we assume it's a pubgraph and not a reversepubgraph
    if g.edgeDirection == "direct":
        transfer_flow(depth)
    elif g.edgeDirection == "reversed":
        print "not yet implemented"
        return
    else:
        print "Must be an oriented publication graph"
        return

#TODO add the depth
def MonoplexFlowUnbounded_dynamic(g, depth = float('inf'), coeffname = "FlowUnboundedCoeff"):
    #Do something to check whether of the coeffs have already been calculated or not
    if g.edgeDirection != "direct" and g.edgeDirection != "reversed":
        print "Must be an oriented publication graph"
        return
    
    coeff = g.graph.getDoubleProperty(coeffname)
    gr = g.addSubGraph(EmptyGraph)
    construction_order,edges = GetConstructionOrder(g)

    k=0
    for n in construction_order: #for all remaining nodes
        #add the node and its out edges
        gr.graph.addNode(n)
        for t in edges[n]:
            e = g.graph.existEdge(n,t)
            gr.graph.addEdge(e) #have to do this otherwise it creates duplicate edges... 
        k+=1
        #update the flow of the graph
        MonoplexFlowUnboundedAddNode(gr, n, depth= depth, coeffname = coeffname)
        if k%100 == 0:
            print k
        
    #now the graph's coeffs should be uptodate. Delete the subgraph, we don't need it anymore
    g.delSubGraph(gr)


def KMonoplexFlowUnbounded_static(g, depth=float('inf'), coeffname="KFlowUnboundedCoeff"): 
    if depth==float('inf'):
        MonoplexFlowUnbounded_static(g)
        return

    if depth==0:
        print "bug"
        return
    
    def get_source_nodes(g):
	l = []	
	for n in g.graph.getNodes():
		if g.graph.indeg(n) == 0 :
			l+=[n]
	return l

    flow_ancestry = {m : {i : 0.0 for i in range(depth)} for m in g.graph.getNodes()} #associates to each node the dict [age of flow => quantity of flow]. By convention the age of the flow coming from a neighbor is 0 and it increases by 1 at each transfer.
    
    def calculate_flow(g):
        metric = g.graph[coeffname]
	if g.nNodes() == 0:	
            return 
	source_list = get_source_nodes(g)
        print g.nNodes()

        if g.nNodes()==1:
            for k in range(depth):
                print flow_ancestry[source_list[0]][k]
        
	for n in source_list:
            for e in g.graph.getOutEdges(n):
                m = g.graph.target(e)
                flowint = 1.0
                flow_ancestry[m][0] += flowint #the new flow generated by the node n
                metric[m] += flowint
                metric[e] += flowint
                
                for k in range(depth-1): #now transfer the flow received in n, updating its age
                    kflow = flow_ancestry[n][k]
                    flow_ancestry[m][k+1] = kflow
                    metric[m] += kflow
                    metric[e] += kflow
                    
                        
            g.graph.delNode(n)
            del flow_ancestry[n]
	return calculate_flow(g)

    if g.edgeDirection != "direct" and g.edgeDirection != "reversed":
        print "Must be an oriented publication graph"
        return
    
    if tlp.AcyclicTest.isAcyclic(g.graph) == False:
        print "Graph not acyclic!"
        return
    
    gr = g.addSubGraph(CloneGraph)
            
    metric = g.graph.getDoubleProperty(coeffname)
    metric.setAllNodeValue(1.0)
    metric.setAllEdgeValue(0.0)
    calculate_flow(gr)
    g.delSubGraph(gr)

#########BALANCE VECTOR MEASURE

#Function which calculates the flow balance for each node    
def KBalanceVector_static(g, metricname, coeffname = "KBalanceVectorCoeff"):
    def round_to_zero(x):
        if abs(x) < 0.000000001:
            return 0.0
        else:
            return x
        
    balance = g.graph.getDoubleProperty(coeffname)
    balance.setAllNodeValue(1.0)
    metric = g.graph[metricname]
    for n in g.graph.getNodes():
        for e in g.graph.getInEdges(n):
            balance[n] += metric[e]
        for e in g.graph.getOutEdges(n):
            balance[n] -= metric[e]
        balance[n] = round_to_zero(balance[n])

    
#####MULTIPLEX MEASURES
#measure which attributes to each node the number of nodes on its layer
#THE SAME AS DEGREE
def MultiplexNumNodes_static(g): 
    if g.__class__.__name__ != "MultiplexGraph" and g.__class__.__name__ != "MultiplexDagGraph":
        print "Must be used with a multiplex graph"
        return

    coeff = g.graph.getDoubleProperty("MultiplexNumNodesCoeff")
    for sgk in g.subGraphs:
        node_id = g.subGraphs[sgk].graph.getName().split("mx")[1]
        coeff[root.id2pub[node_id]] = g.subGraphs[sgk].nNodes()
        
#measure which attributes to each node the number of edges on its layer
def MultiplexNumEdges_static(g): 
    if g.__class__.__name__ != "MultiplexGraph" and g.__class__.__name__ != "MultiplexDagGraph":
        print "Must be used with a multiplex graph"
        return

    coeff = g.graph.getDoubleProperty("MultiplexNumEdgesCoeff")
    for sgk in g.subGraphs:
        node_id = g.subGraphs[sgk].graph.getName().split("mx")[1]
        coeff[root.id2pub[node_id]] = g.subGraphs[sgk].nEdges()

#measure which associates the number of connex parts (if we remove the inducing node) in each layer
def MultiplexNumConnex_static(g):
    pass #du bonus    

#measure which associates the longest path in each layer
def MultiplexLongestPath_static(g):
    pass #du bonus
        


#Calculates the flow in each layer and sums them up
def MultiplexFlowSum_static(g, coeffname="MultiplexFlowSumCoeff"):
    def round_to_one(g, node): #rounds to 1 the pathological cases due to floating point arithmetic errors
        if abs(g.graph[coeffname][node] - 1.0) < 0.0000000001:
            g.graph[coeffname][node] = 1.0
    
    def get_source_nodes(g):
	l = []	
	for n in g.graph.getNodes():
		if g.graph.indeg(n) == 0 :
			l+=[n]
	return l
    def calculate_flow(g, coeffname):
        metric = g.graph[coeffname]
	if g.nNodes() == 0:	
		return 
	source_list = get_source_nodes(g)
	for n in source_list:
		for m in g.graph.getOutNodes(n):
			metric[m] = metric[m]+metric[n]/g.graph.outdeg(n)
		for e in g.graph.getOutEdges(n):
			metric[e] = metric[n]/g.graph.outdeg(n)
		g.graph.delNode(n)	
	return calculate_flow(g, coeffname)


    """if g.__class__.__name__ != "MultiplexGraph" and g.__class__.__name__ != "MultiplexDagGraph":
        print "Must be used with a multiplex graph"
        return"""
            

    
    coeff = g.graph.getDoubleProperty(coeffname)
    coeff.setAllNodeValue(0.0)
    coeff.setAllEdgeValue(0.0)
    tempcoeff = g.graph.getDoubleProperty("tempcoeff")

    for gk in g.subGraphs: #for each subgraph calculate the pondered flow
        print gk
        gr = g[gk].addSubGraph(CloneGraph)
        
        #set up the initial coeff values

        sharedCitedLength = gr.graph["sharedCitedLength"]
        for n in gr.graph.getNodes():
            tempcoeff[n] = 1.0/sharedCitedLength[n] #the contribution of the pub in the layer is 1 divided by the number of layers it belongs to


        calculate_flow(gr, "tempcoeff") #calculate the flow in the layer
        
        
        #for each node add its coeff to the multiplex flow sum
        for n in g[gk].graph.getNodes():
            coeff[n] += tempcoeff[n]

        #for each edge add its coeff to the multiplex flow sum
        for e in g[gk].graph.getEdges():
            coeff[e] += tempcoeff[e]

        g[gk].delSubGraph(gr)

    for n in g.graph.getNodes():
        round_to_one(g,n)    

    for e in g.graph.getEdges():
        round_to_one(g,e)   


#Calculates the flow of the whole graph as if each shared edge is an independant edge. I.e. if an edge belongs to N subgraphs, it is considered to be N edges
def MultiplexFlowAggregated_static(g, coeffname = "MultiplexFlowAggregatedCoeff"):
    #Very similar to the flow function, except you ponder the degree of a node by sharedCitedLength. To simulate the fact that there are "multiple edges".
    #Basically, more flow will be transferred between two similar publications
    def get_source_nodes(g):
	l = []	
	for n in g.graph.getNodes():
		if g.graph.indeg(n) == 0 :
			l+=[n]
	return l
    
    def calculate_flow(g):
	if g.nNodes() == 0:	
            return
        metric = g.graph[coeffname]
        sharedCitedLength = g.graph["sharedCitedLength"]
	source_list = get_source_nodes(g)

        for n in source_list:
            totaloutdeg = 0
            for e in g.graph.getOutEdges(n):
                totaloutdeg+=sharedCitedLength[e] #get the total degree, which is the sum of the shared cited lengths of each outgoing edge
            for e in g.graph.getOutEdges(n):
                m = g.graph.target(e) #get the target of the edge
                metric[m] = metric[m]+metric[n]*sharedCitedLength[e]/totaloutdeg 
                metric[e] = metric[n]*sharedCitedLength[e]/totaloutdeg

            g.graph.delNode(n)	

        return calculate_flow(g)
    
    """if g.__class__.__name__ != "MultiplexGraph" and g.__class__.__name__ != "MultiplexDagGraph":
        print "Must be used with a multiplex graph"
        return"""

    gr = g.addSubGraph(CloneGraph)
            
    coeff = g.graph.getDoubleProperty(coeffname)
    coeff.setAllNodeValue(1.0)
    calculate_flow(gr)
    g.delSubGraph(gr)



def KMultiplexFlowAggregated_static(g, depth = float('inf'),coeffname = "KMultiplexFlowAggregatedCoeff"):
    if depth==float('inf'):
        MultiplexFlowAggregated_static(g)
        return

    if depth==0:
        print "bug"
        return
    

    def get_source_nodes(g):
	l = []	
	for n in g.graph.getNodes():
		if g.graph.indeg(n) == 0 :
			l+=[n]
	return l

    flow_ancestry = {m : {i : 0.0 for i in range(depth)} for m in g.graph.getNodes()}
    
    def calculate_flow(g):
	if g.nNodes() == 0:	
            return
        metric = g.graph[coeffname]
        sharedCitedLength = g.graph["sharedCitedLength"]
	source_list = get_source_nodes(g)

        for n in source_list:
            totaloutdeg = 0
            for e in g.graph.getOutEdges(n):
                totaloutdeg+=sharedCitedLength[e] #get the total degree, which is the sum of the shared cited lengths of each outgoing edge
            for e in g.graph.getOutEdges(n):
                m = g.graph.target(e) #get the target of the edge

                flowint = 1.0*sharedCitedLength[e]/totaloutdeg
                flow_ancestry[m][0] += flowint #the new flow generated by the node n
                metric[m] += flowint
                metric[e] += flowint
                
                for k in range(depth-1): #now transfer the flow received in n, updating its age
                    kflow = flow_ancestry[n][k]*sharedCitedLength[e]/totaloutdeg
                    flow_ancestry[m][k+1] = kflow
                    metric[m] += kflow
                    metric[e] += kflow


                
            g.graph.delNode(n)	
            del flow_ancestry[n]
        return calculate_flow(g)
    
    """if g.__class__.__name__ != "MultiplexGraph" and g.__class__.__name__ != "MultiplexDagGraph":
        print "Must be used with a multiplex graph"
        return"""

    gr = g.addSubGraph(CloneGraph)
            
    coeff = g.graph.getDoubleProperty(coeffname)
    coeff.setAllNodeValue(1.0)
    coeff.setAllEdgeValue(0.0)
    calculate_flow(gr)
    g.delSubGraph(gr)



def MultiplexFlowSelective_static(g, coeffname = "MultiplexFlowSelectiveCoeff"):
    #First things first, segregate between the different types of pipes

    def count_elements(H):
        count = 0
        for h in H:
            for k in H[h]:
                count+=1
        print count
    
    def get_source_nodes(g):
	l = []	
	for n in g.graph.getNodes():
		if g.graph.indeg(n) == 0 :
			l+=[n]
	return l

    flow_ancestry = {m : {} for m in g.graph.getNodes()} #initialize the flow ancestry
    
    def calculate_flow(g):
        #flow_ancestry will always have too much information but nobody cares, as long as it has the information on the source nodes
        
        metric = g.graph[coeffname]
        sharedCitedLength = g.graph["sharedCitedLength"]
        sharedCitedBy = g.graph["sharedCitedBy"]
        
        
	if g.nNodes() == 0:	
            return
        print g.nNodes()
	source_list = get_source_nodes(g)

        #transfer the intrinsic value of the publication, i.e. 1, but also transfer the received flow according to its ancestry
        for n in source_list:                                        
            for e in g.graph.getOutEdges(n):
                m  = g.graph.target(e)
                flowint = 1.0/(g.graph.outdeg(n)*sharedCitedLength[e]) #le flot intrinseque a la publication diffuse sur ce chaque sous arete de l'arete
                    
                for o in sharedCitedBy[e]: #on y ajoute le flot herite dependant de la sous arete et on met a jour flow_ancestry
                    p = tlp.node(o)
                    if p not in flow_ancestry[m]:
                        flow_ancestry[m][p] = 0.0
                    if p not in flow_ancestry[n]:
                        flow_ancestry[n][p] = 0.0

                    #Calculate the "outdegree" of the node for the particular layer concerned
                    #can probably make this more efficient by checking the deg of the node in the associated multiplex subgraph
                    pdeg = 0
                    for e in g.graph.getOutEdges(n):
                        if o in sharedCitedBy[e]:
                            pdeg += 1
                        
                    metric[m] += flowint + flow_ancestry[n][p]/pdeg
                    metric[e] += flowint + flow_ancestry[n][p]/pdeg
                    flow_ancestry[m][p] += flowint + flow_ancestry[n][p]/pdeg #update the flow ancestry of the new node

                        

            g.graph.delNode(n)
		                
	return calculate_flow(g)
    
    """if g.__class__.__name__ != "MultiplexGraph" and g.__class__.__name__ != "MultiplexDagGraph":
        print "Must be used with a multiplex graph"
        return"""

    gr = g.addSubGraph(CloneGraph)
            
    coeff = g.graph.getDoubleProperty(coeffname)
    coeff.setAllNodeValue(1.0)
    coeff.setAllEdgeValue(0.0)
    calculate_flow(gr)
    g.delSubGraph(gr)
    print "the element count is"
    count_elements(flow_ancestry)


def MultiplexFlowSelectiveNew_static(g, coeffname = "MultiplexFlowSelectiveCoeff"):
    #get the association between node id and subgraph id
    id2subgraph = {}
    for gkey in g.subGraphs:
        i = g.subGraphs[gkey].name.split('mx')[1]
        id2subgraph[root.id2pub[i].id] = gkey
    
    def count_elements(H):
        count = 0
        for h in H:
            for k in H[h]:
                count+=1
        print count
    
    def get_source_nodes(g):
	l = []	
	for n in g.graph.getNodes():
		if g.graph.indeg(n) == 0 :
			l+=[n]
	return l

    flow_ancestry = {m : {} for m in g.graph.getNodes()} #initialize the flow ancestry
    
    def calculate_flow(g):
        #flow_ancestry will always have too much information but nobody cares, as long as it has the information on the source nodes
        
        metric = g.graph[coeffname]
        sharedCitedLength = g.graph["sharedCitedLength"]
        sharedCitedBy = g.graph["sharedCitedBy"]
        
        
	if g.nNodes() == 0:	
            return
        print g.nNodes()
	source_list = get_source_nodes(g)
        
        
        #transfer the intrinsic value of the publication, i.e. 1, but also transfer the received flow according to its ancestry
        for n in source_list:                                        
            for e in g.graph.getOutEdges(n):
                m  = g.graph.target(e)
                flowint = 1.0/(g.graph.outdeg(n)*sharedCitedLength[e]) #le flot intrinseque a la publication diffuse sur ce chaque sous arete de l'arete
                    
                for o in sharedCitedBy[e]: #on y ajoute le flot herite dependant de la sous arete et on met a jour flow_ancestry
                    p = tlp.node(o)
                    if p not in flow_ancestry[m]:
                        flow_ancestry[m][p] = 0.0
                    if p not in flow_ancestry[n]:
                        flow_ancestry[n][p] = 0.0

                    #Calculate the "outdegree" of the node for the particular layer concerned
                    #can probably make this more efficient by checking the deg of the node in the associated multiplex subgraph
                    pdeg = g.parent.subGraphs[id2subgraph[o]].graph.outdeg(n)
                            
                    metric[m] += flowint + flow_ancestry[n][p]/pdeg
                    metric[e] += flowint + flow_ancestry[n][p]/pdeg
                    flow_ancestry[m][p] += flowint + flow_ancestry[n][p]/pdeg #update the flow ancestry of the new node

                        

            g.graph.delNode(n)
		                
	return calculate_flow(g)
    
    """if g.__class__.__name__ != "MultiplexGraph" and g.__class__.__name__ != "MultiplexDagGraph":
        print "Must be used with a multiplex graph"
        return"""

    gr = g.addSubGraph(CloneGraph)
            
    coeff = g.graph.getDoubleProperty(coeffname)
    coeff.setAllNodeValue(1.0)
    coeff.setAllEdgeValue(0.0)
    calculate_flow(gr)
    g.delSubGraph(gr)
    print "the element count is"
    count_elements(flow_ancestry)


    
def MultiplexFlowSelectiveNewNew_static(g, coeffname = "MultiplexFlowSelectiveCoeff"):
    def round_to_one(g, node): #rounds to 1 the pathological cases due to floating point arithmetic errors
        if abs(g.graph[coeffname][node] - 1.0) < 0.0000000001:
            g.graph[coeffname][node] = 1.0
    
    def get_source_nodes(g):
	l = []	
	for n in g.graph.getNodes():
		if g.graph.indeg(n) == 0 :
			l+=[n]
	return l
    def calculate_flow(g, coeffname):
        metric = g.graph[coeffname]
	if g.nNodes() == 0:	
		return 
	source_list = get_source_nodes(g)
	for n in source_list:
		for m in g.graph.getOutNodes(n):
			metric[m] = metric[m]+metric[n]/g.graph.outdeg(n)
		for e in g.graph.getOutEdges(n):
			metric[e] = metric[n]/g.graph.outdeg(n)
		g.graph.delNode(n)	
	return calculate_flow(g, coeffname)


    """if g.__class__.__name__ != "MultiplexGraph" and g.__class__.__name__ != "MultiplexDagGraph":
        print "Must be used with a multiplex graph"
        return"""
            

    
    coeff = g.graph.getDoubleProperty(coeffname)
    coeff.setAllNodeValue(0.0)
    coeff.setAllEdgeValue(0.0)
    tempcoeff = g.graph.getDoubleProperty("tempcoeff")

    for gk in g.subGraphs: #for each subgraph calculate the pondered flow
        print gk
        gr = g[gk].addSubGraph(CloneGraph)
        
        #set up the initial coeff values

        sharedCitedLength = g.graph["sharedCitedLength"]

        #Put the correct initial conditions 
        for n in gr.graph.getNodes():
            sharen = 0.0
            for e in g.graph.getOutEdges(n):
                sharen+=1.0/(g.graph.outdeg(n)*sharedCitedLength[e]) #get the contribution of the node in the subgraph
            tempcoeff[n] = sharen


        #calculate_flow(gr, "tempcoeff") #calculate the flow in the layer
        
        
        #for each node add its coeff to the multiplex flow sum
        for n in g[gk].graph.getNodes():
            coeff[n] += tempcoeff[n]

        #for each edge add its coeff to the multiplex flow sum
        for e in g[gk].graph.getEdges():
            coeff[e] += tempcoeff[e]

        g[gk].delSubGraph(gr)

    for n in g.graph.getNodes():
        round_to_one(g,n)    

    for e in g.graph.getEdges():
        round_to_one(g,e)   

    
    


####Unbounded versions of all the previous multiplex flow functions
#Don't know what it's worth


def MultiplexFlowSumUnbounded_static(g, coeffname="MultiplexFlowSumUnboundedCoeff"):
    def round_to_one(g, node): #rounds to 1 the pathological cases due to floating point arithmetic errors
        if abs(g.graph[coeffname][node] - 1.0) < 0.0000000001:
            g.graph[coeffname][node] = 1.0
    
    def get_source_nodes(g):
	l = []	
	for n in g.graph.getNodes():
		if g.graph.indeg(n) == 0 :
			l+=[n]
	return l
    def calculate_flow(g, coeffname):
        metric = g.graph[coeffname]
	if g.nNodes() == 0:	
		return 
	source_list = get_source_nodes(g)
	for n in source_list:
		for m in g.graph.getOutNodes(n):
			metric[m] = metric[m]+metric[n]
		for e in g.graph.getOutEdges(n):
			metric[e] = metric[n]
		g.graph.delNode(n)	
	return calculate_flow(g, coeffname)


    """if g.__class__.__name__ != "MultiplexGraph" and g.__class__.__name__ != "MultiplexDagGraph":
        print "Must be used with a multiplex graph"
        return"""
            

    
    coeff = g.graph.getDoubleProperty(coeffname)
    coeff.setAllNodeValue(0.0)
    coeff.setAllEdgeValue(0.0)
    tempcoeff = g.graph.getDoubleProperty("tempcoeff")

    for gk in g.subGraphs: #for each subgraph calculate the pondered flow
        print gk
        gr = g[gk].addSubGraph(CloneGraph)
        
        #set up the initial coeff values

        sharedCitedLength = gr.graph["sharedCitedLength"]
        for n in gr.graph.getNodes():
            tempcoeff[n] = 1.0/sharedCitedLength[n] #la contribution de la publication dans la suos couche est 1 divise par le nombre de sous couches dans laquelle elle apparait

        calculate_flow(gr, "tempcoeff") #on calcule le flot dans le sous graphe
        
        
        #for each node add its coeff to the multiplex flow sum
        for n in g[gk].graph.getNodes():
            coeff[n] += tempcoeff[n]

        #for each edge add its coeff to the multiplex flow sum
        for e in g[gk].graph.getEdges():
            coeff[e] += tempcoeff[e]

        g[gk].delSubGraph(gr)

    for n in g.graph.getNodes():
        round_to_one(g,n)    

    for e in g.graph.getEdges():
        round_to_one(g,e)   


def MultiplexFlowAggregatedUnbounded_static(g, coeffname = "MultiplexFlowAggregatedUnboundedCoeff"):
    def get_source_nodes(g):
	l = []	
	for n in g.graph.getNodes():
		if g.graph.indeg(n) == 0 :
			l+=[n]
	return l
    
    def calculate_flow(g):
	if g.nNodes() == 0:	
            return
        metric = g.graph[coeffname]
        sharedCitedLength = g.graph["sharedCitedLength"]
	source_list = get_source_nodes(g)

        for n in source_list:
            for e in g.graph.getOutEdges(n):
                m = g.graph.target(e) #get the target of the edge
                metric[m] = metric[m]+metric[n]*sharedCitedLength[e] 
                metric[e] = metric[n]*sharedCitedLength[e]

            g.graph.delNode(n)	

        return calculate_flow(g)
    
    """if g.__class__.__name__ != "MultiplexGraph" and g.__class__.__name__ != "MultiplexDagGraph":
        print "Must be used with a multiplex graph"
        return"""

    gr = g.addSubGraph(CloneGraph)
            
    coeff = g.graph.getDoubleProperty(coeffname)
    coeff.setAllNodeValue(1.0)
    calculate_flow(gr)
    g.delSubGraph(gr)



def KMultiplexFlowAggregatedUnbounded_static(g, depth = float('inf'),coeffname = "KMultiplexFlowAggregatedUnboundedCoeff"):
    if depth==float('inf'):
        MultiplexFlowAggregated_static(g)
        return

    if depth==0:
        print "bug"
        return
    

    def get_source_nodes(g):
	l = []	
	for n in g.graph.getNodes():
		if g.graph.indeg(n) == 0 :
			l+=[n]
	return l

    flow_ancestry = {m : {i : 0.0 for i in range(depth)} for m in g.graph.getNodes()} 
    
    def calculate_flow(g):
	if g.nNodes() == 0:	
            return
        metric = g.graph[coeffname]
        sharedCitedLength = g.graph["sharedCitedLength"]
	source_list = get_source_nodes(g)

        for n in source_list:
            for e in g.graph.getOutEdges(n):
                m = g.graph.target(e) #get the target of the edge

                flowint = 1.0*sharedCitedLength[e]
                flow_ancestry[m][0] += flowint #the new flow generated by the node n
                metric[m] += flowint
                metric[e] += flowint
                
                for k in range(depth-1): #now transfer the flow received in n, updating its age
                    kflow = flow_ancestry[n][k]*sharedCitedLength[e]
                    flow_ancestry[m][k+1] = kflow
                    metric[m] += kflow
                    metric[e] += kflow


                
            g.graph.delNode(n)	
            del flow_ancestry[n]
        return calculate_flow(g)
    
    """if g.__class__.__name__ != "MultiplexGraph" and g.__class__.__name__ != "MultiplexDagGraph":
        print "Must be used with a multiplex graph"
        return"""

    gr = g.addSubGraph(CloneGraph)
            
    coeff = g.graph.getDoubleProperty(coeffname)
    coeff.setAllNodeValue(1.0)
    coeff.setAllEdgeValue(0.0)
    calculate_flow(gr)
    g.delSubGraph(gr)




def MultiplexFlowSelectiveUnbounded_static(g, coeffname = "MultiplexFlowSelectiveUnboundedCoeff"):
    def count_elements(H):
        count = 0
        for h in H:
            for k in H[h]:
                count+=1
        print count
    
    def get_source_nodes(g):
	l = []	
	for n in g.graph.getNodes():
		if g.graph.indeg(n) == 0 :
			l+=[n]
	return l

    flow_ancestry = {m : {} for m in g.graph.getNodes()} #initialize the flow ancestry
    
    def calculate_flow(g):
        #flow_ancestry will always have too much information but nobody cares, as long as it has the information on the source nodes
        
        metric = g.graph[coeffname]
        sharedCitedLength = g.graph["sharedCitedLength"]
        sharedCitedBy = g.graph["sharedCitedBy"]
        
        
	if g.nNodes() == 0:	
            return
        print g.nNodes()
	source_list = get_source_nodes(g)

        #transfer the intrinsic value of the publication, i.e. 1, but also transfer the received flow according to its ancestry
        for n in source_list:                                        
            for e in g.graph.getOutEdges(n):
                m  = g.graph.target(e)
                flowint = 1.0/(sharedCitedLength[e]) #le flot intrinseque a la publication diffuse sur ce chaque sous arete de l'arete
                    
                for o in sharedCitedBy[e]: #on y ajoute le flot herite dependant de la sous arete et on met a jour flow_ancestry
                    p = tlp.node(o)
                    if p not in flow_ancestry[m]:
                        flow_ancestry[m][p] = 0.0
                    if p not in flow_ancestry[n]:
                        flow_ancestry[n][p] = 0.0

                    metric[m] += flowint + flow_ancestry[n][p]
                    metric[e] += flowint + flow_ancestry[n][p]
                    flow_ancestry[m][p] += flowint + flow_ancestry[n][p] #update the flow ancestry of the new node

                        

            g.graph.delNode(n)
		                
	return calculate_flow(g)
    
    """if g.__class__.__name__ != "MultiplexGraph" and g.__class__.__name__ != "MultiplexDagGraph":
        print "Must be used with a multiplex graph"
        return"""

    gr = g.addSubGraph(CloneGraph)
            
    coeff = g.graph.getDoubleProperty(coeffname)
    coeff.setAllNodeValue(1.0)
    coeff.setAllEdgeValue(0.0)
    calculate_flow(gr)
    g.delSubGraph(gr)
    print "the element count is"
    count_elements(flow_ancestry)
    


#####AUTHOR MEASURES
#add more measures. Centrality, size of the associated cluster, eccentricity, etc.
#Author parsing still buggy. Known bugs :
#problem with jr. counted as an author
#some authors are just initials
#juan malcadena not associated with his initials

#measure which attributes the number of coauthors
def MostConnectedAuthors_static(g): #should be used with author graph
    if g.__class__.__name != "AuthorGraph":
        print "Must be used with Author Graph"
        return

    coeff = g.graph.getIntegerProperty("MostConnectedAuthorsCoeff")
    for n in g.graph.getNodes():
        coeff[n] = g.graph.outdeg(n) #no difference if we choose indeg or outdeg
    
    
        
#####Function to display the graph taking into account a coeff
def DisplayMonoplexFlow(g, coeffname):
    metric = g.graph[coeffname]
    
    viewSize = g.graph.getSizeProperty("viewSize")	
    viewLabel = g.graph.getStringProperty("viewLabel")
    viewColor = g.graph.getColorProperty("viewColor")
	
    dataset = tlp.getDefaultPluginParameters("Size Mapping")
    dataset['property'] = metric
    dataset['max size'] = 5
    dataset['min size'] = 0.01
	
    g.graph.applySizeAlgorithm('Size Mapping', viewSize, dataset)

    dataset = tlp.getDefaultPluginParameters('Color Mapping')
    dataset['input property'] = metric
    g.graph.applyColorAlgorithm('Color Mapping', graph["viewColor"], dataset)
	
    #for it to be less ugly..
    baseSizeEdge = tlp.Size(0.125,0.125,0.5	)
    for n in graph.getEdges():
        viewSize[n] = baseSizeEdge
        viewLabel[n] = ""
        viewColor[n] = tlp.Color.Byzantium
	
#Function to reset the display of the graph
def ResetDisplayGraph(g):
    viewSize = g.graph.getSizeProperty("viewSize")	
    viewLabel = g.graph.getStringProperty("viewLabel")
    viewColor = g.graph.getColorProperty("viewColor")
    baseSizeNode = tlp.Size(1,1,1)
    baseSizeEdge = tlp.Size(0.125,0.125,0.5	)

    for n in g.graph.getNodes():
        viewSize[n] = baseSizeNode
        viewLabel[n] = ""
        viewColor[n] = tlp.Color.Red
    for n in g.graph.getEdges():
        viewSize[n] = baseSizeEdge
        viewLabel[n] = ""
        viewColor[n] = tlp.Color.Gray
