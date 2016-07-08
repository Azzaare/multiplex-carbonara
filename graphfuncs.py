from tulip import *
from graphclass import *
import numpy as np


##################GRAPH MODIFICATION AND MANIPULATION RELATED FUNCTIONS

def MakeAcyclic(g): #function to make the graph acyclic
    def remove_loops(g):
	for e in g.graph.getEdges(): #go through all the edges
            if g.graph.source(e) == g.graph.target(e): #if it's a loop delete it
                g.graph.delEdge(e)
                
    def remove_time_travel(g):
	dateOrd = g.graph["ordDate"]
	for e in g.graph.getEdges(): #make sure that there are no edges from the future to the past
            if dateOrd[g.graph.source(e)] < dateOrd[g.graph.target(e)]:
                g.graph.delEdge(e)

    def remove_time_travel_reversed(g):
	dateOrd = g.graph["ordDate"]
	for e in g.graph.getEdges(): #make sure that there are no edges from the future to the past
		if dateOrd[g.graph.source(e)] > dateOrd[g.graph.target(e)]:
			g.graph.delEdge(e)

    def replace_parallel_edges(g): #replace parallel edges to keep the citation flow
	for e in g.graph.getEdges(): #go through all the edges
            if g.graph.isElement(e) == False: #if it's already been deleted
                continue
            er = g.graph.existEdge(g.graph.target(e),g.graph.source(e)) #check if the reverse edge exists
            if er.isValid(): #if it does, delete both edges	
                g.graph.delEdge(e)
                g.graph.delEdge(er)

    def replace_parallel_edges_fusion(g): #replace parallel edges to keep the citation flow
        node_date = g.graph["ordDate"]
	node_authors = g.graph["authors"]
        node_type = g.graph["node_type"]
	node_property = g.graph["node_property"]
        #graph properties to update the properties of the newly created fusion nodes

        
        for e in g.graph.getEdges(): #go through all the edges
            if g.graph.isElement(e) == False: #if it's already been deleted
                continue
            er = g.graph.existEdge(g.graph.target(e),g.graph.source(e)) #check if the reverse edge exists
            if er.isValid(): #if it does, fusion of the two nodes and creation of all edges                
                fusion_node = g.graph.addNode()
                n1 = g.graph.source(e)
                n2 = g.graph.target(e)

                #out of curiosityu, what are the authors?
                print g.graph["authors"][n1]
                print g.graph["authors"][n2]
                print ""
                
                if n1 == n2:
                    print "buuuuuug"
                #update properties of the fusion node
                node_type[fusion_node] = "publication"
                node_property[fusion_node] = node_property[n1]+";"+node_property[n2]

                if node_date[n1] != node_date[n2]:
                    print "BUG DATE!!!!"
                node_date[fusion_node] = node_date[n1]
                node_authors[fusion_node] = list(set(node_authors[n1]+node_authors[n2]))

                
                g.graph.delEdge(e)
                g.graph.delEdge(er)	
                for e in g.graph.getInEdges(n1):
                    if g.graph.existEdge(g.graph.source(e),fusion_node).isValid() == False: #make sure we dont create the edge twice
                        g.graph.setTarget(e,fusion_node)
                for e in g.graph.getOutEdges(n1):
                    if g.graph.existEdge(fusion_node,g.graph.target(e)).isValid() == False:
                        g.graph.setSource(e,fusion_node)
                for e in g.graph.getInEdges(n2):
                    if g.graph.existEdge(g.graph.source(e),fusion_node).isValid() == False:
                        g.graph.setTarget(e,fusion_node)
                for e in g.graph.getOutEdges(n2):
                    if g.graph.existEdge(fusion_node,g.graph.target(e)).isValid() == False:
                        g.graph.setSource(e,fusion_node)
                g.graph.delNode(n1)
                g.graph.delNode(n2)

    if tlp.AcyclicTest.isAcyclic(g.graph) == True:
        print "Graph is already acyclic!"
        return

    print "Making Acyclic..."
    remove_loops(g)
    replace_parallel_edges(g)#_fusion(g)
    if g.edgeDirection == "direct":
        remove_time_travel(g)
    elif g.edgeDirection == "reversed":
        remove_time_travel_reversed(g)


#Returns for a given node the node and all his sons
def GetSons(g,node,direction = "direct"):
    
    node_list = [node]
    
    if g.edgeDirection != direction:
        getDirectNodes = g.graph.getInNodes
    else:
        getDirectNodes = g.graph.getOutNodes
        
    for n in getDirectNodes(node):
        node_list+=[n]
        

    node_list = set(node_list)
    return node_list



#returns the set of nodes in the dag induced by a node of given depth and given direction
def SelectDag(g, node, depth = float('inf'), direction = "direct"):
    if g.edgeDirection != direction:
        getDirectNodes = g.graph.getInNodes
    else:
        getDirectNodes = g.graph.getOutNodes 

    if depth==1:
        sons = set(get_out_nodes(node))
        return sons

    node_list = []
    node_status = {n : 0 for n in g.graph.getNodes()} #1 if visited 0 otherwise
    node_depth = {}
    node_depth[node] = 0
    node_status[node] = 1
    to_visit_nodes = [node]

    def get_out_nodes(node):
        l = []
        for n in getDirectNodes(node):
            if node_status[n] == 0 and node_depth[node] < depth:
                l+=[n]
                node_status[n] = 1
                node_depth[n] = node_depth[node]+1
        return l
    
    while len(to_visit_nodes) != 0:
        node_list+=[to_visit_nodes[0]]
        
        l = get_out_nodes(to_visit_nodes[0])
            
        to_visit_nodes = to_visit_nodes[1:] #remove first element
        to_visit_nodes += l #update the list of nodes to visit
                
    node_list = set(node_list) #convert to set
    
    return node_list

#returns a dictionary whose keys are all the nodes in the subDAG induced by a given node, and whose values are the distances from the given node.
def SelectDagDistance(g, node, depth = float('inf'), direction = "direct"):
    if g.edgeDirection != direction:
        getDirectNodes = g.graph.getInNodes
    else:
        getDirectNodes = g.graph.getOutNodes 

    if depth==1:
        sons = set(get_out_nodes(node))
        return sons

    node_list = {}
    node_status = {n : 0 for n in g.graph.getNodes()} #1 if visited 0 otherwise
    node_depth = {}
    node_depth[node] = 0
    node_status[node] = 1
    to_visit_nodes = [node]

    def get_out_nodes(node):
        l = []
        for n in getDirectNodes(node):
            if node_status[n] == 0 and node_depth[node] < depth:
                l+=[n]
                node_status[n] = 1
                node_depth[n] = node_depth[node]+1
        return l
    
    while len(to_visit_nodes) != 0:
        node_list[to_visit_nodes[0]] = node_depth[to_visit_nodes[0]]
        
        l = get_out_nodes(to_visit_nodes[0])
            
        to_visit_nodes = to_visit_nodes[1:] #remove first element
        to_visit_nodes += l #update the list of nodes to visit
                
    
    return node_list
                

##############Dynamic construction related functinos

#returns a tuple, containing an ordered list of nodes to add, and a dictionary with the edges to add for each node.
def GetConstructionOrder(g): #gives the construction order of a graph given in entry. Must be a dag
    def get_roots(g):
        roots = []
        for n in g.graph.getNodes():
            if g.graph.outdeg(n) == 0:
                roots+=[n]
        return roots
    
    gr = g.addSubGraph(CloneGraph)
    node_order = [] #the order of the nodes to add
    edges = {} #the edges to add when you add a particular node (always out edges by definition)
                
    while gr.nNodes() != 0:
        roots = get_roots(gr)
        for n in roots:
            node_order+=[n]
            le = []
            for e in g.graph.getOutEdges(n):
                le += [g.graph.target(e)]
            edges[n] = le

        for n in roots:
            gr.graph.delNode(n)

    g.delSubGraph(gr)
    return (node_order,edges)
    
#test out the node construction order given by the above function
def TestConstructionOrder(construction_order): 
    node_order = construction_order[0]
    edges = construction_order[1]
    g = root.addSubGraph(EmptyGraph)
    for n in node_order:
        g.graph.addNode(n)
        for t in edges[n]:
            g.graph.addEdge(n,t)

    return g


#Reverse the edges in a graph
def ReverseEdges(g): 
    print "Reversing edges of" #for debugging purposes
    print g
    k=0
    for e in g.graph.getEdges():
        g.graph.reverse(e)
        k+=1
        #print k
    if g.edgeDirection == "reversed":
        g.edgeDirection = "direct"
    elif g.edgeDirection == "direct":
        g.edgeDirection = "reversed"


        
#############Clique functions
#make it into a class
    
def GetCliques(g, depth): #make it better, more efficient
	k=0
	gr = g.graph.addSubGraph(EmptyGraph,"shared cliques")
	
	
	cliques = {} #take edge as key
	
	for n in g.graph.getNodes():
            for m1 in g.graph.getInNodes(n):
                for m2 in g.graph.getInNodes(n):
                    if m1 == m2:
                        continue
                    cliques[(m1,m2)] = 1
		k=k+1
		if k%depth == 0:
                    break
	print "cliques loaded"
	print "len",len(cliques)	
	for c in cliques:
            gr.graph.addNode(c[0])
            gr.graph.addNode(c[1])
            gr.graph.addEdge(c[0],c[1])
        print "clique graph constructed"
        return gr

#do something with this

################Connectivity Maxflow Relabel Function
##the relabel to front maxflow algorithm
def relabel_to_front(graph, source, sink): 
    #create the preflow, excess and height functions
    flow = {}
    height = {}
    excess = {}
	
    capacities={}
    for e in graph.getEdges():	#set all capacities to 1
        capacities[e] = 1
	
    def print_info():
        print "flow"
        print flow	
        #print "height"
        #print height
        #print "excess"
        #print excess
		
    #initialize the preflow,excess and height functions
    for e in graph.getEdges():
        flow[e] = 0
    for n in graph.getNodes():
        height[n] = 0
        excess[n] = 0
    height[source] = int(graph.numberOfNodes())
    for e in graph.getOutEdges(source):
        flow[e] = capacities[e] 
        excess[graph.target(e)] = capacities[e]
        excess[source] -= capacities[e]
	
		
	
    #define the push and relabel operations
    def push(u,v): 
        e = graph.existEdge(u,v)
        if e.isValid() == True:	#if it's an edge in the real graph
            cf = capacities[e]-flow[e]
            df = min(excess[u],cf)
            flow[e] += df
        else: #if it's in the residual graph
            e = graph.existEdge(v,u)
            cf = flow[e]
            df = min(excess[u],cf)
            flow[e]-=df
        excess[u]-=df
        excess[v]+=df
		
    def relabel(u):
        #get the dests of the edges in the residual graph, where u is the source
        l=[]	
        for e in graph.getOutEdges(u):
            if capacities[e]-flow[e] > 0:
                l+=[height[graph.target(e)]]
        #get the source of the edges in the residual graph, where u is the dest
        for e in graph.getInEdges(u):
            if flow[e] > 0:
                l+=[height[graph.source(e)]]
			
			
        height[u] = 1 + min(l)
	
    neighbors = {} #neighbor hashtable
    neighborsIndex = {}	
    for u in graph.getNodes():
        neighbors[u] = []	
        for n in graph.getInNodes(u):
            neighbors[u]+=[n]
        for n in graph.getOutNodes(u):
            neighbors[u]+=[n]		
        neighborsIndex[u] = 0
		
    def discharge(u):
        while excess[u] > 0 :
            if neighborsIndex[u] == len(neighbors[u]):
                relabel(u)
                neighborsIndex[u] = 0
            else:
                v = neighbors[u][neighborsIndex[u]]
                e = graph.existEdge(u,v)
                if e.isValid()==True:
                    cf = capacities[e]-flow[e]
                else:
                    e = graph.existEdge(v,u)
                    cf = flow[e] 
                if cf > 0 and height[u] == height[v]+1:
                    push(u,v)
                else:
                    neighborsIndex[u]+=1
	
    L = []
    for i in graph.getNodes():
        L += [i]
    L.remove(source)
    L.remove(sink)
	
    Lindex = 0	
    while Lindex < len(L):
        u = L[Lindex]	
        oldheight = height[u]
        discharge(u)
        if height[u] > oldheight:
            L.remove(u)
            L = [u] + L
            Lindex = 0
        Lindex+=1	
    print_info()
	
    return flow

#determines the connectivity between two nodes. Returns the list of edge connecting the two graphs as well as the value of the maxflow between the two nodes
def Connectivity(g,source,sink):
    flow = relabel_to_front(g.graph,source,sink)
    edge_list = []
    for e in flow:
        if flow[e] > 0:
            edge_list+=[e]
            		
    #And now determine the maxflow
    i=0	
    for e in g.graph.getInEdges(sink):
        if flow[e] > 0:
            i+=1
    return (edge_list,i)

############Select infection clusters
def InfectClusters(g, node, direction = "direct"): #we only do depth==2
    reverseEdges=False
    if direction != g.edgeDirection:
        reverseEdges = True

    lclusters = []
    if reverseEdges == False: #if we have a direct graph
        for n in g.graph.getInNodes(node):
            for m in g.graph.getInNodes(n):
                lnodes = [node,n,m]
                lclusters+=[lnodes]
    else:
        for n in g.graph.getOutNodes(node):
            for m in g.graph.getOutNodes(n):
                lnodes = [node,n,m]
                lclusters+=[lnodes]

    return lclusters


##########Add author nodes to a pubgraph
def AddAuthorNodes(g):
    #Get author nodes
    autlist = []
    authors = g.graph["authors"]
    for n in g.graph.getNodes():
        for aut in authors[n]:
            autlist+= [root.id2aut[aut]]
    #remove duplicates
    autlist = list(set(autlist))
    
    #Add author nodes
    for aut in autlist:
        g.graph.addNode(aut)

    #Add edges
    for n in g.graph.getNodes():
        for aut in authors[n]:
            g.graph.addEdge(n,root.id2aut[aut])
    

##########EXPERIMENTS WITH EIGENVECTOR CENTRALITY CALCULATIONS
#Worthless

#Function to convert the graph into a numpy matrix
def ConvertToMatrix(g):#, direction = "direct"):
    H = {}
    k = 0
    for n in g.graph.getNodes():
        H[n] = k
        k+=1
        
    a = np.zeros((g.nNodes(),g.nNodes()))
    for n in g.graph.getNodes():
        for m in g.graph.getInNodes(n):
              a[H[n]][H[m]] = 1

        for m in g.graph.getOutNodes(n):
            a[H[n]][H[m]] = -1

    return a
            
#Function which determines the eigenvector coeff for each node in a graph. Returns the eigenvalue
def EigenValueMeasure(g, coeff = "EigenCoeff"):
    M = ConvertToMatrix(g)

    (a,b) = np.linalg.eig(M)
    print a
    print b

def PowerIteration(g, niter = 10002):
    x = np.array([0.,0.,1.]) #starting vector
    M = ConvertToMatrix(g)
    print M
    for n in range(niter):
        tx = x
        x = M.dot(x)
        x += 0.5*tx
        x = x/np.linalg.norm(x)
    print x
    
