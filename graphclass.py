from tulip import *
from graphfuncs import *
import networkx as nx
import citation_graph as cg
import datetime
import abc
import random



#############################CLASS DECLARATIONS


#initialize the global variable base_graph

#The abstract class for the graphclass.
#The goal of this is the abstract method build_graph. This method is called in the constructor, and we can hence create classes with different build_graph methods. Then when we create objects of this class it will automatically build a graph of a certain type.
class AbstractGraph(object):
    __metaclass__ = abc.ABCMeta
    
   
    @abc.abstractproperty
    def _default_name(self): #abstract property for the default name of the type of graph
        pass
    
    def __init__(self, parent = None, name = "", **args): #constructor
        self.parent = parent
        if name != "":
            self.name = name
        else:
            self.name = self._default_name()
            
        self.subGraphs = {} #initialize the list of subgraphs
        self.edgeDirection = "undefined"

        self.build_graph(**args)
        """if len(args) > 0:
            self.build_graph(args)
        else:
            self.build_graph()"""
        
    @abc.abstractmethod
    def build_graph(self): #abstract method to build the graph
        pass

    def getId(self): #get graph ID
        return self.graph.getId()

    def nNodes(self): #get number of nodes
        return self.graph.numberOfNodes()

    def nEdges(self): #get number of edges
        return self.graph.numberOfEdges()
    
    def numberOfSubGraphs(self): #get number of subgraphs
        return len(self.subGraphs)
    
    def printSubGraphs(self): #print subgraphs
        for gr in self.subGraphs:
            print self.subGraphs[gr]

    def printAllSubGraphs(self, tabd = 0): #print also the subgraphs of the subgraphs
        for gr in self.subGraphs:
            print "\t"*tabd,self.subGraphs[gr]
            if self.subGraphs[gr].numberOfSubGraphs() > 0:
                self.subGraphs[gr].printAllSubGraphs(tabd+1)

    def __getitem__(self, id): #get subgraph by ID
        return self.subGraphs[id]

    def __str__(self):
        return str(self.getId())+" : " + self.name

    def addSubGraph(self, GraphClass, name = "", **args):
        g = GraphClass(self,name,**args)
        return g

    def delSubGraph(self,g):
        del self.subGraphs[g.getId()]
        if g.__class__.__name__ !="NetworkxGraph": #networkxgraph being a special case where there is no tulip subgraph
            self.graph.delAllSubGraphs(g.graph)
        
        
    def delSubGraphs(self):
        for key in self.subGraphs:
            if self.subGraphs[key].__class__.__name__ !="NetworkxGraph": #networkxgraph being a special case where there is no tulip subgraph
                self.graph.delAllSubGraphs(self.subGraphs[key].graph)
        self.subGraphs = {}

    def renameGraph(self, name):
        self.name = name
        self.graph.setName(name)


    #Function which replaces current graph with one of its subgraphs in the graph hierarchy. This will not actually delete the graph, so please use sparingly
    def ReplaceWithSubGraph(self, sub): 
        self.parent.subGraphs[sub.getId()] = sub
        sub.parent = self.parent
        del self.parent.subGraphs[self.getId()]

    
#Base class for the graphs.    
class BaseGraph(AbstractGraph): 
    #static data, useful for building the subgraphs
    pub_data = cg.load_data()
    autlist = cg.author_list(pub_data)
    id2days = {}
    id2months = {}
    id2years = {}
    id2pub = {}
    id2aut = {}
    
    def _default_name(self):
        return "Base Graph"

    def build_graph(self):
        self.graph = tlp.newGraph()
        self.graph.setName(self.name)
        #get stringproperty
	node_type = self.graph.getStringProperty("node_type")
        node_property = self.graph.getStringProperty("node_property")	
	
	#Create date nodes
	for i in range(1,32):
            n = self.graph.addNode()
            node_type[n] = "day"
            node_property[n] = str(i)
            BaseGraph.id2days[i] = n

	for i in range(1,13):
            n = self.graph.addNode()
            node_type[n] = "month"
            node_property[n] = str(i)
            BaseGraph.id2months[i] = n	
	
	for i in range(1992,2004):
            n = self.graph.addNode()
            node_type[n] = "year"
            node_property[n] = str(i)
            BaseGraph.id2years[i] = n	

        #create author nodes
	for aut in BaseGraph.autlist:
            n = self.graph.addNode()
            node_type[n] = "author"
            node_property[n] = aut
            BaseGraph.id2aut[aut] = n
	
	#Create publication nodes
	for id in BaseGraph.pub_data:	
            n = self.graph.addNode()
            node_type[n] = "publication"
            node_property[n] = str(id)	
            BaseGraph.id2pub[str(id)] = n

        #Add authors and dates		
	node_date = self.graph.getIntegerProperty("ordDate")
	node_authors = self.graph.getStringVectorProperty("authors")	
	
	for i in BaseGraph.pub_data:
            node_authors[BaseGraph.id2pub[str(i)]] = list(set(BaseGraph.pub_data[i][2]))
            l = BaseGraph.pub_data[i][1]
            d = datetime.date(l[2],l[1],l[0])
            node_date[BaseGraph.id2pub[str(i)]] = d.toordinal()

    #methods to save and load the graph

    def saveGraph(self, filename = "graph.tlp"):
        tlp.saveGraph(self.graph, filename)

    def loadGraph(self, filename = "graph.tlp"): #load the graph and all its subgraphs
        def loadSubGraphs(g):
            for sg in g.graph.getSubGraphs():
                #create the subgraphs, initialize them and call this function on them
                sk = g.addSubGraph(SkeletonGraph, sg.getName())
                sk.append_graph(sg)
            #now call the function for all the subgraphs
            for sg in g.subGraphs:
                loadSubGraphs(g[sg])

        #initialize root
        self.delSubGraphs() #delete all previous subgraphs
        self.graph = tlp.loadGraph(filename)
        #add its subgraphs
        loadSubGraphs(self)
        
        
#global variable root which is the base graph. Will be the root graph of all subsequent graphs
root = BaseGraph() 
            
#class to create an empty graph
class EmptyGraph(AbstractGraph): 
    def _default_name(self):
        return "Empty Graph"
    def build_graph(self):
        if self.parent == None:
            raise Exception("No parent graph specified")
        
        self.graph = self.parent.graph.addSubGraph(self.name) #create the graph
        self.parent.subGraphs[self.graph.getId()] = self #add its id to the parent graph
        self.edgeDirection = self.parent.edgeDirection #take the same edge directino as its parent

#No graph actually created, just a skeleton object        
class SkeletonGraph(AbstractGraph): 
    def _default_name(self):
        return "Skeleton Graph"
    def build_graph(self):
        if self.parent == None:
            raise Exception("No parent graph specified")
    
    def append_graph(self, appgraph): #append the specified graph
        self.graph = appgraph
        self.parent.subGraphs[self.graph.getId()] = self #add its id to the parent graph
        self.edgeDirection = self.parent.edgeDirection #take the same edge direction as its parent

#class to create a clone graph
class CloneGraph(AbstractGraph): 
    def _default_name(self):
        return "Clone Graph of "+self.parent.name
    def build_graph(self):
        if self.parent == None:
            raise Exception("No parent graph specified")
        
        self.graph = self.parent.graph.addCloneSubGraph(self.name) #create the graph
        self.parent.subGraphs[self.graph.getId()] = self #add its id to the parent graph
        self.edgeDirection = self.parent.edgeDirection #take the same edge directino as its parent


#class to create a publications graph
class PubGraph(AbstractGraph):
    def _default_name(self):
        return "Publication Graph"

    def build_graph(self):
        if self.parent == None:
            raise Exception("No parent graph specified")
        
        self.graph = self.parent.graph.addSubGraph(self.name) #create the graph
        self.parent.subGraphs[self.graph.getId()] = self #add its id to the parent graph
        self.edgeDirection = "direct"

        #Add publication nodes
	for i in root.id2pub:
		self.graph.addNode(root.id2pub[i])
	
	
	#Create graph edges
	for i in root.pub_data:
		for j in root.pub_data[i][0]:
			self.graph.addEdge(root.id2pub[str(i)], root.id2pub[str(j)])

        #Make Acyclic
        MakeAcyclic(self)


#class to create a publications graph with reversed edges
class ReversePubGraph(AbstractGraph): 
    def _default_name(self):
        return "Reversed Edges Publication Graph"

    def build_graph(self):
        if self.parent == None:
            raise Exception("No parent graph specified")
        
        self.graph = self.parent.graph.addSubGraph(self.name) #create the graph
        self.parent.subGraphs[self.graph.getId()] = self #add its id to the parent graph
        self.edgeDirection = "reversed" 
        #Add publication nodes
	for i in root.id2pub:
		self.graph.addNode(root.id2pub[i])
	
	
	#Create graph edges
	for i in root.pub_data:
		for j in root.pub_data[i][0]:
			self.graph.addEdge(root.id2pub[str(j)], root.id2pub[str(i)])

        #Make Acyclic
        MakeAcyclic(self)



#class to create an author graph. Uses the author data of the parent graph
class AuthorGraph(AbstractGraph): 
    def _default_name(self):
        return "Author Graph"

    def build_graph(self, metricname = None, TransferFunc = None, coeffname = None):
        if self.parent == None:
            raise Exception("No parent graph specified")

        self.graph = self.parent.graph.addSubGraph(self.name)
        self.parent.subGraphs[self.graph.getId()] = self #add its id to the parent graph
        
	shared_pubs = self.graph.getStringVectorProperty("shared_pubs")
        
        #Get author nodes
        autlist = GetAuthors(self.parent)
        authors = self.parent.graph["authors"]
                    
        #Add author nodes
	for aut in autlist:
            self.graph.addNode(aut)
		
	#Create graph edges
        for n in self.parent.graph.getNodes():
            id = self.parent.graph["node_property"][n]
            for aut1 in authors[n]:
                for aut2 in authors[n]:
                    if aut1 != aut2:
                        e = self.graph.existEdge(root.id2aut[aut1], root.id2aut[aut2])
                        if e.isValid(): #if there already exists an edge binding the two authors
                            er = self.graph.existEdge(root.id2aut[aut2], root.id2aut[aut1]) #update reverse edge too
                            shared_pubs[e] += [id]
                            shared_pubs[er] += [id]
                        else: #if the authors are not yet connected with an edge
                            e = self.graph.addEdge(root.id2aut[aut1],root.id2aut[aut2])
                            er = self.graph.addEdge(root.id2aut[aut2],root.id2aut[aut1])
                            shared_pubs[e] += [id]
                            shared_pubs[er] += [id]

        #now transfer the values, if need be
        if metricname!=None and TransferFunc!=None and coeffname!=None:
            pass #TODO
            
      

##Attention lorsqu'on utilise induceddag avec direction= "reversed" ne selectionne pas les auteurs
#class to create graph with both publications and authors
class PubAuthorGraph(AbstractGraph): 
    def _default_name(self):
        return "Publication and Author Graph"

    def build_graph(self):
        if self.parent == None:
            raise Exception("No parent graph specified")
        
        self.graph = self.parent.graph.addSubGraph(self.name)
        self.parent.subGraphs[self.graph.getId()] = self #add its id to the parent graph
	self.edgeDirection = "direct"
	#Add publication nodes
	for i in root.id2pub:
            self.graph.addNode(root.id2pub[str(i)])
				
	#Add author nodes
	for aut in root.id2aut:
            self.graph.addNode(root.id2aut[aut])
		
	#add pub and aut edges
	for i in root.pub_data:
            for j in root.pub_data[i][0]:
                self.graph.addEdge(root.id2pub[str(i)], root.id2pub[str(j)])
            for aut in root.pub_data[i][2]:
                self.graph.addEdge(root.id2pub[str(i)],root.id2aut[aut])

        #Make Acyclic
        MakeAcyclic(self)


#Class which demultiplexes the parent graph into subgraphs containing the neighbors of a node.
class MultiplexGraph(AbstractGraph): 
    def _default_name(self):
        return "Multiplex Graph"

    def build_graph(self):
        
        if self.parent == None:
            raise Exception("No parent graph specified")

        self.graph = self.parent.graph.addCloneSubGraph(self.name) #create the graph
        self.parent.subGraphs[self.graph.getId()] = self #add its id to the parent graph
        self.edgeDirection = self.parent.edgeDirection

        if (self.edgeDirection != "direct"):
            print "must be direct"
            return
        k = 0

        shared_cited_by = self.graph.getIntegerVectorProperty("sharedCitedBy") #on each edge indicates the papers which cites both nodes of the edge
        shared_cited_length = self.graph.getIntegerProperty("sharedCitedLength")

        node_id = self.graph["node_property"]
        
        for n in self.graph.getNodes():
            gr = self.addSubGraph(SkeletonGraph, name="mx"+node_id[n])
            sons = GetSons(self,n,direction = "reversed")
             
            tgraph = self.graph.inducedSubGraph(sons)
            tgraph.setName("mx"+node_id[n])
            gr.append_graph(tgraph)
            
            for e in gr.graph.getEdges():
                shared_cited_by[e]+=[n.id]
                shared_cited_length[e]+=1

            for m in gr.graph.getNodes():
                shared_cited_length[m]+=1    
            k+=1
            print k

        print "graphs built"
    
#Class which demultiplexes the parent graph into the subdags associated with a node. Default depth of dags is infinity
class MultiplexDagGraph(AbstractGraph): 
    def _default_name(self):
        return "Multiplex Dag Graph"
    
    def build_graph(self, depth = float('inf')):
        
        if self.parent == None:
            raise Exception("No parent graph specified")

        self.graph = self.parent.graph.addCloneSubGraph(self.name) #create the graph
        self.parent.subGraphs[self.graph.getId()] = self #add its id to the parent graph
        self.edgeDirection = self.parent.edgeDirection

        k = 0

        shared_cited_by = self.graph.getIntegerVectorProperty("sharedCitedBy") #on each edge indicates the papers which cites both nodes of the edge
        shared_cited_length = self.graph.getIntegerProperty("sharedCitedLength")
        node_id = self.graph["node_property"]
        
        if depth==1: #if depth==1 we can use getsons instead of selectdag, which is quicker
            for n in self.graph.getNodes():
                gr = self.addSubGraph(SkeletonGraph, name="mx"+node_id[n])
                sons = GetSons(self,n,direction = "reversed")
             
                tgraph = self.graph.inducedSubGraph(sons)
                tgraph.setName("mx"+node_id[n])
                gr.append_graph(tgraph)
                for e in gr.graph.getEdges():
                    shared_cited_by[e]+=[n.id]
                    shared_cited_length[e]+=1

                for m in gr.graph.getNodes():
                    shared_cited_length[m]+=1    
                k+=1
                print k

        else:
            for n in self.graph.getNodes():
                gr = self.addSubGraph(InducedDagGraph, name = "mx"+node_id[n], node = n, depth = depth, direction = "reversed")
                for e in gr.graph.getEdges():
                    shared_cited_by[e]+=[n.id]
                    shared_cited_length[e]+=1

                for m in gr.graph.getNodes():
                    shared_cited_length[m]+=1
                    
                k+=1
                print k


        print "graphs built"

        
#class to create the induced DAG of a node  
class InducedDagGraph(AbstractGraph):
    def _default_name(self):
        return "Induced Dag Graph"

    def build_graph(self, node, depth = float('inf'), direction = "direct"):
        
        if self.parent == None:
            raise Exception("No parent graph specified")
        
        node_list = SelectDag(self.parent, node, depth, direction)
        #now build the subgraph 
        self.graph = self.parent.graph.inducedSubGraph(node_list) #create the graph
        if self.name == self._default_name(): #update the name if not custom name
            self.name += " from node "+str(node.id)+" of depth "+str(depth)
        
        self.graph.setName(self.name) #set name of the graphx
        self.parent.subGraphs[self.graph.getId()] = self #add its id to the parent graph
        self.edgeDirection = self.parent.edgeDirection

#class to create a random graph
class InducedRandomGraph(AbstractGraph):
    def _default_name(self):
        return "Induced Random Graph"

    def build_graph(self, numnodes):
        if self.parent == None:
            raise Exception("No parent graph specified")

        if numnodes > self.parent.nNodes():
            print "numnodes higher than total number of nodes"
            return
        
        visited_nodes = { n : 0 for n in self.parent.graph.getNodes()}
        node_list = []
        while len(node_list) < numnodes:
            n = self.parent.graph.getRandomNode()
            if visited_nodes[n] == 0:
                visited_nodes[n] = 1
                node_list+=[n]
        node_list = set(node_list)
        
        #now build the subgraph 
        self.graph = self.parent.graph.inducedSubGraph(node_list) #create the graph
        if self.name == self._default_name(): #update the name if not custom name
            self.name += " of size "+str(numnodes)
        
        self.graph.setName(self.name) #set name of the graph
        self.parent.subGraphs[self.graph.getId()] = self #add its id to the parent graph
        self.edgeDirection = self.parent.edgeDirection


#class to create a random connex graph        
class InducedConnexRandomGraph(AbstractGraph):
    def _default_name(self):
        return "Induced Random Connex Graph"

    def build_graph(self, numnodes):
        def get_random_neighbor(node):
            l = []
            for n in self.parent.graph.getInOutNodes(node):
                l+=[n]
            return random.choice(l)
        
        if self.parent == None:
            raise Exception("No parent graph specified")

        if numnodes > self.parent.nNodes():
            print "numnodes higher than total number of nodes!"
            return

        visited_nodes = { n : 0 for n in self.parent.graph.getNodes()}
        original_node = self.parent.graph.getRandomNode()
        node_list = [original_node]
        visited_nodes[original_node] = 1

        while len(node_list) < numnodes:
            n = get_random_neighbor(random.choice(node_list)) #get a random neighbor of one of the nodes from the current list chosen randomly 
            if visited_nodes[n] == 0:
                visited_nodes[n] = 1
                node_list+=[n]
            
        node_list = set(node_list)

        #now build the subgraph 
        self.graph = self.parent.graph.inducedSubGraph(node_list) #create the graph
        if self.name == self._default_name(): #update the name if not custom name
            self.name += " of size "+str(numnodes)
        
        self.graph.setName(self.name) #set name of the graph
        self.parent.subGraphs[self.graph.getId()] = self #add its id to the parent graph
        self.edgeDirection = self.parent.edgeDirection


#class which creates a subgraph containing the big connex subgraph of the citation graph        
class InducedConnexMaxGraph(AbstractGraph):
    def _default_name(self):
        return "Induced Connex Max Graph"

    def build_graph(self):

        def get_in_out_nodes(node):
            l = []
            for n in self.parent.graph.getInOutNodes(node):
                l+=[n]
            return l
        
        if self.parent == None:
            raise Exception("No parent graph specified")

        node_list = []
        original_node = tlp.node(25000) #taken randomly from the max connected part
        to_visit_nodes = [original_node]
        node_status = {n : 0 for n in self.parent.graph.getNodes()} #0 if not visited, 1 if visited or to be visited
        node_status[original_node] = 1
    
        while len(to_visit_nodes) != 0:
            
            node_list+=[to_visit_nodes[0]]
        
            l = get_in_out_nodes(to_visit_nodes[0])
        
            to_remove = []
            for i in l:
                if node_status[i] == 1:
                    to_remove+=[i]

            l = list(set(l) - set(to_remove)) #remove the nodes already visited or already in the to visit list
                
            for i in l:
                node_status[i] = 1#update status list
            
            to_visit_nodes = to_visit_nodes[1:] #remove first element
            to_visit_nodes += l #update the list of nodes to visit
                
        node_list = set(node_list) #convert to set
        
        #now build the subgraph 
        self.graph = self.parent.graph.inducedSubGraph(node_list) #create the graph
        if self.name == self._default_name(): #update the name if not custom name
            self.name += " of node "+self.parent.graph["node_property"][original_node]
        
        self.graph.setName(self.name) #set name of the graph
        self.parent.subGraphs[self.graph.getId()] = self #add its id to the parent graph
        self.edgeDirection = self.parent.edgeDirection

        
#creates a graph which contains the distinct paths connecting two nodes in the parent graph. Uses the maxflow relabel algorithm
class InducedConnectivityGraph(AbstractGraph):
    def _default_name(self):
        return "Induced Connectivity Graph"

    def build_graph(self, source, sink):
        if self.parent == None:
            raise Exception("No parent graph specified")

        connect_result = Connectivity(self.parent, source, sink)
        edge_list = connect_result[0]
        self.flow = connect_result[1]
        #now build the subgraph 
        if self.name == self._default_name(): #update the name if not custom name
            self.name += " from source "+str(source.id)+" to sink "+str(sink.id)
            
        self.graph = self.parent.graph.addSubGraph(self.name) #create the graph
        #add the edges to the graph
        for e in edge_list:
            self.graph.addNode(self.parent.graph.source(e))
            self.graph.addNode(self.parent.graph.target(e))
            self.graph.addEdge(e)
            
        self.parent.subGraphs[self.graph.getId()] = self #add its id to the parent graph
        self.edgeDirection = self.parent.edgeDirection


#Creates a graph identical to the parent graph but networkx compatible
class NetworkxGraph(AbstractGraph):
    def _default_name(self):
        return "Networkx compatible graph"
    
    def build_graph(self):
        if self.parent == None:
            raise Exception("No parent graph specified")
        
        #now build the subgraph 
        if self.name == self._default_name(): #update the name if not custom name
            self.name += " of "+self.parent.name
            
        self.graph = nx.DiGraph() #create the graph

        #fill in the nodes
        for n in self.parent.graph.getNodes():
            self.graph.add_node(n)

        #add the edges
        for e in self.parent.graph.getEdges():
            self.graph.add_edge(self.parent.graph.source(e),self.parent.graph.target(e))
            
        self.parent.subGraphs["nx"+str(self.parent.getId())] = self #add its id to the parent graph
        self.edgeDirection = self.parent.edgeDirection

    def getId(self): #Have to overload these functions
        return "nx"+str(self.parent.graph.getId())

    def nNodes(self): #get number of nodes
        return self.graph.number_of_nodes()

    def nEdges(self): #get number of edges
        return self.graph.number_of_edges()
    
    def numberOfSubGraphs(self): #get number of subgraphs
        pass

    ##Function to "transfer" the value of calculated coeffs to the parent graph
    def transferCoeff(self, H, coeffname):
        coeff = self.parent.graph.getDoubleProperty(coeffname) #by default doubleproperty, perhaps make it a choice by passing an argument?
        for n in H:
            coeff[n] = H[n]



#Create the partially constructed graph in the correct constructed order, with numnodes nodes
class PartiallyConstructedGraph(AbstractGraph):
    def _default_name(self):
        return "Partially Constructed Graph"

    def build_graph(self, numnodes, order = (0,0)):
        if self.parent == None:
            raise Exception("No parent graph specified")

        #now build the subgraph 
        if self.name == self._default_name(): #update the name if not custom name
            self.name += " with "+str(numnodes)+" nodes"
            
        self.graph = self.parent.graph.addSubGraph(self.name) #create the graph

        if order == (0,0):
            construction_order,edges = GetConstructionOrder(self.parent)
        else:
            construction_order,edges = order
        for i in range(numnodes): #for all remaining nodes
            #add the node and its out edges
            n = construction_order[i]
            self.graph.addNode(n)
            for t in edges[n]:
                e = self.parent.graph.existEdge(n,t) #have to do this otherwise it creates duplicate edges...
                self.graph.addEdge(e) 
            
            
        self.parent.subGraphs[self.graph.getId()] = self #add its id to the parent graph
        self.edgeDirection = self.parent.edgeDirection
