# Powered by Python 2.7

# To cancel the modifications performed by the script
# on the current graph, click on the undo button.

# Some useful keyboards shortcuts : 
#   * Ctrl + D : comment selected lines.
#   * Ctrl + Shift + D  : uncomment selected lines.
#   * Ctrl + I : indent selected lines.
#   * Ctrl + Shift + I  : unindent selected lines.
#   * Ctrl + Return  : run script.
#   * Ctrl + F  : find selected text.
#   * Ctrl + R  : replace selected text.
#   * Ctrl + Space  : show auto-completion dialog.

from tulip import *

# the updateVisualization(centerViews = True) function can be called
# during script execution to update the opened views

# the pauseScript() function can be called to pause the script execution.
# To resume the script execution, you will have to click on the "Run script " button.

# the runGraphScript(scriptFile, graph) function can be called to launch another edited script on a tlp.Graph object.
# The scriptFile parameter defines the script name to call (in the form [a-zA-Z0-9_]+.py)

# the main(graph) function must be defined 
# to run the script on the current graph

#Function to get the source nodes of a graph (i.e. nodes with no incoming edges)
def get_source_nodes(graph):
	l = []	
	for n in graph.getNodes():
		if graph.indeg(n) == 0 :
			l+=[n]
	return l

#Function to reverse all edges in a graph
def reverse_graph(graph):
	for e in graph.getEdges():	
		graph.reverse(e)

#Function to delete all subgraphs of a graph
def del_sub_graphs(graph):
	for i in graph.getSubGraphs():
		graph.delSubGraph(i)

#Function to calculate the "direct flow" in a graph
def flow_direct(graph,metric):
	if graph.numberOfNodes() == 0:	
		return 
	source_list = get_source_nodes(graph)
	for n in source_list:
		for m in graph.getOutNodes(n):
			metric.setNodeValue(m,(metric.getNodeValue(m)+metric.getNodeValue(n)/graph.outdeg(n)))
		for e in graph.getOutEdges(n):
			metric.setEdgeValue(e,(metric.getNodeValue(n)/graph.outdeg(n)))
		graph.delNode(n)	
	return flow_direct(graph, metric)

#Function to calculate the "reversed flow" in a graph
def flow_reversed(graph,metric):
	reverse_graph(graph)
	return flow_direct(graph,metric)
	
#Function to display the graph taking into account a metric (coeff_direct or coeff_reversed)
def display_graph(graph,metric):
	viewSize = graph.getSizeProperty("viewSize")	
	viewLabel = graph.getStringProperty("viewLabel")
	viewColor = graph.getColorProperty("viewColor")
	
	dataset = tlp.getDefaultPluginParameters("Size Mapping")
	dataset['property'] = metric
	dataset['max size'] = 5
	dataset['min size'] = 0.01
	
	graph.applySizeAlgorithm('Size Mapping', viewSize, dataset)

	dataset = tlp.getDefaultPluginParameters('Color Mapping')
	dataset['input property'] = metric
	graph.applyColorAlgorithm('Color Mapping', graph["viewColor"], dataset)
	
	#for it to be less ugly..
	baseSizeEdge = tlp.Size(0.125,0.125,0.5	)
	for n in graph.getEdges():
		viewSize[n] = baseSizeEdge
		viewLabel[n] = ""
		viewColor[n] = tlp.Color.Byzantium
	
#Function to reset the display of the graph
def reset_display_graph(graph):
	viewSize = graph.getSizeProperty("viewSize")	
	viewLabel = graph.getStringProperty("viewLabel")
	viewColor = graph.getColorProperty("viewColor")
	baseSizeNode = tlp.Size(1,1,1)
	baseSizeEdge = tlp.Size(0.125,0.125,0.5	)

	for n in graph.getNodes():
		viewSize[n] = baseSizeNode
		viewLabel[n] = ""
		viewColor[n] = tlp.Color.Red
	for n in graph.getEdges():
		viewSize[n] = baseSizeEdge
		viewLabel[n] = ""
		viewColor[n] = tlp.Color.Gray
	

def main(graph): 
	#ID = graph.getStringProperty("ID") 
	viewBorderColor = graph.getColorProperty("viewBorderColor")
	viewBorderWidth = graph.getDoubleProperty("viewBorderWidth")
	viewColor = graph.getColorProperty("viewColor")
	viewFont = graph.getStringProperty("viewFont")
	viewFontAwesomeIcon = graph.getStringProperty("viewFontAwesomeIcon")
	viewFontSize = graph.getIntegerProperty("viewFontSize")
	viewLabel = graph.getStringProperty("viewLabel")
	viewLabelBorderColor = graph.getColorProperty("viewLabelBorderColor")
	viewLabelBorderWidth = graph.getDoubleProperty("viewLabelBorderWidth")
	viewLabelColor = graph.getColorProperty("viewLabelColor")
	viewLabelPosition = graph.getIntegerProperty("viewLabelPosition")
	viewLayout = graph.getLayoutProperty("viewLayout")
	viewMetaGraph = graph.getGraphProperty("viewMetaGraph")
	viewMetric = graph.getDoubleProperty("viewMetric")
	viewRotation = graph.getDoubleProperty("viewRotation")
	viewSelection = graph.getBooleanProperty("viewSelection")
	viewShape = graph.getIntegerProperty("viewShape")
	viewSize = graph.getSizeProperty("viewSize")
	viewSrcAnchorShape = graph.getIntegerProperty("viewSrcAnchorShape")
	viewSrcAnchorSize = graph.getSizeProperty("viewSrcAnchorSize")
	viewTexture = graph.getStringProperty("viewTexture")
	viewTgtAnchorShape = graph.getIntegerProperty("viewTgtAnchorShape")
	viewTgtAnchorSize = graph.getSizeProperty("viewTgtAnchorSize")
	
        #reset the display
	reset_display_graph(graph)
	
        #clone graphs, as to not delete nodes in the real graph
        subgraph_reversed = graph.addCloneSubGraph()
	subgraph_direct = graph.addCloneSubGraph()
	
        #get metrics
	coeff_direct = graph.getDoubleProperty("coeff_direct")
	coeff_reversed = graph.getDoubleProperty("coeff_reversed")
	
        #set all the initial metric values to 1
	coeff_direct.setAllNodeValue(1.0)	
	coeff_reversed.setAllNodeValue(1.0)	

        #calculate flow
	flow_direct(subgraph_direct, coeff_direct)
	flow_reversed(subgraph_reversed, coeff_reversed)
	
        #reverse graph back to original
	reverse_graph(graph)

        #display graph
	display_graph(graph,coeff_direct)	
	
        print "Done"
	
        #delete clone subgraphs
	del_sub_graphs(graph)
	
	#regarder plus grand chemin entre 2 sommets
	#determiner quels ss ensembles utiliser
