from tulip import *
#import citation_graph as cg
from graphclass import *
import graphfuncs as gf
import measures as me
import rankings as rk
import numpy as np
import networkx as nx


"""
#TODO
Fusion noeuds seulement si memes auteurs, et dans ce cas faire une vraie fusion. Sinon on supprime les aretes tout simplement


"""


#determiner mathematiquement comment le sampling random permet d'approximer les mesures reeles
#determiner une "g-mesure" grace a l'analyse factorielle.
#HYPOTHESE : il existe une "vraie" mesure de l'influence, la g-mesure, correlee avec nos mesures calculees.
#pratique : mesurer en pratique le nbe de calculs necessaires en incremental et comparer au predictions theoriques
#determiner des caracteristiques structurelles du graphe genre maxdepth, mindeg maxdeg, averagedeg, etc.
#FAIRE DU MULTIPLEX. Caracteriser les sous-couches multiplex. Les caracteriser qualitativement (qu'est ce que ca veut dire?) mais aussi les caracteriser quantitativement (profondeur du dag, etc...)
#Mesure multiplex : somme des flots. Encore mieux : flot pondere par l'inverse du nombre de citations communes
#IMPORTANT attribuer a chaque arrete les neouds dont elle appartient au ss dag induit de depth infini IMPORTANT

#TODO :faire en sorte quon puisse reprendre les calculs de mesure en cours
#TODO : utiliser les bibliotheques numpy pour calculer les corrcoeff
#TODO : reflechir sur comment utiliser les edges dans le calcul de mesure

#TODO : h-index resultats interessants. Tester pour voir si ces resultats sont accurate

#TODO : essayer de voir si on ne peut pas definir un mesure en termes "dynamiques", cad en analysant l'impact de l'ajout d'un noeud sur la mesure

#Ameliorer InducedConnectivityGraph en ne prenant deja que le ss dag induit par la source, et en commanceant a verifier si le node n'est pas inclus dans le ss dag. Dailleurs il doit y avoir un moyen beaucoup plus efficace que le maxflow relabel vu quon a un dag

#Creer une fonction print_metadata qui imprime correctement le metadata
#Et puis ajouter le nom des publications au nodes publications!

#PROBLEME JUAN MALCADENA


#Dire qu'il y a fondamentalement deux problemes avec la mesure flot monoplexe : on fait tout remonter jusqu'en haut, cad qu'on ne peut pas determiner a partir de quelle profondeur l'influence d'une publication devient negligeable, et deuxieme probleme on ne peut pas determiner, lorsque plusieurs publications sont citees, laquelle est plus importante vis a vis de l'autre. C'est a dire si une publication cite 5 publications, comment determiner laquelle de ces publications participe plus au transfert de l'information? En fait le probleme fondamental c'est la nature de graphe, c'est a dire qu'en ne considerant que la monoplexite on ne peut pas apprecier des differences entre deux aretes "equivalentes" 




#probleme de calcul flottant, ca donne des valeurs legerement differentes...


"""

Bugs 
-multiplex selectif avait nu bug, modifier le multiplex selectif unbuonded

Todo:

Medium priority :
-faire les fonctions de mesure de temps
-sortir toutes les mesures de temps
-faire version dynamique des mesures multiplex
-faire version dynamique de toutes les mesures en fait

Low priority:
-determiner a partir de combien de noeuds a ajouter dynamiquement cela devient il plus rentable de les ajouter statiquement
-tracer les courbes tpsdynamic/tpsstatic
-faire marcher flowselectivenewnew

Todo:
-sortir tous les tableaux de mesure du temps
-essayer de fitter des courbes theoriques pour ces tableaux de mesure du temps
-sortir toutes les correlations pour tous les types de coeff ranks, pub ou aut, possibles
-reecrire rapport un peu mieux
-changer l'introduction pour vraiment montrer le fait qu'on va pas revolutioner les classements et dire que le but c'est vraiment de voir quelles publications ressortent plus avec nos mesures, pas dire que ed witten est le meilleur
-determiner les plus grands "sauts" entre deux mesures

Todo sur le long terme:
-interpreter les coefficients

Observations : 
On va pas revolutioner le classement, si une pub a 1000 citations elle aura aussi un gros score

En revanche, si une mesure a peu de citations mais ses citations en ont beaucoup, alors elle fera remonter la publication


Observations : le probleme de FlowSelectionNewNew c'est qu'en fait les sous graphes incluent les sous graphes qui ne contiennent qu'un noeud. Y'a toujours une fraction du flot qui n'est pas envoye parce que c'est le flot correspondant au ss graphe induit par ce noeud qui n'est transmis a personne.

Donc ca peut pas marcher. Au pire on s'en fout. Mais si en fait ca devrait marcher si on met 0 cmome coeff dans les ss grapehs qui servent a rien. Bon, TODO eventuellement

"""

def testselective():
    reloading()

    #g = root.addSubGraph(PubGraph)
    #g = g.addSubGraph(InducedConnexMaxGraph)
    #g = g.addSubGraph(MultiplexGraph)
    g = root[1][2][3]
    
    #me.MultiplexFlowSelectiveNew_static(g, coeffname = "MultiplexFlowSelectiveNewCoeff")
    me.MultiplexFlowSelectiveNewNew_static(g, coeffname = "MultiplexFlowSelectiveNewNewCoeff")

    #rk.CompareMeasures(g,["MultiplexFlowSelectiveNewCoeff","MultiplexFlowSelectiveNewNewCoeff"])
    #rk.PrintTop20(g,["MultiplexFlowSelectiveNewCoeff","MultiplexFlowSelectiveNewNewCoeff"])
    for n in g.graph.getNodes():
        print g.graph["MultiplexFlowSelectiveNewNewCoeff"][n]
        

def testtruc():
    reloading()
    g = root.addSubGraph(PubGraph)
    g = g.addSubGraph(InducedConnexMaxGraph)
    gf.SamplePartiallyConstructedMultiplexSubGraphs(g,100,27000,1000)
    #g = root[27384][27385]
    #rk.MeasureTime_static(g, me.Degree_static,"degree")
    rk.MeasureTime_static(g, me.MonoplexFlow_static,"flow")
    rk.MeasureTime_dynamic(g, me.MonoplexFlow_static, me.MonoplexFlow_addNode, "flowaddnode2")
    #rk.MeasureTime_static(g, me.MultiplexFlowSum_static,"mxflowsum")
    #rk.MeasureTime_static(g, me.MultiplexFlowAggregated_static,"mxflowaggregated")
    #rk.MeasureTime_static(g, MonoplexFlow_static,"mxflowselective")
    
def testhindex():
    reloading()
    g = root.addSubGraph(PubAuthorGraph)
    me.HIndexPubAut_static(g)
    l = rk.RankAuthorsByCoeff(g, g.graph["HIndexCoeff"])
    for i in range(20):
        print g.graph["node_property"][l[i]]

    gr = root.addSubGraph(AuthorGraph)
    gf.TransferCoeff(g,gr,"HIndexCoeff")
    l = rk.RankAuthorsByCoeff(gr, gr.graph["HIndexCoeff"])
    for i in range(20):
        print gr.graph["node_property"][l[i]]

def testauthorcoeffs():
    reloading()
    gaut = root[1][2].addSubGraph(AuthorGraph)
    g = root[1][2][3]
    rk.PubCoeffsToAuthorCoeffs(g,gaut,["DegreeCoeff","FlowCoeff", "MultiplexFlowSumCoeff","MultiplexFlowAggregatedCoeff","MultiplexFlowSelectiveCoeff"], rk.PubToAuthorPonderedSum)

    rk.AuthorCoeffToRank(gaut, ["authorDegreeCoeff","authorFlowCoeff", "authorMultiplexFlowSumCoeff","authorMultiplexFlowAggregatedCoeff","authorMultiplexFlowSelectiveCoeff"])

    tlp.saveGraph(gaut.graph,"authorgraph.tlp")
    
    
def correlations():
    reloading()
    g = root[1][2][3]
    """g = root.addSubGraph(PubGraph)
    g = g.addSubGraph(InducedConnexMaxGraph)
    g = g.addSubGraph(MultiplexGraph)

    me.MonoplexFlow_static(g)
    me.MultiplexFlowAggregated_static(g)
    me.MultiplexFlowSum_static(g)
    me.MultiplexFlowSelective_static(g)
    me.Degree_static(g)

    names = ["DegreeCoeff","FlowCoeff","MultiplexFlowSumCoeff","MultiplexFlowAggregatedCoeff","MultiplexFlowSelectiveCoeff"]
    
    rk.ExportCoeffs_ordered(g, names, rk.RankPubsByCoeff, rk.KendallTauDistance, "coeffskendall")

    rk.CoeffToRank(g, g.graph["FlowCoeff"], "FlowCoeffRank")
    rk.CoeffToRank(g, g.graph["DegreeCoeff"], "DegreeCoeffRank")
    rk.CoeffToRank(g, g.graph["MultiplexFlowSumCoeff"], "MultiplexFlowSumCoeffRank")
    rk.CoeffToRank(g, g.graph["MultiplexFlowAggregatedCoeff"], "MultiplexFlowAggregatedCoeffRank")
    rk.CoeffToRank(g, g.graph["MultiplexFlowSelectiveCoeff"], "MultiplexFlowSelectiveCoeffRank")"""
    rk.PubCoeffToRank(g, ["FlowCoeff","DegreeCoeff","MultiplexFlowSumCoeff","MultiplexFlowAggregatedCoeff","MultiplexFlowSelectiveCoeff"])
    
    #l = rk.BestAuthorProgressionRank(g, g.graph["DegreeCoeff"],g.graph["FlowCoeff"])
    #for i in range(20):
     #   print i, " : ", root.graph["node_property"][l[i]]
    tlp.saveGraph(g.graph, "rapportgraph.tlp")

def newselective():
    reloading()
    """g = root.addSubGraph(PubGraph)
    g = g.addSubGraph(MultiplexGraph)
    me.MonoplexFlow_static(g)
    me.MultiplexFlowSelective_static(g)
    me.MultiplexFlowSelectiveNew_static(g)"""
    g = root[3][4]
    
    rk.CompareMeasures(g,["FlowCoeff","MultiplexFlowSelectiveCoeff","MultiplexFlowSelectiveNewCoeff"],rk.SpearmanRankCoeff)
    rk.PrintTopAuthors(g,["FlowCoeff","MultiplexFlowSelectiveCoeff","MultiplexFlowSelectiveNewCoeff"],rk.PubToAuthorPonderedSum)

def inducedcoefftest():
    reloading()
    H = rk.InducedCoeffs(["FlowCoeff","5FlowCoeff","20FlowCoeff"], "corcoef","newcorcef")
    print H

def timetest():
    reloading()
    #g = root[27776]
    g = root.addSubGraph(PubGraph)
    g = g.addSubGraph(InducedConnexMaxGraph)
    print "plop"
    gf.SamplePartiallyConstructedSubGraphs(g, 100, 5000, 500)
    g.printSubGraphs()
    rk.MeasureTime_static(g, me.MonoplexFlow_dynamic, "testtime")

def partialgraphtest():
    reloading()
    g = root.addSubGraph(PubGraph)
    print "plop"
    g = g.addSubGraph(PartiallyConstructedGraph, numnodes = 15000)
    return g
    
#An example of how to calculate the measures:

def timedynamictest():
    reloading()
    g = root.addSubGraph(PubGraph)
    g = g.addSubGraph(InducedConnexMaxGraph)
    g = g.addSubGraph(InducedConnexRandomGraph, numnodes = 15000)
    rk.MeasureTime_dynamic(g, me.MonoplexFlow_static, "prout")

def test():
    reloading()
    g = root.addSubGraph(PubGraph) #add a publication graph
    g = g.addSubGraph(InducedConnexMaxGraph) #get the big connex graph contained in the publication graph
    g = g.addSubGraph(MultiplexGraph)#demultiplex the graph
    me.MonoplexFlow_static(g) #calculate the monoplex flow measure
    me.Degree_static(g) #calculate the degree measure
    me.MultiplexFlowSum_static(g) #calculate the multiplex flow sum measure
    me.MultiplexFlowAggregated_static(g)#calculate the multiplex flow aggregated measure
    me.MultiplexFlowSelective_static(g)#calculate the multiplex flow selective measure

    #get the names of the coeffs
    names = ["FlowCoeff","DegreeCoeff","MultiplexFlowSumCoeff","MultiplexFlowAggregatedCoeff","MultiplexFlowSelectiveCoeff"]
    #print the top authors of the measures, using a reduction function which divides the coeff of a pub and distributes it among its authors.
    rk.PrintTopAuthors(g,names,rk.PubToAuthorPonderedSum)
    #calculate the spearman correlations of the measures and export them to a file
    rk.ExportCoeffs_ordered(g, names, rk.RankPubsByCoeff, rk.SpearmanRankCoeff, "test")
    

def rapportCitationGraph():
    reloading()
    n = 10468 # node id : 2123
    g = root.addSubGraph(PubGraph)
    g = g.addSubGraph(InducedDagGraph, node = tlp.node(n), depth = 2, direction = "reversed")

    g = g.addSubGraph(MultiplexGraph)
        
    
    """pub = root[1].addSubGraph(InducedDagGraph, node = tlp.node(n), depth = 2, direction = "reversed")
    pubaut = root[1].addSubGraph(InducedDagGraph, node = tlp.node(n), depth = 2, direction = "reversed")"""
    tlp.saveGraph(g.graph, "rapport.tlp")
    


def reloading():
    reload(gf)
    reload(me)
    reload(rk)
    print "reloaded"


