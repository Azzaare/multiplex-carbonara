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


Todo:


Low priority:
-faire marcher flowselectivenewnew

Todo:
-sortir tous les tableaux de mesure du temps
-essayer de fitter des courbes theoriques pour ces tableaux de mesure du temps
-sortir toutes les correlations pour tous les types de coeff ranks, pub ou aut, possibles
-determiner les plus grands "sauts" entre deux mesures V Essayer de plutot regarder les sauts en coefficients plutot que les sauts en classement

Todo sur le long terme:
-interpreter les coefficients

Observations : 
On va pas revolutioner le classement, si une pub a 1000 citations elle aura aussi un gros score

En revanche, si une mesure a peu de citations mais ses citations en ont beaucoup, alors elle fera remonter la publication


Observations : le probleme de FlowSelectionNewNew c'est qu'en fait les sous graphes incluent les sous graphes qui ne contiennent qu'un noeud. Y'a toujours une fraction du flot qui n'est pas envoye parce que c'est le flot correspondant au ss graphe induit par ce noeud qui n'est transmis a personne.

Donc ca peut pas marcher. Au pire on s'en fout. Mais si en fait ca devrait marcher si on met 0 cmome coeff dans les ss grapehs qui servent a rien. Bon, TODO eventuellement

Observations : le caractere unbounded on s'en fout pour le calcul du temps, par contre le cote k est interessant

Observation : peut etre que dans les progressions comparer le rank n'est pas la meilleure idee due a la facilite dont on augmente dans les ranks avec un tout petit increase de coeff

Observation : les fonctinos de mesures dynamiques ne servent pas a grand chose, sauf peut etre a effectuer des tests...


Les mesures de temps a calculer :
Pour commencer:
-degre (pour le lol)
-flot, kflot avec k variant en statique et en dynamique
-addnode, kaddnode avec k variant
-en gros faire varier les k. Voir comment ca impacte le temps de calcul.

En gros, le k va avoir un impact sur le temps de calcul, et le k et la divisibilite de l'information va avoir un impact sur les resultats (coeff de correlation)

Ensuite:
-les flots multiplexes en statique
-les flots multiplexes addnode 
-les flots multiplexes en dynamique 

Enfin:
-fitter des courbes polynomiales
-expliquer un peu les resultats obtenus
-comparer les valeurs THEORIQUES avec les valeurs experimentales

Pour terminer : 
-construire "le pire DAG possible" et tester les mesures la dessus


Observation : en fait nous on fait varier nos mesures flot selon 2 dimensions : la divisibilite de l'information, et la distance de propagation de l'information




TODO :
-pour tous les coefficients, ne garder que x chiffres significatifs


TODO WITH VERY HIGH PRIORITY :
debugger le bordel avec MeasureTime et les fonctions multiplexes


TODO WITH EXTREMELY HIGH PRIORITY : 
voir les correlations des percolations auteur pour les differentes mesures

TODO WITH SO MUCH HIGH PRIORITY IT FUCKS YOUR BRAIN :
WHY THE FUCK IS NOTHING WORKING? AGGREGATED DIVISION BY 0 WTF??
YE COUMPRI, completemet changer la maniere dont on mesure le time. Enlever samplemultiplexsubgraphs, faire juste samplesubgraphs. Et ensuite, a chaque iteration, recalculer le ss graphe

Ou faire une fonction "updatemultiplex"

Ouais faire ca en fait
Ce sera plus facile
et comme ca ya rien a changer

Bref se demerder

Pour l'instant stop la dessus. On se preoccupe d'autres choses
Faire le alpha Monoplex Flow. Tester les correlations
"""


def testtot():
    reloading()
    totalcorrelationmeasures()
    timemeasures()


def alphatest():
    reloading()
    g = root.addSubGraph(PubGraph)
    g = g.addSubGraph(InducedConnexMaxGraph)
    #g = root[1][2]
    #gaut = root[1][2][3]
    gaut = root[1].addSubGraph(AuthorGraph)

    root.printAllSubGraphs()
    
    me.AlphaMonoplexFlow_static(g,alpha = 0.2)
    me.MonoplexFlow_static(g)
    me.AlphaMonoplexFlow_static(g,alpha = 0.5)


    names = ["0.2FlowCoeff","0.5FlowCoeff","FlowCoeff"]
    
    rk.CompareMeasures(g,names)
    gf.PubCoeffsToAuthorCoeffs(g,gaut,names)
    rk.PrintTop20(gaut,[i+"PTASumPA" for i in names])

    
def timetest():
    reloading()
    g = root.addSubGraph(PubGraph)
    g = g.addSubGraph(InducedConnexRandomGraph, numnodes = 100)
    g = g.addSubGraph(MultiplexGraph)
    me.MultiplexFlowSum_static(g)
    me.Measure_dynamic(g, me.MultiplexFlowSum_addNode, coeffname="MultiplexFlowSumDynamicCoeff")

    rk.CompareMeasures(g,["MultiplexFlowSumCoeff","MultiplexFlowSumDynamicCoeff"])
    rk.PrintDifferences(g,"MultiplexFlowSumCoeff","MultiplexFlowSumDynamicCoeff")

def timemeasures():
    reloading()
    g = root.addSubGraph(PubGraph)
    g = g.addSubGraph(InducedConnexMaxGraph)
    gf.SamplePartiallyConstructedSubGraphs(g,3000,27000,1000) #have to start at 3000 because otherwise shit.

    
    #first the statics
    rk.MeasureTime_static(g, me.Degree_static,"TM_degreestatic")
    rk.MeasureTime_static(g, me.MonoplexFlow_static,"TM_monoplexflowstatic")
    """rk.MeasureTime_static(g, me.MultiplexFlowSum_static,"TM_mxflowsumstatic")
    rk.MeasureTime_static(g, me.MultiplexFlowAggregated_static,"TM_mxflowaggstatic")
    rk.MeasureTime_static(g, me.MultiplexFlowSelective_static,"TM_mxflowselstatic")"""

    for k in [1,2,5,10,20]:
        KMonoplexFlow_static = lambda g, depth = k, coeffname="FlowCoeff" : me.KMonoplexFlow_static(g,depth,coeffname)
        #KMultiplexFlowAggregated_static = lambda g, depth = k, coeffname="MultiplexFlowAggregatedCoeff" : me.KMultiplexFlowAggregated_static(g,depth,coeffname)
        rk.MeasureTime_static(g, KMonoplexFlow_static,"TM_"+str(k)+"monoplexflowstatic")
        #rk.MeasureTime_static(g, KMultiplexFlowAggregated_static,"TM_"+str(k)+"mxflowaggstatic")

    
    #now the addnodes
    rk.MeasureTime_addNode(g, me.Degree_static, me.Degree_addNode, "TM_degreeaddnode")
    rk.MeasureTime_addNode(g, me.MonoplexFlow_static,me.MonoplexFlow_addNode,"TM_monoplexflowaddnode")
    """rk.MeasureTime_addNode(g, me.MultiplexFlowSum_static,me.MultiplexFlowSum_addNode,"TM_mxflowsumaddnode")
    rk.MeasureTime_addNode(g, me.MultiplexFlowAggregated_static,me.MultiplexFlowAggregated_addNode,"TM_mxflowaggaddnode")
    rk.MeasureTime_addNode(g, me.MultiplexFlowSelective_static,me.MultiplexFlowSelective_addNode,"TM_mxflowseladdnode")"""
    
    for k in [1,2,5,10,20]:
        KMonoplexFlow_static = lambda g, depth = k, coeffname="FlowCoeff" : me.KMonoplexFlow_static(g,depth,coeffname)
        #KMultiplexFlowAggregated_static = lambda g, depth = k, coeffname="MultiplexFlowAggregatedCoeff" : me.KMultiplexFlowAggregated_static(g,depth,coeffname)
        KMonoplexFlow_addNode = lambda g, node, depth = k, coeffname="FlowCoeff" : me.MonoplexFlow_addNode(g,node,depth,coeffname)
        #KMultiplexFlowAggregated_addNode = lambda g, node, depth = k, coeffname="MultiplexFlowAggregatedCoeff" : me.MultiplexFlowAggregated_addNode(g,node,depth,coeffname)
        rk.MeasureTime_addNode(g, KMonoplexFlow_static, KMonoplexFlow_addNode,"TM_"+str(k)+"monoplexflowaddnode")
        #rk.MeasureTime_addNode(g, KMultiplexFlowAggregated_static, KMultiplexFlowAggregated_addNode,"TM_"+str(k)+"mxflowaggaddnode")


    #now the dynamics, threshold = 60 seconds
    rk.MeasureTime_dynamic(g, me.Degree_addNode,"TM_degreedynamic",threshold=60.0)
    rk.MeasureTime_dynamic(g, me.MonoplexFlow_addNode,"TM_monoplexflowdynamic",threshold=60.0)
    """rk.MeasureTime_dynamic(g, me.MultiplexFlowSum_addNode,"TM_mxflowsumdynamic",threshold=60.0)
    rk.MeasureTime_dynamic(g, me.MultiplexFlowAggregated_addNode,"TM_mxflowaggdynamic",threshold=60.0)
    rk.MeasureTime_dynamic(g, me.MultiplexFlowSelective_addNode,"TM_mxflowseldynamic",threshold=60.0)"""

    for k in [1,2,5,10,20]:
        KMonoplexFlow_addNode = lambda g, node, depth = k, coeffname="FlowCoeff" : me.MonoplexFlow_addNode(g,node,depth,coeffname)
        #KMultiplexFlowAggregated_addNode = lambda g, node, depth = k, coeffname="MultiplexFlowAggregatedCoeff" : me.MultiplexFlowAggregated_addNode(g,node,depth,coeffname)
        rk.MeasureTime_dynamic(g, KMonoplexFlow_addNode,"TM_"+str(k)+"monoplexflowdynamic",threshold=60.0)
        #rk.MeasureTime_dynamic(g, KMultiplexFlowAggregated_addNode,"TM_"+str(k)+"mxflowaggdynamic",threshold=60.0)"""

    
   
    
def totalcorrelationmeasures(): 
    reloading()
    """g = root.addSubGraph(PubGraph)
    g = g.addSubGraph(InducedConnexMaxGraph)
    g = g.addSubGraph(MultiplexGraph)

    me.HIndexPub_static(g)
    me.ExForce_static(g)

    me.Degree_static(g, coeffname = "1DegreeCoeff")
    me.DagDepth_static(g)
    me.KDegree_static(g, coeffname = "InfDegreeCoeff")
    
    
    me.MonoplexFlow_static(g)
    me.MonoplexFlowUnbounded_static(g)
    
    me.MultiplexFlowAggregated_static(g)
    me.MultiplexFlowAggregatedUnbounded_static(g)
    me.MultiplexFlowSum_static(g)
    me.MultiplexFlowSumUnbounded_static(g)
    me.MultiplexFlowSelective_static(g)
    me.MultiplexFlowSelectiveUnbounded_static(g)

    for k in [2,5,10,20]:
        me.KDegree_static(g, depth = k, coeffname = str(k)+"DegreeCoeff")
    for k in [1,2,5,10,20]:
        me.KMonoplexFlow_static(g, depth = k, coeffname = str(k)+"FlowCoeff")
        me.KMonoplexFlowUnbounded_static(g, depth=k, coeffname = str(k)+"FlowUnboundedCoeff")
        me.KMultiplexFlowAggregated_static(g, depth = k, coeffname = str(k)+"MultiplexFlowAggregatedCoeff")
        me.KMultiplexFlowAggregatedUnbounded_static(g, depth=k, coeffname = str(k)+"MultiplexFlowAggregatedUnboundedCoeff")
        
    for a in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]:
        me.AlphaMonoplexFlow_static(g,a)
        me.AlphaMonoplexFlowUnbounded_static(g,a)
        me.AlphaMultiplexFlowAggregated_static(g,a)
        me.AlphaMultiplexFlowAggregatedUnbounded_static(g,a)

    tlp.saveGraph(g.graph,"bestgraph.tlp")"""

    #set the list of names
    g=root[1][2][3]
    
    names = ["FlowCoeff","FlowUnboundedCoeff","MultiplexFlowSumCoeff","MultiplexFlowSumUnboundedCoeff","MultiplexFlowAggregatedCoeff","MultiplexFlowAggregatedUnboundedCoeff","MultiplexFlowSelectiveCoeff","MultiplexFlowSelectiveUnboundedCoeff","HIndexPubCoeff","ExForceCoeff","DagDepthCoeff","InfDegreeCoeff"]

    names+= [str(k)+"DegreeCoeff" for k in [1,2,5,10,20]]
    names+= [str(k)+"FlowCoeff" for k in [1,2,5,10,20]]
    names+= [str(k)+"FlowUnboundedCoeff" for k in [1,2,5,10,20]]
    names+= [str(k)+"MultiplexFlowAggregatedCoeff" for k in [1,2,5,10,20]]
    names+= [str(k)+"MultiplexFlowAggregatedUnboundedCoeff" for k in [1,2,5,10,20]]

    names+= [str(a)+"FlowCoeff" for a in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]]
    names+= [str(a)+"FlowUnboundedCoeff" for a in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]]
    names+= [str(a)+"MultiplexFlowAggregatedCoeff" for a in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]]
    names+= [str(a)+"MultiplexFlowAggregatedUnboundedCoeff" for a in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]]

    #export the pub coeffs
    #only use spearman coeffs
    rk.ExportCoeffs_ordered(g,names,"coeffspubs")

    #now transfer this to an author graph and export all the coeffs
    gaut = root[1][2].addSubGraph(AuthorGraph)
    #hindex
    me.HIndex_static(gaut)
    gf.PubCoeffsToAuthorCoeffs(g,gaut,names,ReductionFunc = gf.PTASum)
    gf.PubCoeffsToAuthorCoeffs(g,gaut,names,ReductionFunc = gf.PTASumPA)
    gf.PubCoeffsToAuthorCoeffs(g,gaut,names,ReductionFunc = gf.PTAAverage)
    gf.PubCoeffsToAuthorCoeffs(g,gaut,names,ReductionFunc = gf.PTAAveragePA)

    #still use only spearman
    autnames = [n+"PTASum" for n in names]+[n+"PTASumPA" for n in names]+[n+"PTAAverage" for n in names]+[n+"PTAAveragePA" for n in names]
    rk.ExportCoeffs_ordered(gaut,autnames,"coeffsauts")
    


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


