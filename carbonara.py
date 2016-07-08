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


import time                                                

def timeit(method):

    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()

        print '%r (%r, %r) %2.2f sec' % \
              (method.__name__, args, kw, te-ts)
        return result

    return timed



#An example of how to calculate the measures:
@timeit
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
