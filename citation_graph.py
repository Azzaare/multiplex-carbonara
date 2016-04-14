import os
from py2neo import Graph, Node, Relationship, authenticate

#lists and tables required to parse the date

months = {
        "jan": 1,
        "january": 1,
        "feb": 2,
        "february": 2,
        "mar": 3,
        "march": 3,
        "apr": 4,
        "april": 4,
        "may": 5,
        "jun": 6,
        "june": 6,
        "jul": 7,
        "july": 7,
        "aug": 8,
        "august": 8,
        "sep": 9,
        "september": 9,
        "oct": 10,
        "october": 10,
        "nov": 11,
        "november": 11,
        "dec": 12,
        "december":12
}

days = ["mon","tue","wed","thu","fri","sat","sun"]

dates = ["1","01","2","02","3","03","4","04","5","05","6","06","7","07","8","08","9","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31"]

years = ["1991","1992","1993","1994","1995","1996","1997","1998","1999","2000","2001","2002","2003"]

years_short = {"91":1991,"92":1992,"93":1993,"94":1994,"95":1995,"96":1996,"97":1997,"98":1998,"99":1999,"00":2000,"01":2001,"02":2002,"03":2003}


#function used in parsing authors list
def remove_text_inside_brackets(text, brackets="()[]"): #taken from http://stackoverflow.com/questions/14596884/remove-text-between-and-in-python
        count = [0] * (len(brackets) // 2) # count open/close brackets
	saved_chars = []
	for character in text:
		for i, b in enumerate(brackets):
			if character == b: # found bracket
				kind, is_close = divmod(i, 2)
				count[kind] += (-1)**is_close # `+1`: open, `-1`: close
				if count[kind] < 0: # unbalanced bracket
					count[kind] = 0
				break
		else: # character is not a bracket
			if not any(count): # outside brackets
				saved_chars.append(character)
	return ''.join(saved_chars)

#function used to determine the publication at which to start push_neo_graph() in function of the total number of citations already loaded
def citation_no(pub_data,l):
        k=0
        h=0
        for i in pub_data:
                for j in pub_data[i][0]:
                        if h == l:
                                return k
                        h=h+1
                k=k+1
                                
#Parsing functions for parsing the date and the author list

def parse_date(line): 
	l =" ".join(line.split()) #remove extra spaces
	l = l.lower() #remove capitals
	l = l.split(' ') #split the sentence into words
	
	j=0
	m=0
	a=0

	for i in l :
		if i in dates:
			j = int(i)
		if i in months:
			m = months[i]
		if i in years :
			a = int(i)
		if i in years_short:
			a = years_short[i]

	return [j,m,a]




def parse_author(line): # Can be better 
	l = line.strip() #remove special chars
	l = remove_text_inside_brackets(l)
        #remove all instances of special accents (\\)
	l = l.replace("\\'","") 
	l = l.replace("\\\"","")
	#l = l.replace("\\","") there are still special accents to remove
	l =" ".join(l.split()) #remove extra spaces
	l = l.split(' ',2) #delete the "Authors:"
	l = " ".join(l[1:])
	
	l = l.split('and ') #remove the "and"s and commas
	lp = []
	for i in l:
		lp += i.split(',')
	lp = [x.lstrip().rstrip() for x in lp if x.lstrip().rstrip() != ""] #remove the spaces at the beginning and end of authors name
	return lp
    
#Function for loading the data structure which associates for each publication the other publications which it cites, its publication date and its list of authors

def load_data():
        pub_data = {} #Data structure for our program. Associates to an id (int) a list of 3 lists : the list of citations, the date and the list of authors 
        print "Loading data..."

        #First we will load the file with the citation data to add the citations to the data structure

        f = open('/home/vivek/prog/multiplex-carbonara/Cit-HepTh.txt','r') 
        for i in range(4): #first four lines are useless
                line = f.readline() 
        
        for line in f : #read lines
                l = line.strip().split('\t') 

                i1 = int(l[0])
                if i1 not in pub_data:
                        pub_data[i1] = [[],[],[]] #if the entry for that publication doesn't exit, initialize it
                i2 = int(l[1])
                if i2 not in pub_data:
                        pub_data[i2] = [[],[],[]] #if the entry for that publication doesn't exit, initialize it

                pub_data[i1][0].append(i2) #add citation 

        #Secondly we will load the files with the metadata to add the dates and authors of the publications to the data structure

        for root,dirs,fns in os.walk("/home/vivek/prog/multiplex-carbonara/cit-HepTh-abstracts/") :
                for fn in fns :
                        if fn.endswith(".abs") :	
                                f = open(os.path.join(root, fn),'r')
                                id = int(fn.split('.')[0]) #the ID of the publication is its filename
                                if id in pub_data: #if the publication is in our citations data
                                        lauthors = [] #list of authors for the publication
                                        ldate = [] #date for the publication, in the format [day,month,year] (int)
    
                                        line=f.readline()
                                        while line != "" :
                                                if line.split(' ')[0] == "Date:" :
                                                        ldate=parse_date(line)
                                                if line.split(' ')[0] == "Authors:" or line.split(' ')[0] == "Author:" : #Authors can be written over several lines...
                                                        laut = line
                                                        line = f.readline()
                                                        while (line.split(' ')[0] != "Comments:" and line.split(' ')[0] != "Report-no:" and
                                                        line.split(' ')[0] != "Subj-class:" and line.split(' ')[0] != "Journal-ref:" and 
                                                        line.split(' ')[0].strip() != "\\\\") : #we read until we reach another section
                                                                laut+=line
                                                                line = f.readline()
                                                        lauthors = parse_author(laut)
                                                line = f.readline()
                                
                                                pub_data[id][1] = ldate #add the metadata to the data structure
                                                pub_data[id][2] = lauthors 
        print "Data loaded"
        return pub_data


#Function which creates a graph in neo4j from the data in pub_data
""" 
Structure of the graph
Nodes : authors, publications, dates
Relations : CITED_BY (between two pubs) WROTE (between pubs and authors) PUB_DAY/MONTH/YEAR (between pubs and days/months/years)
"""

def push_neo_graph(pub_data):
        #connect to the server
        authenticate("localhost:7474", "neo4j", "azerty01") #put your login and password here
        graph = Graph("http://localhost:7474/db/data/")
        
        print "Checking current graph status"
        #checker si les dates/nodes/auteurs ont ete crees
        #checker ou on en est avec les relations
        
        
        record = graph.cypher.execute("MATCH ()-[r:CITED_BY]->() RETURN count(r) AS count")
        l = record[0].count
        
        if l>=352807 : 
                print "Graph is already loaded"
                return
        
        if l == 0 : #if there is no graph, we initialize the other nodes
                #starting to create graph
                del_neo_graph() #delete whatever already exists
                print "Starting graph creation"
                tx = graph.cypher.begin() #tx will be our cypher query object throughout this function

                #add nodes for the dates
                print "Creating date nodes..."
                for i in range(1,32):
                        tx.append("CREATE (n:Day {name:{i}})",{"i": i})
                for i in range(1,13):
                        tx.append("CREATE (n:Month {name:{i}})",{"i": i})
                for i in range(13):
                        tx.append("CREATE (n:Year {name:{i}})",{"i": 1991+i})
                tx.commit()
                print "Done"

                #add all the publication nodes and connect them to the date nodes
                print "Creating publication nodes..."
                tx = graph.cypher.begin()
                for i in pub_data:
                        tx.append("MATCH (d:Day), (m:Month), (y:Year) "
                                  "WHERE d.name={d} AND m.name = {m} AND y.name = {y} "
                                  "CREATE (n:Publication {id:{i}}), (n)-[:PUB_DAY]->(d), (n)-[:PUB_MONTH]->(m), (n)-[:PUB_YEAR]->(y)"
                                  , {"i": i, "d": pub_data[i][1][0], "m": pub_data[i][1][1], "y": pub_data[i][1][2]})
                tx.commit()
                print "Done"
        
                #add all the author nodes
                print "Creating author nodes..."
                #determine the unique list of authors (no duplicates)
                autlist=[]
                for i in pub_data:
                        autlist+=pub_data[i][2]
                autlist = list(set(autlist))
        
                tx = graph.cypher.begin()
                for i in autlist:
                        tx.append("CREATE (a:Author {name:{i}})",{"i":i})
                tx.commit()
                print "Done"

        #create citations and authors connections
        print "Creating citations and authors connections, starting at : "
        tx = graph.cypher.begin()
        #determine where to start
        l = citation_no(pub_data,l)
        k=0
        for i in pub_data:
                if k<l:
                        k=k+1
                        continue
                print k
                for j in pub_data[i][0]:
                        tx.append("MATCH (p1:Publication {id:{id1}}),(p2:Publication {id:{id2}})  CREATE (p2)-[:CITED_BY]->(p1)",{"id1":i,"id2":j})            
                for j in pub_data[i][2]:
                        tx.append("MATCH (a:Author {name:{j}}),(p:Publication {id:{i}}) CREATE (a)-[:WROTE]->(p)",{"i":i,"j":j})
                k=k+1
                tx.commit()
                tx = graph.cypher.begin()
                
        print "Done"


#A function to only update the edges for the nodes included in the list
def selected_push_neo_graph(pub_data, l): #list of ids to update
        authenticate("localhost:7474","neo4j","azerty01")
        graph= Graph("http://localhost:7474/db/data")
        print "Creating selected citations and authors connections, starting at : "
        tx = graph.cypher.begin()

        for i in l:
                print i
                tx.append("MATCH ()-[r:CITED_BY]->(p1:Publication {id:{id}}) DELETE r",{"id":i})
                tx.commit()
                tx = graph.cypher.begin()
                tx.append("MATCH ()-[r:WROTE]->(p1:Publication {id:{id}}) DELETE r",{"id":i})
                tx.commit()
                tx = graph.cypher.begin()

                for j in pub_data[i][0]:
                        tx.append("MATCH (p1:Publication {id:{id1}}),(p2:Publication {id:{id2}})  CREATE (p2)-[:CITED_BY]->(p1)",{"id1":i,"id2":j})   
                        tx.commit()
                        tx = graph.cypher.begin()
                for j in pub_data[i][2]:
                        tx.append("MATCH (a:Author {name:{j}}),(p:Publication {id:{i}}) CREATE (a)-[:WROTE]->(p)",{"i":i,"j":j})
                tx.commit()
                tx = graph.cypher.begin()
                
        print "Done"


#Function to load the graph from the neo4j database
"""def pull_neo_graph(): 
        pub_data = {} #the data structure to hold the graph
        
        authenticate("localhost:7474", "neo4j", "azerty01") #put your login and password here
        graph = Graph("http://localhost:7474/db/data/")
        
        k=1
        for id in graph.cypher.execute("MATCH (n:Publication) RETURN n.id AS id"):
                print k
                if k == 104 :#for testing purposes
                        break
                pub_data[id.id]=[[],[],[]]
                for cits in graph.cypher.execute("MATCH (n:Publication {id:{id}}),(c:Publication) WHERE  (c)-[:CITED_BY]->(n) RETURN c.id AS id", {"id":id.id}):
                        pub_data[id.id][0].append(cits.id)
                for date in graph.cypher.execute("MATCH (n:Publication {id:{id}}),(d:Day),(m:Month),(y:Year) WHERE (n)-[:PUB_DAY]->(d) AND (n)-[:PUB_MONTH]->(m) AND (n)-[:PUB_YEAR]->(y) "
                                                 "RETURN d.name AS dn, m.name AS mn, y.name AS yn",{"id":id.id}):
                        pub_data[id.id][1] = [date.dn,date.mn,date.yn]
                for auts in graph.cypher.execute("MATCH (n:Publication {id:{id}}),(a:Author) WHERE (a)-[:WROTE]->(n) RETURN a.name AS name",{"id":id.id}):
                        pub_data[id.id][2].append(auts.name)
                k=k+1
        return pub_data"""

#Function to delete neo4j graph
def del_neo_graph():
        authenticate("localhost:7474", "neo4j", "azerty01") #put your login and password here
        graph = Graph("http://localhost:7474/db/data/")
        graph.delete_all()
        
#Verify that the graph has the correct number of edges, to ensure it has been loaded into neo4j properly
def verif_graph(pub_data):
        authenticate("localhost:7474","neo4j","azerty01")
        graph= Graph("http://localhost:7474/db/data")
        tx = graph.cypher.begin()
        l=0
        ll = []
        for i in pub_data:
                tx = graph.cypher.begin()
                record1 = graph.cypher.execute("MATCH ()-[r:CITED_BY]->(P:Publication {id:{i}}) RETURN count(r)={l} AS name",{"i":i,"l":len(pub_data[i][0])})
                record2 = graph.cypher.execute("MATCH ()-[r:WROTE]->(P:Publication {id:{i}}) RETURN count(r)={l} AS name",{"i":i,"l":len(pub_data[i][2])})
                
                if record1[0].name != True or record2[0].name != True: 
                        print i
                        ll+=[i]
                l=l+1
                if l%100 == 0: #just to keep track of the advancement
                        print l
        return ll



