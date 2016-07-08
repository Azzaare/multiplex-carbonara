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


def adjust_initials(aut):
        l = aut.split(".")
        ll = []
        lll = []
        for i in l:
                ll+= [i.lstrip().rstrip()]
        ll = ". ".join(ll)

        ll   = ll.split(" ")
        for i in ll:
                if len(i)==1: #if it's an initial
                        lll += [i+"."]
                else:
                        lll += [i]
        lll = " ".join(lll)
        return lll

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
	lp = [adjust_initials(x.lstrip().rstrip()).lower() for x in lp if x.lstrip().rstrip() != ""] #remove the spaces at the beginning and end of authors name, and add spaces between initials
        
        return lp

    
#Function for loading the data structure which associates for each publication the other publications which it cites, its publication date and its list of authors

#function to return list of unique authors
def author_list(pub_data):
        autlist = []
        for i in pub_data:
		autlist+=pub_data[i][2]
        autlist = list(set(autlist))
        return autlist

#function to return count of authors
def count_authors(pub_data):
        return len(author_list(pub_data))

#function which adjusts the initials to the correct format
def author_initials(name):
        tnames = name.lower().split(" ")
        tname = ""
        for s in tnames[:len(tnames)-1]:
                if s[len(s)-1]!='.':
                        tname += s[0]+'.'
                else:
                        tname+=s
        return tname+tnames[len(tnames)-1]

#function which checks if there are conflicts between different authors sharing the same initials
def check_author_initials_conflict(pub_data):
        autlist = author_list(pub_data)
        initial_table = {}
        for a in autlist:
                initial_table[author_initials(a)] = []
        for a in autlist:
                #if "".join(a.lower().split()) != author_initials(a):
                initial_table[author_initials(a)] += [a]

        #corrections
        #remove singletons
        to_delete = []
        for i in initial_table:
                if len(initial_table[i]) <= 1:
                        to_delete+=[i]
        for i in to_delete:
                del initial_table[i]
            
        k=0
        for i in initial_table:
                 print i,initial_table[i]
                 if len(initial_table[i])>2:
                         k+=1

        print k

#function to reduce the number of authors by fusioning authors according to whether one authors is just the initials of another author
def reduce_authors(pub_data): #PROBLEMATIC if the authors have the same initials especially if one of the authors only appears with his initials and the other authors has both initials and full name
        #First get lists of all authors, then classify authors by initials. If two (and only two) authors share the same initials, and if one of them is equal to the initials, then mark the change to use the other author name

        #######BUGGGGGGG with jr.

        
        autlist = author_list(pub_data)
        initial_table = {}
        change_table = {}
        for a in autlist: #build initials tables
                initial_table[author_initials(a)] = []
        for a in autlist:
                initial_table[author_initials(a)] += [a]

        #if one author corresponds to one initial, nothing to do. If two authors correspond to one initial check if we can reduce. If 3 or more authors correspond to the same initial too complicated to do anything
        for i in initial_table:
                if len(initial_table[i]) == 2:
                        if "".join(initial_table[i][0].lower().split()) == author_initials(initial_table[i][0]):
                                change_table[initial_table[i][0]] = initial_table[i][1]
                        elif "".join(initial_table[i][1].lower().split()) == author_initials(initial_table[i][1]):
                                change_table[initial_table[i][1]] = initial_table[i][0]
        #now we reduce
        for id in pub_data:
                for i in range(len(pub_data[id][2])):
                        if pub_data[id][2][i] in change_table:
                                pub_data[id][2][i] = change_table[pub_data[id][2][i]]
        

#Function which loads the data into the data structure
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
        reduce_authors(pub_data) #reduce the number of authors (check if some different authors are the same author but with name written differently 
        print "Data loaded"
        return pub_data



