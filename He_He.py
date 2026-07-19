import numpy as np
import math
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
from string import ascii_uppercase

# Implementation of He & He paper "Optimal Monotone Drawings of Trees", https://arxiv.org/pdf/1604.03921v1.pdf
# Instructions: Use PathDrawAlgorithm() function

def invert(tup):
	return tuple([tup[1],tup[0]])

def FareySeq(n, descending=False):
	"""Print the nth Farey sequence, either ascending or descending."""
	result=[]
	a, b, c, d = 0, 1, 1, n
	if descending: 
		a, c = 1, n-1
	result.append(tuple([a,b]))

	while (c <= n and not descending) or (a > 0 and descending):
		k = int((n + b) / d)
		a, b, c, d = c, d, (k*c-a), (k*d-b)
		result.append(tuple([a,b]))

	result.pop(0)	
	for i in range(len(result)):
		if invert(result[i]) not in result:
			result.append(invert(result[i]))		

	return sorted(result,key=lambda x: x[1]/float(x[0]))	


def FareySeq2(n, descending=False):
	"""Print the nth Farey sequence, either ascending or descending."""
	result=[]
	a, b, c, d = 0, 1, 1, n
	if descending: 
		a, c = 1, n-1
	result.append(tuple([a,b]))
	while (c <= n and not descending) or (a > 0 and descending):
		k = int((n + b) / d)
		a, b, c, d = c, d, (k*c-a), (k*d-b)
		result.append(tuple([a,b]))

	result.pop(0)	
	for i in range(len(result)):
		if invert(result[i]) not in result:
			result.append(invert(result[i]))		

	return [(1,0)]+sorted(result,key=lambda x: x[1]/float(x[0]))+[(0,1)]	


def PathDecomposition(tree):

	leafs=[]
	graph=dict()

	for line in nx.generate_adjlist(tree):
		ln=line.split(" ")
		graph.update({ ln[0]: ln[1:] })

	for v in graph.keys():
		if graph[v]==[]:
			leafs.append(v)

	b = list(nx.shortest_path(tree,leafs[0],'A'))
	T = set()
	T.update(b)
	bset = []
	leafs.pop(0)
	bset.append(b)
	for v in leafs:
		b = list(nx.shortest_path(tree,v,'A'))

		for u in b:
			if u in T:
				bset.append(list(nx.shortest_path(tree,v,u)))
				T.update(b)
				break;

	for i in range(0,len(bset)):
		del bset[i][-1]
				

	return bset

def LDPD(tree):
	'''Length Decreasing Path Decomposition of a Tree'''

	r = 'A'
	leafs = []
	graph = dict()
	b=[]
	T = set()
	for line in nx.generate_adjlist(tree):
		ln = line.split(" ")
		graph.update({ ln[0]: ln[1:] })

	for v in graph.keys():
		if graph[v] == []:
			leafs.append(v)
	leafs_paths = []

	for u in leafs:
		leafs_paths.append(list(nx.shortest_path(tree,u,r)))

	maxlen_path = max(leafs_paths,key=len)
	b.append(maxlen_path)
	T.update(maxlen_path)
	usedleafs = []
	usedleafs.append(b[0][0])

	while len(leafs) != len(usedleafs) :
		lens = []
		for v in leafs:
			if v not in usedleafs:
				b1 = list(nx.shortest_path(tree,v,'A'))
				for u in b1:
					if u in T:
						l=list(nx.shortest_path(tree,v,u))
						lens.append(l)
						break
		l = max(lens,key=len)
		b.append(l)
		T.update(l)
		usedleafs.append(l[0][0])

	return b


def c_partition(tree,ldpd,c):
	'''c-partition of an LDPD of a tree'''

	n = len(list(tree.nodes()))
	K = math.ceil(math.log(n,c))
	D = []
	D1 = []

	for b in ldpd:
		if len(b)-1 >= (n-1)/c and len(b)-1<=(n-1):
			D1.append(b)
	D.append(D1)

	for j in range(2,K+1):
		d = []
		for b in ldpd:
			if len(b)-1 >= (n-1)/(c**j) and len(b)-1<(n-1)/(c**(j-1)):
				d.append(b)
		D.append(d)		

	return D


def construct_prim_vectors(f,d,n):

	c = f+1
	K = math.ceil(math.log(n,c))

	def find_elements_between(source,start,end,howmany):
		begin = source.index(start)
		end = source.index(end)
		between = source[begin+1:end]
		result = between[0:howmany]
		return result

	P = FareySeq(d)

	P1 = P.index((1,1))
	S1 = P[0:f]
	S2 = P[P1+1:P1+f+1]
	R1 = [x for x in S1]+[(1,1)]+[y for y in S2]

	Pd2 = FareySeq2(d**2)
	R2 = []
	newR1 = [(1,0)]+R1+[(0,1)]

	for first,second in zip(newR1,newR1[1:]):
		elements_between=find_elements_between(Pd2,first,second,f)
		R2+=elements_between

	All_Rs = [R1,R2]
	union_of_R = [*R1 , *R2]
	union_of_R = sorted(union_of_R,key=lambda x: x[1]/float(x[0]))

	for j in range(3,K+1):
		if (1,0) not in union_of_R and (0,1) not in union_of_R:
			union_of_R.insert(0,(1,0))
			union_of_R.append((0,1))

		R = []
		jsource = FareySeq2(d**j)

		for first,second in zip(union_of_R , union_of_R[1:]):
			elements_between=find_elements_between(jsource,first,second,f)
			R+=elements_between

		union_of_R.remove((0,1))
		union_of_R.remove((1,0))
		union_of_R+=R
		union_of_R = sorted(union_of_R,key=lambda x: x[1]/float(x[0]))
		All_Rs.append(R)

	return [union_of_R,All_Rs]	


def getleafs(tree,start):

	leafs = []
	for node in tree.keys():
		if tree[node] == []:
			leafs.append(node)

	return leafs	


def PathDrawAlgorithm(g, root='A', display=False, save=False, filename="He_He.png", f=3, d=3, w_labels=False, n_size=25, n_color='black', l_color='black'):

	tree = g.copy()
	graph = dict()
	for line in nx.generate_adjlist(tree):
		ln = line.split(" ")
		graph.update({ ln[0]: ln[1:] })

	leafs = getleafs(graph,root)
	t = len(leafs)
	n = len(list(tree.nodes()))
	V = construct_prim_vectors(f,d,n)
	v = V[0]
	all_Rs = V[1]
	n = len(list(tree.nodes()))
	letters = list(ascii_uppercase) + [char1+char2 for char1 in ascii_uppercase for char2 in ascii_uppercase]
	integers = [x for x in range(0,n+1)]
	for_relabel = dict(zip(letters,integers))

	B = LDPD(tree)
	B = sorted(B,key=lambda x: for_relabel[x[0]] )
	D = c_partition(tree,B,f+1)

	R_levels = dict()
	i = 0
	for r in all_Rs:
		for rs in r:
			R_levels.update({tuple(rs):i})
		i+=1	


	for b in B:
		b.remove(b[-1])

	b_levels = dict()
	i = 0
	for b in D:
		for bs in b:
			b_levels.update({tuple(bs):i})
		i+=1	

	assignedEdges=dict()
	

	lastindex = 0
	for l in range(0,t):

		bl = B[l]
		bl_level = b_levels[tuple(bl)]
		for vector in v:
			if R_levels[vector] <= bl_level and v.index(vector)>lastindex:
				assignedEdges.update({ tuple(bl) : vector })
				lastindex=v.index(vector)
				break


	absCoords = dict()
	for path in B:
		for node in path:
			absCoords.update({ node:assignedEdges[tuple(path)] })
	
	gridPoints = dict()
	Root = root
	gridPoints.update({Root:(0,0)})
	ccwNodes = list(graph.keys())

	for node in ccwNodes:
		for child in graph[node]:
			if child not in gridPoints.keys():
				(x,y) = absCoords[child]
				(xp,yp) = gridPoints[node]
				gridPoints.update({child:(x+xp,y+yp)})
	
	data = pd.DataFrame.from_dict(gridPoints).to_numpy()
	x = data[0]
	y = data[1]
	maxx = max(x)
	maxy = max(y)
	n = len(list(graph.keys()))
	gridArea = n**2
	graphArea = maxx*maxy

  #------------------------------------------------------------------------------------------------------------------------
  #                                           Plot Graph usign matplotlib
  #------------------------------------------------------------------------------------------------------------------------

	if (display == True) or (save == True):
		plot = plt.figure(figsize=(7,5))
		nx.draw_networkx(tree, pos=gridPoints, with_labels=w_labels, node_size = n_size,arrows=False,node_color=n_color,font_color=l_color)
		plt.grid(color="gray")
		ax = plt.gca()
		ax.set_ylim([-1,maxy+1])
		ax.set_xlim([-1,maxx+1])
		ax.set_xticklabels([])
		ax.set_yticklabels([])
		plt.title('He & He Optimal Algorithm\n\nGrid Size: '+str(maxx+1)+' x '+str(maxy+1)+' ('+str(n)+' nodes)')
		ax.set_aspect(1)
		plt.xticks(np.arange(0,maxx+1,1)) 
		plt.yticks(np.arange(0,maxy+1,1)) 
		plot = plt.gcf()
		if save == True:
			plt.savefig(filename,bbox_inches='tight')
		if display == True:
			plt.show()

		plt.close()

	 
	return [graphArea,maxx,maxy,gridPoints] 

