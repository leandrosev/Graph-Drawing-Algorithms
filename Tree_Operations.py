#------------------------------------------------------------------------------------------------------------------------
# @author: Evangelidakis Leandros
# @School of Applied Mathematics and Physical Sciences 
# -----------------------------------------------------------------------------------------------------------------------

from itertools import permutations
import networkx as nx
import random
from string import ascii_uppercase
import matplotlib.pyplot as plt


def isvalid(parenth):

	''' Checks if a string of parenthesis is balanced '''

	left = right = 0
	if parenth[0] == ')':
		return False
	for char in parenth:
		if char == '(':
			left += 1
		elif char == ')':
			right += 1
		if right > left:
			return False
	if right == left:	
		return True
	else:
		return False	

#print(isvalid( '()()()())(') )
#Result : False
#------------------------------------------------------------------------------------------------------------------------

def genParenLazy(n):
	''' 
	Creates every valid combination of balanced parenthesis
	input: number of nodes minus one for the fixed root	
	Complexity: Brute force method, takes O(2^(2n)) to create all combinations and O(n) to check them so O(n*2*(2n)) time and total space
	''' 
	if n<0:
			err = 'Input must be a positive integer'
			return err			
	perms = set(''.join(p) for p in permutations('(' * n + ')' * n))
	return [s for s in perms if isvalid(s)]
  

#print(genParenLazy(3))  
#Result : ['(())()', '()()()', '((()))', '(()())', '()(())']
#------------------------------------------------------------------------------------------------------------------------

def genParenFast(n):
	''' 
	Creates every valid combination of balanced parenthesis
	input: number of nodes minus one for the fixed root	
	In every iteration we check if the next action results into a valid set or not.
	We know that a parenthesis string is valid if in begins with an open parenthesis and at any point, the number of open parenthesis 
	is the same or larger than the number of closed parenthesis.

	Complexity: O(4^n / n^1.5 ) (See Catalan numbers for better understanding)

	'''
	n = n-1

	if n <= 0:
		err = 'Input must be a positive integer'
		return err

	res = []
	def process(string = '', open = 0, close = 0):
		if len(string) == 2 * n:
			res.append(string)
			return
		if open < n:
			process(string+'(', open+1, close)
		if close < open:
			process(string+')', open, close+1)

	process()
	return res

#print(genParenFast(4))
#Result: ['((()))', '(()())', '(())()', '()(())', '()()()']

#------------------------------------------------------------------------------------------------------------------------

def random_tree(n):
# Creates a random tree with n nodes considering a fixed root
# First create a random shuffle of parenthesis and check each one until a valid one  is found
#  Input:              n : integer (n+1 nodes of the tree)
#  Output: parenthesis representation of the tree	
# Complexity: O(n^1.5)

	i = 0
	def random_gen(n):
		chars = ['(', ')'] * (n-1)
		random.shuffle( chars )
		return ''.join(chars)

	found = False
	while ( found == False ):
		parenth = random_gen(n)
		if ( isvalid(parenth) == True ) :
			#print('Number of attempts: ',i)
			return parenth
		else:
			i = i+1
			continue	

#print(random_tree(4))

#------------------------------------------------------------------------------------------------------------------------
def paren_to_nxgraph(parenthesis):
# Transforms a string of balanced parenthesis into a networkx graph data type
# For node labels, letters of english alphabet is used
# Input : string 
# Output: networkx Graph

	letters = list(ascii_uppercase) + [char1+char2 for char1 in ascii_uppercase for char2 in ascii_uppercase]

	if  isvalid(parenthesis) == False:
		message = 'Not a valid string of parenthesis'
		return message
	else:	
		g = nx.Graph()
		parent = letters.pop(0)
		child = 0
		index = 0
		g.add_node('A')
		nodes = []
		nodes.append('A')
		open = close = 0
		parentgraph = dict()
		parentgraph.update({'A':'A'})
		for paren in parenthesis:
			if paren == '(':
				index+=1
				child = letters.pop(0)
				g.add_node(child)
				nodes.append(child)
				g.add_edge(parent,child)
				parentgraph.update({child:parent})
				parent = nodes[index]
				open+=1
			else:
				parent = parentgraph[child]
				child  = parentgraph[child]
				close+=1
			if open == close:
				open   = 0
				close  = 0
				parent = 'A'
		return g

#------------------------------------------------------------------------------------------------------------------------

#paren_to_nxgraph("()()()")


#------------------------------------------------------------------------------------------------------------------------
#                                           Plot NetworkX Graph 
#------------------------------------------------------------------------------------------------------------------------
def draw(g):
	'''Simple graph visualization using NetworkX drawing functions'''	
	plt.figure(figsize=(7,7))
	nx.draw_spectral(g, node_size = 20)
	plt.title('NetworkX drawing functions')
	plt.show()

#g = paren_to_nxgraph("(()()())(()()())(()()())(()()())(()()())(()()())(()()())(()()())(()()())")
#draw(g)

#------------------------------------------------------------------------------------------------------------------------
#                                           Plot NetworkX Graph usign Monotone Drawing Algorithms
#------------------------------------------------------------------------------------------------------------------------

#g = paren_to_nxgraph("(()()())(()()())(()()())(()()())(()()())(()()())(()()())(()()())(()()())")
#ang.getGridArea(g,root_node='A',display=True,save=False,algo="bfs",view=None,filename="tree1.png",w_labels=False,n_size=25,n_color='black',l_color='black')

#------------------------------------------------------------------------------------------------------------------------
